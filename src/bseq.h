
/**
 * @file bseq.h
 * @brief SIMD-parallel FASTA/Q parser (<- minialign.c <- bseq.h in minimap but everything has changed)
 */

/* include all */
#include "utils/utils.h"


/* FASTA/Q parser configurations */
typedef struct {
	size_t batch_size;						/* buffer (block) size */
	size_t head_margin;						/* head space before bseq_batch_t object */
	uint64_t keep_qual, keep_comment;		/* 1 to keep quality string */
	uint8_t conv_table[16];					/* can be '\0' */
} bseq_conf_t;

/* packed record container */
typedef struct {
	uint32_t nlen, hlen;					/* name and header (name + comment) length */
	uint8_t *seq;							/* 4-bit packed seq, raw qual string, tags (comment) */
	size_t slen;							/* sequence length */

	/* for index construction and multithreaded alignment computation */
	union {
		void *ptr;
		uint64_t flag;
	} u;
} bseq_meta_t;
_static_assert(sizeof(bseq_meta_t) == 32);


/* getters */
#define bseq_seq(_x)				( (uint8_t *)(_x)->seq )
#define bseq_seq_len(_x)			( (size_t)(_x)->slen )

#define bseq_is_qual_set(_x)		( bseq_qual(_x)[0] != '\0' )
#define bseq_qual(_x)				( (uint8_t *)&(_x)->seq[(_x)->slen + 2] )
#define bseq_qual_len(_x)			( bseq_is_qual_set(_x) ? (size_t)(_x)->slen : 0 )

#define bseq_name(_x)				( (char const *)&(_x)->seq[-((size_t)(_x)->hlen)] )
#define bseq_name_len(_x)			( (size_t)(_x)->nlen )
#define bseq_comment(_x)			( (char const *)&(_x)->seq[-((size_t)(_x)->hlen) + (_x)->nlen + 1] )
#define bseq_comment_len(_x)		( (size_t)(_x)->hlen - (size_t)(_x)->nlen - 2 )

#define bseq_circ(_x)				( (_x)->u.flag )
#define bseq_aln(_x)				( (_x)->u.ptr )

/* memory management */

#define BSEQ_BATCH_SIZE				( 2ULL * 1024 * 1024 )
#define BSEQ_MGN					( 64 )			/* buffer margin length */
_static_assert((BSEQ_MGN % 32) == 0);

typedef struct {
	void *base;
	size_t size;
} bseq_bin_t;

/* sequence container */
typedef struct {
	/* block information */
	bseq_bin_t bin;
	uint32_t head_margin, base_id;

	/* sequences */
	size_t cnt;					/* #sequences */
	bseq_meta_t meta[];			/* array of bseq_meta_t */
} bseq_batch_t;
_static_assert(sizeof(bseq_batch_t) == sizeof(bseq_meta_t));
typedef kvecm_t(bseq_meta_t) bseq_metam_v;			/* margined */

#define bseq_bin(_x)				( (_x)->bin )
#define bseq_bin_ptr(_x)			( _add_offset((_x)->base, BSEQ_MGN) )
#define bseq_bin_size(_x)			( (_x)->size - 2 * BSEQ_MGN )

#define bseq_bin_base_ptr(_x)		( (_x)->base )
#define bseq_bin_base_size(_x)		( (_x)->size )

#define bseq_meta_ptr(_x)			( (_x)->meta )
#define bseq_meta_cnt(_x)			( (_x)->cnt )

typedef struct {
	size_t cnt;					/* #sequences read */
	uint64_t status;			/* non-zero if error occurred */
} bseq_close_t;


typedef struct {
	/* conversion table (constant) */
	uint8_t conv[32];

	/* internal buffer */
	uint8_t *buf;
	size_t batch_size;

	/* fetcher */
	uint8_t *p, *t;					/* base, current pointer, and tail pointer */
	size_t acc;						/* accumulator to parse qual string */
	size_t scnt;					/* #sequences read, minimum sequence length (shorter is discarded) */
	uint8_t state, is_eof;

	/* constants */
	uint8_t delim[2], keep_qual, keep_comment;
	uint32_t head_margin;

	/* decompresser context */
	rbread_t rb;
} bseq_file_t;

#define bseq_is_eof(_fp)			( (_fp)->is_eof > 1 )
#define bseq_is_broken(_fp)			( (_fp)->is_eof > 2 )

static _force_inline
void bseq_buf_init(bseq_file_t *fp, size_t batch_size)
{
	size_t const malloc_size = _roundup(sizeof(uint8_t) * batch_size + 2 * BSEQ_MGN, ARCH_HUGEPAGE_SIZE);
	uint8_t *ptr = aligned_malloc(malloc_size);
	_memset_blk_u(ptr, '\n', BSEQ_MGN);		/* clear head margin so that parsing will not depend on uninitialized values */

	/* save */
	fp->batch_size = batch_size;
	fp->buf = ptr + BSEQ_MGN;

	/* make sure p >= t for the initial state */
	fp->p = fp->buf;
	fp->t = NULL;
	return;
}

static _force_inline
void bseq_buf_destroy(bseq_file_t *fp)
{
	if(fp->buf != NULL) { free(fp->buf - BSEQ_MGN); }
	return;
}

static _force_inline
size_t bseq_fetch_block(bseq_file_t *fp)
{
	/* read a block */
	uint64_t const back = fp->p == fp->t;		/* '\n' remains at the last byte of the buffer */
	size_t const bytes = rbread_bulk(&fp->rb, fp->buf, fp->batch_size);		/* always bypass internal buffer; up to batch_size */

	/* save pointers */
	fp->p = fp->buf - back;
	fp->t = fp->buf + bytes;

	/* add padding */
	_memset_blk_u(fp->t, '\n', BSEQ_MGN);		/* fill margin of the input buffer */

	/* update EOF */
	fp->is_eof = MAX2(fp->is_eof, rbeof(&fp->rb));

	debug("bytes, state(%u), len(%zu), is_eof(%u), rb(%u, %zu, %zu), next(%x, %x, %x, %x)", fp->state, bytes, fp->is_eof, fp->rb.eof, fp->rb.head, fp->rb.tail, fp->p[0], fp->p[1], fp->p[2], fp->p[3]);
	return(bytes);
}


/* bseq_read_fasta components */
#define _match(_v1, _v2)		( ((v32_masku_t){ .mask = _mask_v32i8(_eq_v32i8(_v1, _v2)) }).all )

typedef struct bseq_work_s {
	v32i8_t cvl, cvh, dv, sv, tv, lv;
	bseq_meta_t *s;
	uint8_t *p, *t, *q;
	size_t n;
} bseq_work_t;

static _force_inline
void bseq_init_work(bseq_work_t *w, bseq_file_t *fp, bseq_metam_v const *seq, uint8m_v const *bin)
{
	w->cvl = _from_v16i8_v32i8(_loadu_v16i8(&fp->conv[0]));
	w->cvh = _from_v16i8_v32i8(_loadu_v16i8(&fp->conv[16]));

	w->dv = _set_v32i8(fp->delim[1]);
	w->sv = _set_v32i8(' ');
	w->tv = _set_v32i8('\t');
	w->lv = _set_v32i8('\n');

	w->s = &seq->a[seq->n - 1];					/* restore previous states */
	w->p = fp->p;
	w->t = fp->t;

	/* load pointers from output memory block */
	w->q = bin->a;
	w->n = bin->n;
	return;
}

static _force_inline
uint64_t bseq_finish_work(bseq_work_t *w, bseq_file_t *fp, uint8m_v *bin)
{
	debug("return, p(%p), t(%p)", w->p, w->t);
	fp->p = w->p;

	/* write back pointers to the output memory block */
	bin->n = w->n;
	debugblock({ if(bin->n > bin->m) { debug("buffer overrun, bin->n(%lu), bin->m(%lu)", bin->n, bin->m); trap(); } });
	return(w->p >= w->t);
}

static _force_inline
size_t bseq_seq_transvp(bseq_work_t *w, v32i8_t v, size_t l) {
	/* apply lower and upper tables */
	v32i8_t const nl = _shuf_v32i8(w->cvl, v);
	v32i8_t const nh = _shuf_v32i8(w->cvh, v);

	/* merge two */
	v32i8_t const selv = _xor_v32i8(_add_v32i8(v, v), _shl_v32i8(v, 3));
	v32i8_t const nucl = _sel_v32i8(selv, nh, nl);

	_storeu_v32i8(&w->q[w->n], nucl);
	return(l);
}
/* without conversion */
static _force_inline
size_t bseq_idv(bseq_work_t *w, v32i8_t v, size_t l) {
	_storeu_v32i8(&w->q[w->n], v); 
	return(l);
}
/* escape tabs */
static _force_inline
size_t bseq_escv(bseq_work_t *w, v32i8_t v, size_t l) {
	_storeu_v32i8(&w->q[w->n], _sel_v32i8(v, w->sv, _eq_v32i8(v, w->tv)));
	return(l);
}
/* skip (nothing to do) */
static _force_inline
size_t bseq_skipv(bseq_work_t *w, v32i8_t v, size_t l) {
	_unused(w); _unused(v); _unused(l);
	return(0);
}

static _force_inline
size_t bseq_readline(
	bseq_work_t *w,
	size_t (*const vector)(bseq_work_t *w, v32i8_t v, size_t l))
{
	uint8_t const *b = w->p;
	debug("vector(%p), p(%p), t(%p), n(%lu), c(%x)", vector, w->p, w->t, w->n, w->q[w->n]);

	/* bulk parse loop */
	size_t len;
	do {
		v32i8_t r = _loadu_v32i8(w->p);
		uint64_t m = _match(r, w->lv);
		ZCNT_RESULT size_t l = _tzc_u32(m);	/* ZCNT_RESULT is gcc-4.7 workaround */
		len = l;							/* we do not need to consider the remaining input length, (size_t)(w->t - w->p), because the tail margin of the buffer is filled with '\n'. */
		w->n += vector(w, r, len);
		w->p += 32;
	} while(len >= 32);
	w->p += len - 32;

	/* remove '\r' for '\r\n' */
	size_t const back  = w->p[-1] == '\r';
	size_t const bytes = (w->p - b) - back;
	w->n -= (vector != bseq_skipv) & back;

	/* forward pointer (skip an obvious delimiter at the end) */
	w->p++;

	/* fixup last byte for the next readline */
	w->q[w->n] = '\0';				/* mark empty; invariant condition */
	debug("p(%p), t(%p), n(%lu), next(%x, %x, %x, %x)", w->p, w->t, w->n, w->p[-1], w->p[0], w->p[1], w->p[2]);
	return(bytes);
}

static _force_inline
void bseq_fixup_name(bseq_work_t *w, bseq_file_t *fp)
{
	size_t nlen = fp->acc, p = w->n - nlen, t = w->n;
	v32i8_t const sv = _set_v32i8(' ');

	debug("nlen(%zu), p(%zu), t(%zu)", nlen, p, t);

	/* strip head spaces */ {
		uint64_t const m = _match(sv, _loadu_v32i8(&w->q[p]));
		ZCNT_RESULT size_t l = _tzc_u32(~m);
		size_t const len = l; p += len; nlen -= len;
		debug("nlen(%zu), len(%zu), p(%zu)", nlen, len, p);
	}

	/* strip tail spaces */ {
		uint64_t const m = _match(sv, _loadu_v32i8(&w->q[t - 32]));
		ZCNT_RESULT size_t l = _lzc_u32(~m);
		size_t const len = l; t -= len; nlen -= len;
		debug("nlen(%zu), len(%zu), t(%zu)", nlen, len, t);
	}

	/* split at the first space */
	size_t q = p;
	while(q < t) {
		uint64_t const m = _match(sv, _loadu_v32i8(&w->q[q]));
		size_t const len = MIN2(32, t - q);
		q += len;
		if(m) {
			ZCNT_RESULT size_t u = _tzc_u32(m);
			size_t const v = u, x = MIN2(v, nlen);	/* clip split position at the tail of name-comment section */
			q += x - len;
			debug("p(%zu), q(%zu), t(%zu), len(%zu)", p, q, t, x);
			break;
		}
	}

	/* remove comment if needed */
	if(fp->keep_comment == 0) { t = q + 1; }
	t += t == q;					/* add space if keep_comment == 1 and comment is missing */

	/* fill terminators */
	w->s->nlen = q - p;
	w->s->hlen = (t + 1) - p;		/* might have further spaces at the head of comment */
	debug("nlen(%u), hlen(%u), p(%zu), q(%zu), t(%zu), keep(%u)", w->s->nlen, w->s->hlen, p, q, t, fp->keep_comment);

	w->q[q] = '\0';
	w->q[t] = '\0';					/* string terminator; keep q[n] initialized (invariant condition) */
	w->n = t;

	debug("name(%lu, %s), comment(%lu, %s), p(%zu), q(%zu), t(%zu), n(%zu)",
		bseq_name_len(w->s), &w->q[w->n + 1 - w->s->hlen],
		bseq_comment_len(w->s), &w->q[w->n + 1 - w->s->hlen + w->s->nlen + 1],
		p, q, t, w->n
	);
	return;
}

static _force_inline
void bseq_fixup_seq(bseq_work_t *w, bseq_file_t *fp)
{
	w->s->seq   = (uint8_t *)w->n - fp->acc;
	w->s->slen  = fp->acc;
	w->s->u.ptr = NULL;
	debug("done, nlen(%u), hlen(%u), slen(%u), seq(%p)", w->s->nlen, w->s->hlen, w->s->slen, w->s->seq);
	return;
}

#undef _match

#define BSEQ_STATE_HEADER			( 0x01 )
#define BSEQ_STATE_BODY				( 0x02 )
#define BSEQ_STATE_SEQ				( 0x04 )
#define BSEQ_STATE_QUAL				( 0x08 )
#define BSEQ_STATE_SKIP				( 0x10 )

/**
 * @fn bseq_read_fasta
 * @brief parse one sequence, returns 0 for success, 1 for buffer starvation, 2 for broken format
 * bin must have enough space (e.g. 2 * buffer)
 */
static _force_inline
uint64_t bseq_read_fasta(bseq_file_t *fp, bseq_metam_v *seq, uint8m_v *bin)
{
	#define _state(_n, _cond, _body, _post) { \
		debug("forward state(%u), n(%zu)", _n, w.n); \
		fp->state = _n; \
		fp->acc = 0; \
		case _n: /* debug("enter state(%u)", _n); */ \
		while(1) { \
			_body; \
			uint64_t const flag = (_cond); \
			if(w.p >= w.t) { goto _refill; } \
			if(!flag) { break; } \
		} \
		_post; \
		w.q[++w.n] = 0;		/* invariant condition; open new empty byte for the next parsing */ \
	}

	/* error if buffer is empty */
	if(fp->p >= fp->t) {
		return(1);
	}

	/* keep them on registers */
	bseq_work_t w __attribute__(( aligned(64) ));
	bseq_init_work(&w, fp, seq, bin);

	debug("enter, state(%u), p(%p), t(%p), eof(%u)", fp->state, w.p, w.t, fp->is_eof);
	switch(fp->state) {							/* dispatcher for the first iteration */
		default:								/* idle or broken */
		if(*w.p++ != fp->delim[0]) { break; }	/* '>' */
		w.s = kvm_pushp(bseq_meta_t, *seq);
		_state(BSEQ_STATE_HEADER | BSEQ_STATE_SEQ, ({ 0; }),	/* sequence name */
			{ fp->acc += bseq_readline(&w, bseq_idv); },
			{ bseq_fixup_name(&w, fp); }
		);
		_state(BSEQ_STATE_BODY | BSEQ_STATE_SEQ,/* sequence */
			({ (*w.p != fp->delim[1]); }),
			{ fp->acc += bseq_readline(&w, bseq_seq_transvp); },
			{ bseq_fixup_seq(&w, fp); }
		);
		if(fp->delim[0] != '>') {				/* parse qual if fastq */
			_state(BSEQ_STATE_HEADER | BSEQ_STATE_QUAL, ({ 0; }),
				{ bseq_readline(&w, bseq_skipv); }, {}
			);
			if(fp->keep_qual) {					/* save if keep_qual */
				_state(BSEQ_STATE_BODY | BSEQ_STATE_QUAL,
					({ fp->acc < w.s->slen; }),
					{ fp->acc += bseq_readline(&w, bseq_idv); }, {}
				);
			} else {
				_state(BSEQ_STATE_BODY | BSEQ_STATE_QUAL | BSEQ_STATE_SKIP,
					({ fp->acc < w.s->slen; }),
					{ fp->acc += bseq_readline(&w, bseq_skipv); }, {}
				);
			}
			// w.q[++w.n] = 0;		/* for qual name */
		} else {
			/* we need two spaces for dummy qual string */
			w.q[++w.n] = 0;		/* for qual name */
			w.q[++w.n] = 0;		/* for qual body */
		}
		fp->state = 0;
	}
_refill:;
	debug("break, state(%u), eof(%u), p(%p), t(%p), n(%zu)", fp->state, fp->is_eof, fp->p, fp->t, w.n);
	if(fp->is_eof == 1 && (fp->state & BSEQ_STATE_BODY)) {
		if(fp->state & BSEQ_STATE_SEQ) { bseq_fixup_seq(&w, fp); }
		fp->state = 0;					/* reset state so that it is not regarded as an error */
		fp->is_eof = 2;					/* already reached the end in the previous fetch; mark as all done */
	}
	return(bseq_finish_work(&w, fp, bin));

	#undef _state
	#undef _init
	#undef _term
}

static _force_inline
uint64_t bseq_append(bseq_file_t *fp, bseq_metam_v *seq, uint8m_v *bin)
{
	debugblock({ if(fp->is_eof > 1) { trap(); } });
	while(bin->n < fp->batch_size) {				/* fetch-and-parse loop */
		while(bseq_read_fasta(fp, seq, bin)) {		/* buffer starved and sequence block continues */
			if(bseq_fetch_block(fp) == 0) {			/* fetch next */
				fp->is_eof = 2;
				break;
			}
			
			/* reserve room for the next parsing unit */
			kvm_reserve(uint8_t, *bin, kvm_cnt(*bin) + 2 * fp->batch_size);
			debug("continue, state(%u), is_eof(%u), seq(%p, %lu), bin(%p, %lu)", fp->state, fp->is_eof, kvm_ptr(*seq), kvm_cnt(*seq), kvm_ptr(*bin), kvm_cnt(*bin));
		}
		debug("finished seq, state(%u), is_eof(%u), seq(%p, %lu), bin(%p, %lu)", fp->state, fp->is_eof, kvm_ptr(*seq), kvm_cnt(*seq), kvm_ptr(*bin), kvm_cnt(*bin));

		if(_unlikely(fp->is_eof == 2)) { break; }	/* reached end of file */
		if(_likely(fp->state == 0)) { continue; }

		/* an error occurred */
		fp->is_eof = 3;								/* mark error occurred */
		if(_likely(kvm_cnt(*seq) > 0)) {
			(void)kvm_pop(*seq);					/* remove last (broken) element */
		}
		if(_unlikely(kvm_cnt(*seq) == 0)) {
			fp->is_eof = MAX2(fp->is_eof, 2);		/* no sequence is remaining */
			return(1);
		}
		break;
	}
	return(0);
}

static _force_inline
void bseq_init_bin(bseq_file_t const *fp, bseq_metam_v *seq, uint8m_v *bin)
{
	/* reserve bseq_batch_t space at the head of bseq_meta_t array */
	kvm_init(*seq, sizeof(bseq_batch_t) + fp->head_margin, 0);
	// debug("h(%u), margin(%u, %u)", sizeof(bseq_batch_t) + fp->head_margin, seq->h, seq->t);

	/* margined sequence bin */
	kvm_init(*bin, BSEQ_MGN, BSEQ_MGN);
	kvm_reserve(uint8_t, *bin, 2 * fp->batch_size);
	kvm_fill_margin(uint8_t, *bin, 0);

	/* clear the entire memory (for debugging) */
	debugblock({ memset(kvm_ptr(*bin), 0, 2 * fp->batch_size); });

	/* invariant condition: q[n] is always initialized */
	uint8_t *q = kvm_ptr(*bin);
	*q = 0;
	return;
}

static _force_inline
bseq_bin_t bseq_finalize_bin(bseq_file_t *fp, bseq_metam_v *seq, uint8m_v *bin)
{
	/* shrink memory if vacancy is large */
	if(kvm_max(*bin) - kvm_cnt(*bin) > 2 * fp->batch_size) {
		kvm_resize(uint8_t, *bin, kvm_cnt(*bin));
	}

	/* fill zero margin */
	kvm_fill_margin(uint8_t, *bin, 0);

	/* add offsets to validate pointers */
	ptrdiff_t const b = (ptrdiff_t)kvm_ptr(*bin);
	kvm_foreach(bseq_meta_t, *seq, { p->seq += b; });

	// debug("state(%u), is_eof(%u), cnt(%zu)", fp->state, fp->is_eof, kvm_cnt(*seq));

	return((bseq_bin_t){
		.base = kvm_base_ptr(*bin),
		.size = kvm_base_size(*bin)
	});
}

static _force_inline
bseq_batch_t *bseq_read(bseq_file_t *fp)
{
	if(fp == NULL || fp->is_eof > 1) { return(NULL); }	/* check the context is valid */

	/* allocate memory */
	bseq_metam_v seq = { 0 };
	uint8m_v buf = { 0 };
	bseq_init_bin(fp, &seq, &buf);	

	/* fetch block */
	if(bseq_append(fp, &seq, &buf) != 0) { goto _bseq_read_error; }
	size_t const base_id = fp->scnt;
	fp->scnt += kvm_cnt(seq);				/* update sequence counter */

	/* shrink memory and adjust pointers */
	bseq_bin_t bin = bseq_finalize_bin(fp, &seq, &buf);

	/* store meta to the bseq_batch_t object, allocated at the head of the bseq_meta_t array */
	bseq_batch_t *s = (bseq_batch_t *)kvm_ptr(seq);
	// debug("s(%p, %p)", kvm_ptr(seq), kvm_base_ptr(seq));

	s[-1] = (bseq_batch_t){
		.bin = bin,
		.head_margin = fp->head_margin,
		.base_id = base_id,
		.cnt = kvm_cnt(seq)
	};
	return(&s[-1]);

_bseq_read_error:;
	kvm_destroy(seq);
	kvm_destroy(buf);
	return(NULL);
}

static _force_inline
void *bseq_ptr(bseq_batch_t *batch)
{
	return(batch->bin.base);
}

static _force_inline
void bseq_bin_free(bseq_batch_t *batch) {
	free(batch->bin.base);
}
static _force_inline
void bseq_meta_free(bseq_batch_t *batch) {
	if(batch != NULL) { free(_sub_offset(batch, batch->head_margin)); }
}
static _force_inline
void bseq_free(bseq_batch_t *batch)
{
	bseq_bin_free(batch);
	bseq_meta_free(batch);
	return;
}

static _force_inline
bseq_close_t bseq_close(bseq_file_t *fp)
{
	if(fp == NULL) { return((bseq_close_t){ .cnt = 0, .status = 0 }); }

	size_t const scnt     = fp->scnt;
	uint64_t const status = bseq_is_broken(fp);

	rbclose_static(&fp->rb);
	bseq_buf_destroy(fp);
	free(fp);
	return((bseq_close_t){
		.cnt    = scnt,
		.status = status
	});
}

static _force_inline
void bseq_load_conv(uint8_t *q, uint8_t const *p)
{
	#define _c(x)		( ((x) & 0x1f) ^ (((x) & 0x40)>>2) )

	/* ASCII -> 4bit encoding conversion table */
	static uint8_t const conv[32] __attribute__(( aligned(32) )) = {
		/* '\0' and ' ' (as EOF) */
		[0] = 0xff,

		/* simple nucleotides */
		[_c('A')] = A,
		[_c('C')] = C,
		[_c('G')] = G,
		[_c('T')] = T,
		[_c('U')] = T,
		[_c('N')] = N,

		/* ambiguous bases */
		[_c('R')] = R,
		[_c('Y')] = Y,
		[_c('S')] = S,
		[_c('W')] = W,
		[_c('K')] = K,
		[_c('M')] = M,

		[_c('B')] = B,
		[_c('D')] = D,
		[_c('H')] = H,
		[_c('V')] = V,

		/* end-of-line marker */
		[_c('\r')] = 0xff,
		[_c('\n')] = 0xff
	};

	/* copy default conversion table to working buffer */
	uint8_t w[32] __attribute__(( aligned(32) )) = { 0 };
	_store_v32i8(w, _load_v32i8(conv));

	/* use default table if provided table is all zero */
	v16i8_t const cv = _loadu_v16i8(p), eqv = _eq_v16i8(cv, _zero_v16i8());
	if(((v16_masku_t){ .mask = _mask_v16i8(eqv) }).all != 0xffff) {

		/* 0x00 and 0x0f are excluded */
		for(size_t i = 0x00; i < 0x0f; i++) {
			uint8_t const ch = "NACMGRSVTWYHKDBN"[i];

			/* overwrite */
			w[_c(ch)] = p[i];
		}
	}

	/* save */
	_storeu_v32i8(q, _load_v32i8(w));
	return;

	#undef _c
}

static _force_inline
void bseq_load_conf(bseq_file_t *fp, bseq_conf_t const *conf)
{
	bseq_conf_t const def = { 0 };
	bseq_conf_t const *c = conf == NULL ? &def : conf;
	size_t const batch_size = c->batch_size == 0 ? BSEQ_BATCH_SIZE : c->batch_size;

	/* copy flags */
	*fp = (bseq_file_t){
		.batch_size   = batch_size,
		.keep_qual    = c->keep_qual,
		.keep_comment = c->keep_comment,
		.head_margin  = _roundup(c->head_margin, sizeof(bseq_batch_t))
	};

	/* build conversion table */
	bseq_load_conv(fp->conv, c->conv_table);
	return;
}

static _force_inline
bseq_file_t *bseq_open(bseq_conf_t const *conf, char const *fn)	/* file path; "-" to use stdin */
{
	/* exclude trivial failure */
	if(fn == NULL || fn[0] == '\0') {
		error("empty file path `'.");
		goto _bseq_open_fail;
	}

	/* create instance */
	bseq_file_t *fp = (bseq_file_t *)calloc(1, sizeof(bseq_file_t));
	bseq_load_conf(fp, conf);

	/* open stream */
	if(rbopen_bulk_static(&fp->rb, fn, fp->batch_size) != 0) {
		error("failed to open file `%s'.", fn);
		goto _bseq_open_fail;
	}

	/* init buffer; buffer size shoule be the same to the rb internal; fetch the first block */
	bseq_buf_init(fp, batch_size);
	if(bseq_fetch_block(fp) == 0) {
		fp->is_eof = 2;		/* empty file; we can't read any sequence but it is not an error on opening file */
		return(fp);
	}

	/* file has content; determine format. if format was not determined it is regarded as an error */
	v32i8_t const  v = _loadu_v32i8(fp->p), da = _set_v32i8('>'), dq = _set_v32i8('@');
	uint64_t const mask = ((v32_masku_t){
		.mask = _mask_v32i8(_or_v32i8(_eq_v32i8(v, da), _eq_v32i8(v, dq)))
	}).all;

	if(mask == 0) {
		error("failed to determine file format for `%s'.", fn);
		goto _bseq_open_fail;
	}
	ZCNT_RESULT uint64_t fwd = _tzc_u32(mask);
	fp->p += fwd;
	fp->delim[0] = fp->p[fwd];
	fp->delim[1] = fp->p[fwd] == '@' ? '+' : '>';
	return(fp);

_bseq_open_fail:;
	bseq_close(fp);
	return(NULL);
}


unittest( .name = "bseq.fasta" ) {
	char const *filename = "./bseq.unittest.tmp";
	char const *content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\r\n\r\n"
		">  test2\n\nAAAA\n"
		">test3 comment comment  \nACGT\n\n";

	FILE *fp = fopen(filename, "w");
	ut_assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	for(size_t i = 0; i < 4; i++) {
		bseq_conf_t const conf = {
			.batch_size   = BSEQ_BATCH_SIZE,
			.head_margin  = 64,
			.keep_qual    = (i & 0x01) != 0,
			.keep_comment = (i & 0x02) != 0,
			.conv_table = {
				[A] = 0x00,
				[C] = 0x01,
				[G] = 0x02,
				[T] = 0x03
			}
		};
		bseq_file_t *b = bseq_open(&conf, filename);
		ut_assert(b != NULL);

		bseq_batch_t *s = bseq_read(b);

		ut_assert(s != NULL);
		ut_assert(s->cnt == 4, "cnt(%u)", s->cnt);
		ut_assert(s->bin.base != NULL);
		ut_assert(s->bin.size > 0, "size(%lu)", s->bin.size);

		bseq_meta_t const *c = bseq_meta_ptr(s);
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu), hlen(%zu), name(%s)", bseq_name_len(c), c[0].hlen, bseq_name(c));
		ut_assert(bseq_seq_len(c) == 4, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == 0, "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == 0, "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == 0, "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test0") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x0\x0\x0") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), "") == 0, "tag(%s)", bseq_comment(c));

		c++;
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu)", bseq_name_len(c));
		ut_assert(bseq_seq_len(c) == 8, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == 0, "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == 0, "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == 0, "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test1") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), "") == 0, "tag(%s)", bseq_comment(c));

		c++;
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu)", bseq_name_len(c));
		ut_assert(bseq_seq_len(c) == 4, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == 0, "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == 0, "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == 0, "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test2") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x0\x0\x0") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), "") == 0, "tag(%s)", bseq_comment(c));

		c++;
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu)", bseq_name_len(c));
		ut_assert(bseq_seq_len(c) == 4, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == 0, "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == 0, "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == ((i & 0x02) != 0 ? 15 : 0), "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test3") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x1\x2\x3") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), (i & 0x02) != 0 ? "comment comment" : "") == 0, "tag(%s)", bseq_comment(c));

		bseq_bin_free(s);
		bseq_close_t const e = bseq_close(b);
		ut_assert(e.cnt == 4);
		ut_assert(e.status == 0);
	}
	remove(filename);
}

unittest( .name = "bseq.fastq" ) {
	char const *filename = "./bseq.unittest.tmp";
	char const *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\r\n+123\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n\r\n"
		"@  test2      \n\nAAAA\n+  test2\n\n\n12\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\r\n@\r\n1\r\n2\r\n3\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";

	FILE *fp = fopen(filename, "w");
	ut_assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	for(size_t i = 0; i < 4; i++) {
		debug("i(%zu)", i);
		bseq_conf_t const conf = {
			.batch_size   = BSEQ_BATCH_SIZE,
			.head_margin  = 64,
			.keep_qual    = (i & 0x01) != 0,
			.keep_comment = (i & 0x02) != 0,
			.conv_table = {
				[A] = 0x00,
				[C] = 0x01,
				[G] = 0x02,
				[T] = 0x03
			}
		};
		bseq_file_t *b = bseq_open(&conf, filename);
		ut_assert(b != NULL);

		bseq_batch_t *s = bseq_read(b);

		ut_assert(s != NULL);
		ut_assert(s->cnt == 4, "cnt(%u)", s->cnt);
		ut_assert(s->bin.base != NULL);
		ut_assert(s->bin.size > 0, "size(%lu)", s->bin.size);

		bseq_meta_t const *c = bseq_meta_ptr(s);
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu)", bseq_name_len(c));
		ut_assert(bseq_seq_len(c) == 4, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == ((i & 0x01) != 0 ? 4 : 0), "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == ((i & 0x01) != 0), "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == 0, "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test0") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x0\x0\x0") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), (i & 0x01) != 0 ? "NNNN" : "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), "") == 0, "tag(%s)", bseq_comment(c));

		c++;
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu)", bseq_name_len(c));
		ut_assert(bseq_seq_len(c) == 8, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == ((i & 0x01) != 0 ? 8 : 0), "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == ((i & 0x01) != 0), "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == 0, "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test1") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), (i & 0x01) != 0 ? "12+3+123" : "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), "") == 0, "tag(%s)", bseq_comment(c));

		c++;
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu)", bseq_name_len(c));
		ut_assert(bseq_seq_len(c) == 4, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == ((i & 0x01) != 0 ? 4 : 0), "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == ((i & 0x01) != 0), "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == 0, "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test2") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x0\x0\x0") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), (i & 0x01) != 0 ? "12@3" : "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), "") == 0, "tag(%s)", bseq_comment(c));

		c++;
		ut_assert(bseq_name_len(c) == 5, "nlen(%zu)", bseq_name_len(c));
		ut_assert(bseq_seq_len(c) == 4, "slen(%zu)", bseq_seq_len(c));
		ut_assert(bseq_qual_len(c) == ((i & 0x01) != 0 ? 4 : 0), "qlen(%zu)", bseq_qual_len(c));
		ut_assert(bseq_is_qual_set(c) == ((i & 0x01) != 0), "is_qual(%u)", bseq_is_qual_set(c));
		ut_assert(bseq_comment_len(c) == ((i & 0x02) != 0 ? 16 : 0), "clen(%zu)", bseq_comment_len(c));
		ut_assert(strcmp(bseq_name(c), "test3") == 0, "name(%s)", bseq_name(c));
		ut_assert(strcmp((char const *)bseq_seq(c), "\x0\x1\x2\x3") == 0, "seq(%s)", bseq_seq(c));
		ut_assert(strcmp((char const *)bseq_qual(c), (i & 0x01) != 0 ? "@123" : "") == 0, "qual(%s)", bseq_qual(c));
		ut_assert(strcmp((char const *)bseq_comment(c), (i & 0x02) != 0 ? " comment comment" : "") == 0, "tag(%s)", bseq_comment(c));

		bseq_close_t const e = bseq_close(b);
		ut_assert(e.cnt == 4);
		ut_assert(e.status == 0);
	}
	remove(filename);
}

typedef struct {
	char *seq;
	size_t len;
} bseq_unittest_seq_t;

static _force_inline
uint8_t bseq_unittest_random_char(void)
{
	switch(rand() % 4) {
		case 0:  return('A');
		case 1:  return('C');
		case 2:  return('G');
		case 3:  return('T');
		default: return('A');
	}
}

static _force_inline
bseq_unittest_seq_t bseq_unittest_rand_seq(size_t len)
{
	char *seq = malloc(sizeof(char) * (len + 1));
	for(size_t i = 0; i < len; i++) {
		seq[i] = bseq_unittest_random_char();
	}
	seq[len] = '\0';

	return((bseq_unittest_seq_t){
		.seq = seq,
		.len = len
	});
}

static _force_inline
void bseq_unittest_padding(FILE *fp, size_t min, size_t max, char c)
{
	size_t const x = (rand() % (max - min)) + min;
	for(size_t i = 0; i < x; i++) {
		fprintf(fp, "%c", c);
	}
	return;
}

static _force_inline
void bseq_unittest_dump_seq(char const *filename, bseq_unittest_seq_t const *s, size_t cnt)
{
	FILE *fp = fopen(filename, "w");
	for(size_t i = 0; i < cnt; i++) {
		fprintf(fp, ">");
		bseq_unittest_padding(fp, 0, 31, ' ');
		fprintf(fp, "%06zu comment:COMMENT comment:COMMENT;a;b;c;d-e;xyz-pqr@195+360>3", i);
		bseq_unittest_padding(fp, 0, 31, ' ');
		bseq_unittest_padding(fp, 1, 4, '\n');

		size_t rem = s[i].len;
		while(rem > 0) {
			size_t const x = (rand() % 80) + 1, frac = MIN2(rem, x);

			fprintf(fp, "%.*s", (int)frac, &s[i].seq[s[i].len - rem]);
			bseq_unittest_padding(fp, 1, 4, '\n');

			rem -= frac;
		}
	}
	fclose(fp);
	return;
}

unittest( .name = "bseq.fasta.large" ) {
	char const *filename = "./bseq.unittest.tmp";

	size_t const cnt = 2048, len = 1024;
	bseq_unittest_seq_t *arr = malloc(sizeof(bseq_unittest_seq_t) * cnt);
	for(size_t i = 0; i < cnt; i++) {
		arr[i] = bseq_unittest_rand_seq(len + (rand() % len));
	}
	bseq_unittest_dump_seq(filename, arr, cnt);

	size_t const batches[] = {
		128, 256, 801, 1000, 1024, 12345, 1ULL * 1024 * 1024, 10000000, 12345678, 16ULL * 1024 * 1024, 0
	};

	for(size_t k = 0; batches[k] != 0; k++) {
		for(size_t i = 0; i < 4; i++) {
			kvec_t(bseq_batch_t *) batch;
			kvec_t(bseq_meta_t) meta;

			kv_init(batch);
			kv_init(meta);

			bseq_conf_t const conf = {
				.batch_size   = batches[k],
				.head_margin  = 64,
				.keep_qual    = (i & 0x01) != 0,
				.keep_comment = (i & 0x02) != 0,
				.conv_table = {
					[A] = 'A',
					[C] = 'C',
					[G] = 'G',
					[T] = 'T'
				}
			};
			bseq_file_t *b = bseq_open(&conf, filename);
			bseq_batch_t *s = NULL;
			while((s = bseq_read(b)) != NULL) {
				kv_push(bseq_batch_t, batch, s);

				bseq_meta_t *m = bseq_meta_ptr(s);
				for(size_t j = 0; j < bseq_meta_cnt(s); j++) {
					kv_push(bseq_meta_t, meta, m[j]);
				}
			}

			bseq_close_t const e = bseq_close(b);
			ut_assert(e.cnt == cnt, "cnt(%zu, %zu)", e.cnt, cnt);
			ut_assert(e.status == 0);

			ut_assert(kv_cnt(meta) == cnt, "cnt(%zu, %zu)", kv_cnt(meta), cnt);
			for(size_t j = 0; j < kv_cnt(meta); j++) {
				char buf[256] = { 0 };
				sprintf(buf, "%06zu", j);

				bseq_meta_t const *c = &kv_ptr(meta)[j];
				ut_assert(bseq_name_len(c) == 6, "nlen(%zu)", bseq_name_len(c));
				ut_assert(bseq_seq_len(c)  == arr[j].len, "slen(%zu, %zu), seq(\n%s\n%s\n)", bseq_seq_len(c), arr[j].len, bseq_seq(c), arr[j].seq);
				ut_assert(bseq_qual_len(c) == 0, "qlen(%zu)", bseq_qual_len(c));
				ut_assert(bseq_is_qual_set(c) == 0, "is_qual(%u)", bseq_is_qual_set(c));
				ut_assert(bseq_comment_len(c) == ((i & 0x02) != 0 ? 59 : 0), "clen(%zu)", bseq_comment_len(c));
				ut_assert(strcmp(bseq_name(c), buf) == 0, "name(%s)", bseq_name(c));
				ut_assert(strcmp((char const *)bseq_seq(c), arr[j].seq) == 0, "seq(%s)", bseq_seq(c));
				ut_assert(strcmp((char const *)bseq_qual(c), "") == 0, "qual(%s)", bseq_qual(c));
				ut_assert(strcmp((char const *)bseq_comment(c), (i & 0x02) != 0 ? "comment:COMMENT comment:COMMENT;a;b;c;d-e;xyz-pqr@195+360>3" : "") == 0, "tag(%s)", bseq_comment(c));
			}

			for(size_t j = 0; j < kv_cnt(batch); j++) {
				bseq_bin_free(kv_ptr(batch)[j]);
			}
			kv_destroy(batch);
			kv_destroy(meta);
		}
	}
	// remove(filename);
}


/* printer */

/* custom printer for comment field */
typedef struct bseq_dump_s bseq_dump_t;
typedef size_t (*bseq_print_attr_t)(bseq_dump_t *self, bseq_meta_t const *meta);

/* context */
struct bseq_dump_s {
	struct {
		uint8_t *p, *t;		/* margined */
		size_t size;
	} buf;

	/* configurations */
	uint32_t comment;		/* print the original comment */
	uint32_t fast;			/* use fast converter (enabled when the conversion table is small) */
	uint32_t qual;			/* print qual */

	bseq_print_attr_t callback;		/* custom comment formatter */

	/* conversion table for entire ASCII space */
	uint8_t mask[16];
	uint8_t conv[256];
};

static _force_inline
void bseq_dump_init_buf(bseq_dump_t *self, size_t batch_size)
{
	size_t const size = batch_size;
	size_t const margined_size = size + 3 * BSEQ_MGN;

	/* allocate (or die) */
	void *base = aligned_malloc(margined_size);
	uint8_t *p = _add_offset(base, BSEQ_MGN);

	/* save all */
	self->buf.p = p;
	self->buf.t = p + size;
	self->buf.size = size;
	return;
}

static _force_inline
void bseq_dump_destroy_buf(bseq_dump_t *self)
{
	uint8_t *p = self->buf.t - self->buf.size;
	void *base = _sub_offset(p, BSEQ_MGN);
	free(base);
	return;
}

static _force_inline
uint8_t *bseq_dump_reserve_buf(bseq_dump_t *self, size_t len)
{
	if(_unlikely(self->buf.p + len > self->buf.t + BSEQ_MGN)) {
		uint8_t *p = self->buf.t - self->buf.size;
		uint8_t const dump_size = MIN2((size_t)(self->buf.p - p), self->buf.size);

		/* until everything is dumped */
		size_t sum = 0;
		while((sum += fwrite(p + sum, 1, dump_size - sum, stdout)) < dump_size) {}

		/* copy remaining (if there is) and reset pointer */
		_memcpy_blk_uu(p, p + dump_size, BSEQ_MGN);
		self->buf.p = p;
	}

	/* allocate buffer */
	uint8_t *q = self->buf.p;
	self->buf.p += len;
	return(q);
}


static _force_inline
void bseq_dump_load_conv(uint8_t *conv, uint8_t const *table)
{
	for(size_t i = 0x00; i < 0x0f; i++) {
		conv[table[i]] = "NACMGRSVTWYHKDBN"[i];
	}
	return;
}

static _force_inline
uint64_t bseq_dump_is_fast(uint8_t const *table)
{
	v16i8_t const tv = _loadu_v16i8(table);
	v16i8_t const fv = _set_v16i8(0x0f);

	/* if there is a cell larger than 0x0f it becomes 0xff */
	v16i8_t const gv = _gt_v16i8(tv, fv);

	/* if there is at least one cell larger than 0x0f mask becomes nonzero */
	uint64_t const gm = ((v16_masku_t){ .mask = _mask_v16i8(gv) }).all;
	return(gm == 0);
}

static _force_inline
uint64_t bseq_dump_init_static(bseq_dump_t *self, bseq_conf_t const *conf, bseq_print_attr_t callback)
{
	memset(self, 0, sizeof(bseq_dump_t));

	/* init buffer */
	size_t const batch_size = conf->batch_size == 0 ? BSEQ_BATCH_SIZE : conf->batch_size;
	bseq_dump_init_buf(self, batch_size);

	/* init conversion table */
	bseq_dump_load_conv(self->conv, conf->conv_table);

	/* configurations */
	self->comment = conf->keep_comment;
	self->qual = conf->keep_qual;
	self->fast = bseq_dump_is_fast(conf->conv_table);
	self->callback = callback;
	return(0);
}

static _force_inline
void bseq_dump_destroy_static(bseq_dump_t *self)
{
	bseq_dump_destroy_buf(self);
	return;
}



static _force_inline
size_t bseq_dump_putchar(bseq_dump_t *self, char c)
{
	uint8_t *q = bseq_dump_reserve_buf(self, 1);
	*q = c;
	return(1);
}

static _force_inline
size_t bseq_dump_putsn(bseq_dump_t *self, char const *str, size_t len)
{
	uint8_t const *p = (uint8_t const *)str, *t = p + len;
	while(p < t) {
		/* allcoate buffer */
		size_t const l = MIN2(t - p, 64);
		uint8_t *q = bseq_dump_reserve_buf(self, l);

		/* copy */
		v32i8_t const v1 = _loadu_v32i8(p);
		v32i8_t	const v2 = _loadu_v32i8(p + 32);
		_storeu_v32i8(q,      v1);
		_storeu_v32i8(q + 32, v2);

		/* forward pointer */
		p += 64;
	}
	return(len);
}

static _force_inline
size_t bseq_dump_puts(bseq_dump_t *self, char const *str)
{
	size_t const len = mm_strlen(str);
	return(bseq_dump_putsn(self, str, len));
}

static _force_inline
void bseq_dump_putchar_wrap(char c, bseq_dump_t *self)
{
	bseq_dump_putchar(self, c);
	return;
}

static _force_inline
size_t bseq_dump_printf(bseq_dump_t *self, char const *format, ...)
{
	va_list va;
	va_start(va, format);
	size_t const cnt = xvfctprintf((xp_putc_t)bseq_dump_putchar_wrap, (void *)self, format, va);
	va_end(va);
	return(cnt);
}

static _force_inline
size_t bseq_dump_header(bseq_dump_t *self, bseq_meta_t const *meta, uint8_t token)
{
	size_t bytes = 0;

	/* base name */
	bytes += bseq_dump_putchar(self, token);
	bytes += bseq_dump_putsn(self,
		bseq_name(meta),
		bseq_name_len(meta)
	);

	/* preserved comment */
	if(self->comment & (bseq_comment_len(meta) > 0)) {
		bytes += bseq_dump_putchar(self, ' ');
		bytes += bseq_dump_putsn(self,
			bseq_comment(meta),
			bseq_comment_len(meta)
		);
	}

	/* custom comment */
	if(self->callback) {
		bytes += bseq_dump_putchar(self, ' ');
		bytes += self->callback(self, meta);
	}

	/* done */
	bytes += bseq_dump_putchar(self, '\n');
	return(bytes);
}


static _force_inline
size_t bseq_dump_core_fast(bseq_dump_t *self, uint8_t const *seq, size_t slen)
{
	v32i8_t const conv = _loadu_v32i8(self->conv);
	uint8_t const *p = (uint8_t const *)seq, *t = p + slen;

	while(p < t) {
		/* allcoate buffer */
		size_t const l = MIN2(t - p, 64);
		uint8_t *q = bseq_dump_reserve_buf(self, l);

		/* convert and save */
		v32i8_t const v1 = _loadu_v32i8(p);
		v32i8_t const v2 = _loadu_v32i8(p + 32);
		v32i8_t const w1 = _shuf_v32i8(conv, v1);
		v32i8_t const w2 = _shuf_v32i8(conv, v2);
		_storeu_v32i8(q,      w1);
		_storeu_v32i8(q + 32, w2);

		/* forward pointer */
		p += 64;
	}
	return(slen);
}

static _force_inline
size_t bseq_dump_core_slow(bseq_dump_t *self, uint8_t const *seq, size_t slen)
{
	uint8_t const *conv = self->conv;
	uint8_t const *p = (uint8_t const *)seq, *t = p + slen;

	/* really slow... */
	while(p < t) {
		uint8_t const c = *p++;
		uint8_t const x = conv[c];

		uint8_t *q = bseq_dump_reserve_buf(self, 1);
		*q = x;
	}
	return(slen);
}

static _force_inline
size_t bseq_dump_seq(bseq_dump_t *self, bseq_meta_t const *meta)
{
	size_t bytes = 0;

	/* determine record type first */
	uint64_t const qual = self->qual & bseq_is_qual_set(meta);
	uint8_t const token = qual ? '@' : '>';

	/* header */
	bytes += bseq_dump_header(self, meta, token);

	/* body */
	bytes += (self->fast ? bseq_dump_core_fast : bseq_dump_core_slow)(self,
		bseq_seq(meta),
		bseq_seq_len(meta)
	);
	if(qual == 0) { return(bytes); }

	/* qual */
	bytes += bseq_dump_putchar(self, '+');
	bytes += bseq_dump_putchar(self, '\n');
	bytes += bseq_dump_putsn(self,
		(char const *)bseq_qual(meta),
		bseq_qual_len(meta)
	);
	bytes += bseq_dump_putchar(self, '\n');
	return(bytes);
}

static _force_inline
size_t bseq_dump_batch(bseq_dump_t *self, bseq_batch_t const *batch)
{
	/* load seq array */
	bseq_meta_t const *seq = bseq_meta_ptr(batch);
	size_t const scnt = bseq_meta_cnt(batch);

	/* dump all */
	size_t bytes = 0;
	for(size_t i = 0; i < scnt; i++) {
		bytes += bseq_dump_seq(self, &seq[i]);
	}
	return(bytes);
}


/**
 * end of bseq.h
 */
