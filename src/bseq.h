
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

#define bseq_is_qual_set(_x)		( bseq_qual(_x)[0] == '\0' )
#define bseq_qual(_x)				( (uint8_t *)&(_x)->seq[(_x)->slen] )
#define bseq_qual_len(_x)			( bseq_is_qual_set(_x) ? 0 : (size_t)(_x)->slen )

#define bseq_name(_x)				( (char const *)&(_x)->seq[-((size_t)(_x)->hlen)] )
#define bseq_name_len(_x)			( (size_t)(_x)->nlen )
#define bseq_comment(_x)			( (char const *)&(_x)->seq[-((size_t)(_x)->hlen) + (_x)->nlen + 1] )
#define bseq_comment_len(_x)		( (size_t)(_x)->hlen - (size_t)(_x)->nlen - 2 )

#define bseq_circ(_x)				( (_x)->u.flag )
#define bseq_aln(_x)				( (_x)->u.ptr )

/* memory management */

#define BSEQ_MGN					( 64 )			/* buffer margin length */

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
	size_t malloc_size = _roundup(sizeof(uint8_t) * batch_size + 2 * BSEQ_MGN, ARCH_HUGEPAGE_SIZE);
	uint8_t *ptr = aligned_malloc(malloc_size);
	_memset_blk_u(ptr, 0, BSEQ_MGN);		/* clear head margin so that parsing will not depend on uninitialized values */

	/* save */
	fp->batch_size = batch_size;
	fp->buf = ptr + BSEQ_MGN;
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
	size_t read = rbread_bulk(&fp->rb, fp->buf, fp->batch_size);		/* always bypass internal buffer; up to batch_size */

	/* save pointers */
	fp->p = fp->buf;
	fp->t = fp->buf + read;

	/* add padding */
	_memset_blk_u(fp->t, '\n', BSEQ_MGN);		/* fill margin of the input buffer */

	/* update EOF */
	fp->is_eof = MAX2(fp->is_eof, rbeof(&fp->rb));

	// debug("read, state(%u), len(%lu), is_eof(%u)", fp->state, read, fp->is_eof);
	return(read);
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
	// debug("return, p(%p), t(%p)", w->p, w->t);
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
	// debug("single(%p), vector(%p), p(%p), t(%p), n(%lu), c(%x)", single, vector, w->p, w->t, w->n, w->q[w->n]);

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
	size_t back = w->p[-1] == '\r';
	size_t read = (w->p - b) - back;

	/* forward pointer (skip an obvious delimiter at the end) */
	w->p++;

	/* fixup last byte for the next readline */
	w->q[w->n] = '\0';				/* mark empty; invariant condition */
	return(read);
}

static _force_inline
void bseq_fixup_name(bseq_work_t *w, bseq_file_t *fp)
{
	size_t nlen = fp->acc, p = w->n - nlen, t = w->n;
	v32i8_t const sv = _set_v32i8(' ');

	/* strip head spaces */ {
		uint64_t m = _match(sv, _loadu_v32i8(&w->q[p]));
		ZCNT_RESULT size_t l = _tzc_u32(~m);
		size_t len = l; p += len; nlen -= len;
	}

	/* strip tail spaces */ {
		uint64_t m = _match(sv, _loadu_v32i8(&w->q[t - 32]));
		ZCNT_RESULT size_t l = _lzc_u32(~m);
		size_t len = l; t -= len; nlen -= len;
	}

	/* split at the first space */
	size_t q = p;
	while(q < t) {
		uint64_t m = _match(sv, _loadu_v32i8(&w->q[q]));
		size_t len = MIN2(32, t - q);
		q += len;
		if(m) {
			ZCNT_RESULT size_t u = _tzc_u32(m);
			q += u - len; break;
		}
	}

	/* remove comment if needed */
	if(fp->keep_comment == 0) { t = q + 1; }

	/* fill terminators */
	w->s->nlen = q - p;
	w->s->hlen = (t + 1) - p;				/* might have further spaces at the head of comment */
	// debug("nlen(%u), hlen(%u)", w->s->nlen, w->s->hlen);

	w->q[q] = '\0';
	w->q[t] = '\0';							/* string terminator; keep q[n] initialized (invariant condition) */
	w->n = t;

	// debug("name(%lu, %s), comment(%lu, %s), n(%lu)", bseq_name_len(w->s), &w->q[p], bseq_comment_len(w->s), &w->q[q + 1], w->n);
	return;
}

static _force_inline
void bseq_fixup_seq(bseq_work_t *w, bseq_file_t *fp)
{
	w->s->seq   = (uint8_t *)w->n - fp->acc;
	w->s->slen  = fp->acc;
	w->s->u.ptr = NULL;
	// debug("seq(%p), slen(%lu)", w->s->seq, w->s->slen);
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
		/* debug("forward state(%u)", _n); */ \
		fp->state = _n; \
		fp->acc = 0; \
		case _n: /* debug("enter state(%u)", _n); */ \
		while(1) { \
			_body; \
			uint64_t flag = (_cond); \
			if(w.p >= w.t) { goto _refill; } \
			if(!flag) { break; } \
		} \
		_post; \
		w.q[++w.n] = 0;		/* invariant condition; open new empty byte for the next parsing */ \
	}

	/* error if buffer is empty */
	if(fp->p >= fp->t) { return(1); }

	/* keep them on registers */
	bseq_work_t w __attribute__(( aligned(64) ));
	bseq_init_work(&w, fp, seq, bin);

	// debug("enter, state(%u), p(%p), t(%p), eof(%u)", fp->state, w.p, w.t, fp->is_eof);
	switch(fp->state) {							/* dispatcher for the first iteration */
		default:								/* idle or broken */
		if(*w.p++ != fp->delim[0]) { return(0); }	/* '>' */
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
		}
		fp->state = 0;
		// debug("done");
	}
_refill:;
	// debug("break, state(%u), eof(%u)", fp->state, fp->is_eof);
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
			if(bseq_fetch_block(fp) == 0) { break; }/* fetch next */
			
			/* reserve room for the next parsing unit */
			kvm_reserve(uint8_t, *bin, kvm_cnt(*bin) + 2 * fp->batch_size);
			// debug("continue, state(%u), is_eof(%u), seq(%p, %lu), bin(%p, %lu)", fp->state, fp->is_eof, kvm_ptr(*seq), kvm_cnt(*seq), kvm_ptr(*bin), kvm_cnt(*bin));
		}
		// debug("finished seq, state(%u), is_eof(%u), seq(%p, %lu), bin(%p, %lu)", fp->state, fp->is_eof, kvm_ptr(*seq), kvm_cnt(*seq), kvm_ptr(*bin), kvm_cnt(*bin));

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
	ptrdiff_t b = (ptrdiff_t)kvm_ptr(*bin);
	kvm_foreach(bseq_meta_t, *seq, { p->seq += b; });

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

	size_t scnt = fp->scnt;
	uint64_t status = bseq_is_broken(fp);

	rbclose_static(&fp->rb);
	bseq_buf_destroy(fp);
	free(fp);
	return((bseq_close_t){
		.cnt = scnt,
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
		for(size_t i = 0x01; i < 0x0f; i++) {
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
bseq_file_t *bseq_open(bseq_conf_t const *conf, char const *fn)	/* file path; "-" to use stdin */
{
	/* create instance */
	bseq_file_t *fp = (bseq_file_t *)calloc(1, sizeof(bseq_file_t));
	*fp = (bseq_file_t){
		.keep_qual = conf->keep_qual,
		.keep_comment = conf->keep_comment,
		.head_margin = _roundup(conf->head_margin, sizeof(bseq_batch_t))
	};

	/* exclude trivial failure */
	if(fn == NULL || fn[0] == '\0') {
		error("empty file path `'.");
		goto _bseq_open_fail;
	}

	/* open stream */
	if(rbopen_bulk_static(&fp->rb, fn, conf->batch_size) != 0) {
		error("failed to open file `%s'.", fn);
		goto _bseq_open_fail;
	}

	/* build conversion table */
	bseq_load_conv(fp->conv, conf->conv_table);

	/* init buffer; buffer size shoule be the same to the rb internal; fetch the first block */
	bseq_buf_init(fp, conf->batch_size);
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

#if 0
unittest( .name = "bseq.fasta" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\r\n\r\n"
		">  test2\n\nAAAA\n"
		">test3 comment comment  \nACGT\n\n";

	FILE *fp = fopen(filename, "w");
	ut_assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 64, 1, 0, 1, tags);
	ut_assert(b != NULL);
	bseq_meta_t *s = q(b);

	ut_assert(s != NULL);
	ut_assert(s->scnt == 4, "scnt(%u)", s->scnt);
	ut_assert(s->base != NULL);
	ut_assert(s->size > 0, "size(%lu)", s->size);

	ut_assert(s->seq[0].nlen == 5, "nlen(%u)", s->seq[0].nlen);
	ut_assert(s->seq[0].slen == 4, "slen(%u)", s->seq[0].slen);
	ut_assert(s->seq[0].n_tag == 0, "n_tag(%u)", s->seq[0].n_tag);
	ut_assert(strcmp((const char*)s->seq[0].name, "test0") == 0, "name(%s)", s->seq[0].name);
	ut_assert(strcmp((const char*)s->seq[0].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[0].seq);
	ut_assert(strcmp((const char*)s->seq[0].qual, "") == 0, "qual(%s)", s->seq[0].qual);
	ut_assert(strcmp((const char*)s->seq[0].tag, "") == 0, "tag(%s)", s->seq[0].tag);

	ut_assert(s->seq[1].nlen == 5, "nlen(%u)", s->seq[1].nlen);
	ut_assert(s->seq[1].slen == 8, "slen(%u)", s->seq[1].slen);
	ut_assert(s->seq[1].n_tag == 0, "n_tag(%u)", s->seq[1].n_tag);
	ut_assert(strcmp((const char*)s->seq[1].name, "test1") == 0, "name(%s)", s->seq[1].name);
	ut_assert(strcmp((const char*)s->seq[1].seq, "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", s->seq[1].seq);
	ut_assert(strcmp((const char*)s->seq[1].qual, "") == 0, "qual(%s)", s->seq[1].qual);
	ut_assert(strcmp((const char*)s->seq[1].tag, "") == 0, "tag(%s)", s->seq[1].tag);

	ut_assert(s->seq[2].nlen == 5, "nlen(%u)", s->seq[2].nlen);
	ut_assert(s->seq[2].slen == 4, "slen(%u)", s->seq[2].slen);
	ut_assert(s->seq[2].n_tag == 0, "n_tag(%u)", s->seq[2].n_tag);
	ut_assert(strcmp((const char*)s->seq[2].name, "test2") == 0, "name(%s)", s->seq[2].name);
	ut_assert(strcmp((const char*)s->seq[2].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[2].seq);
	ut_assert(strcmp((const char*)s->seq[2].qual, "") == 0, "qual(%s)", s->seq[2].qual);
	ut_assert(strcmp((const char*)s->seq[2].tag, "") == 0, "tag(%s)", s->seq[2].tag);

	ut_assert(s->seq[3].nlen == 5, "nlen(%u)", s->seq[3].nlen);
	ut_assert(s->seq[3].slen == 4, "slen(%u)", s->seq[3].slen);
	ut_assert(s->seq[3].n_tag == 1, "n_tag(%u)", s->seq[3].n_tag);
	ut_assert(strcmp((const char*)s->seq[3].name, "test3") == 0, "name(%s)", s->seq[3].name);
	ut_assert(strcmp((const char*)s->seq[3].seq, "\x0\x1\x2\x3") == 0, "seq(%s)", s->seq[3].seq);
	ut_assert(strcmp((const char*)s->seq[3].qual, "") == 0, "qual(%s)", s->seq[3].qual);
	ut_assert(strcmp((const char*)s->seq[3].tag, "COZcomment comment") == 0, "tag(%s)", s->seq[3].tag);

	uint64_t scnt = bseq_close(b);
	ut_assert(scnt == 4);
	remove(filename);
}

unittest( .name = "bseq.fastq" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\n+123\r\n"
		"@  test2\n\nAAAA\n+  test2\n\n\n12@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\n@123";

	FILE *fp = fopen(filename, "w");
	ut_assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 256, 1, 0, 1, tags);
	ut_assert(b != NULL);
	bseq_meta_t *s = q(b);

	ut_assert(s != NULL);
	ut_assert(s->scnt == 4, "scnt(%u)", s->scnt);
	ut_assert(s->base != NULL);
	ut_assert(s->size > 0, "size(%lu)", s->size);

	ut_assert(s->seq[0].nlen == 5, "nlen(%u)", s->seq[0].nlen);
	ut_assert(s->seq[0].slen == 4, "slen(%u)", s->seq[0].slen);
	ut_assert(s->seq[0].n_tag == 0, "n_tag(%u)", s->seq[0].n_tag);
	ut_assert(strcmp((const char*)s->seq[0].name, "test0") == 0, "name(%s)", s->seq[0].name);
	ut_assert(strcmp((const char*)s->seq[0].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[0].seq);
	ut_assert(strcmp((const char*)s->seq[0].qual, "NNNN") == 0, "qual(%s)", s->seq[0].qual);
	ut_assert(strcmp((const char*)s->seq[0].tag, "") == 0, "tag(%s)", s->seq[0].tag);

	ut_assert(s->seq[1].nlen == 5, "nlen(%u)", s->seq[1].nlen);
	ut_assert(s->seq[1].slen == 8, "slen(%u)", s->seq[1].slen);
	ut_assert(s->seq[1].n_tag == 0, "n_tag(%u)", s->seq[1].n_tag);
	ut_assert(strcmp((const char*)s->seq[1].name, "test1") == 0, "name(%s)", s->seq[1].name);
	ut_assert(strcmp((const char*)s->seq[1].seq, "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", s->seq[1].seq);
	ut_assert(strcmp((const char*)s->seq[1].qual, "12+3+123") == 0, "qual(%s)", s->seq[1].qual);
	ut_assert(strcmp((const char*)s->seq[1].tag, "") == 0, "tag(%s)", s->seq[1].tag);

	ut_assert(s->seq[2].nlen == 5, "nlen(%u)", s->seq[2].nlen);
	ut_assert(s->seq[2].slen == 4, "slen(%u)", s->seq[2].slen);
	ut_assert(s->seq[2].n_tag == 0, "n_tag(%u)", s->seq[2].n_tag);
	ut_assert(strcmp((const char*)s->seq[2].name, "test2") == 0, "name(%s)", s->seq[2].name);
	ut_assert(strcmp((const char*)s->seq[2].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[2].seq);
	ut_assert(strcmp((const char*)s->seq[2].qual, "12@3") == 0, "qual(%s)", s->seq[2].qual);
	ut_assert(strcmp((const char*)s->seq[2].tag, "") == 0, "tag(%s)", s->seq[2].tag);

	ut_assert(s->seq[3].nlen == 5, "nlen(%u)", s->seq[3].nlen);
	ut_assert(s->seq[3].slen == 4, "slen(%u)", s->seq[3].slen);
	ut_assert(s->seq[3].n_tag == 1, "n_tag(%u)", s->seq[3].n_tag);
	ut_assert(strcmp((const char*)s->seq[3].name, "test3") == 0, "name(%s)", s->seq[3].name);
	ut_assert(strcmp((const char*)s->seq[3].seq, "\x0\x1\x2\x3") == 0, "seq(%s)", s->seq[3].seq);
	ut_assert(strcmp((const char*)s->seq[3].qual, "@123") == 0, "qual(%s)", s->seq[3].qual);
	ut_assert(strcmp((const char*)s->seq[3].tag, "COZcomment comment") == 0, "tag(%s)", s->seq[3].tag);

	uint64_t scnt = bseq_close(b);
	ut_assert(scnt == 4);
	remove(filename);
}

unittest( .name = "bseq.fastq.skip" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\n+123\n"
		"@  test2\n\nAAAA\n+  test2\n\n\n12@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\n@123";

	FILE *fp = fopen(filename, "w");
	ut_assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 256, 0, 0, 1, tags);
	ut_assert(b != NULL);
	bseq_meta_t *s = q(b);

	ut_assert(s != NULL);
	ut_assert(s->scnt == 4, "scnt(%u)", s->scnt);
	ut_assert(s->base != NULL);
	ut_assert(s->size > 0, "size(%lu)", s->size);

	ut_assert(s->seq[0].nlen == 5, "nlen(%u)", s->seq[0].nlen);
	ut_assert(s->seq[0].slen == 4, "slen(%u)", s->seq[0].slen);
	ut_assert(s->seq[0].n_tag == 0, "n_tag(%u)", s->seq[0].n_tag);
	ut_assert(strcmp((const char*)s->seq[0].name, "test0") == 0, "name(%s)", s->seq[0].name);
	ut_assert(strcmp((const char*)s->seq[0].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[0].seq);
	ut_assert(strcmp((const char*)s->seq[0].qual, "") == 0, "qual(%s)", s->seq[0].qual);
	ut_assert(strcmp((const char*)s->seq[0].tag, "") == 0, "tag(%s)", s->seq[0].tag);

	ut_assert(s->seq[1].nlen == 5, "nlen(%u)", s->seq[1].nlen);
	ut_assert(s->seq[1].slen == 8, "slen(%u)", s->seq[1].slen);
	ut_assert(s->seq[1].n_tag == 0, "n_tag(%u)", s->seq[1].n_tag);
	ut_assert(strcmp((const char*)s->seq[1].name, "test1") == 0, "name(%s)", s->seq[1].name);
	ut_assert(strcmp((const char*)s->seq[1].seq, "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", s->seq[1].seq);
	ut_assert(strcmp((const char*)s->seq[1].qual, "") == 0, "qual(%s)", s->seq[1].qual);
	ut_assert(strcmp((const char*)s->seq[1].tag, "") == 0, "tag(%s)", s->seq[1].tag);

	ut_assert(s->seq[2].nlen == 5, "nlen(%u)", s->seq[2].nlen);
	ut_assert(s->seq[2].slen == 4, "slen(%u)", s->seq[2].slen);
	ut_assert(s->seq[2].n_tag == 0, "n_tag(%u)", s->seq[2].n_tag);
	ut_assert(strcmp((const char*)s->seq[2].name, "test2") == 0, "name(%s)", s->seq[2].name);
	ut_assert(strcmp((const char*)s->seq[2].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[2].seq);
	ut_assert(strcmp((const char*)s->seq[2].qual, "") == 0, "qual(%s)", s->seq[2].qual);
	ut_assert(strcmp((const char*)s->seq[2].tag, "") == 0, "tag(%s)", s->seq[2].tag);

	ut_assert(s->seq[3].nlen == 5, "nlen(%u)", s->seq[3].nlen);
	ut_assert(s->seq[3].slen == 4, "slen(%u)", s->seq[3].slen);
	ut_assert(s->seq[3].n_tag == 1, "n_tag(%u)", s->seq[3].n_tag);
	ut_assert(strcmp((const char*)s->seq[3].name, "test3") == 0, "name(%s)", s->seq[3].name);
	ut_assert(strcmp((const char*)s->seq[3].seq, "\x0\x1\x2\x3") == 0, "seq(%s)", s->seq[3].seq);
	ut_assert(strcmp((const char*)s->seq[3].qual, "") == 0, "qual(%s)", s->seq[3].qual);
	ut_assert(strcmp((const char*)s->seq[3].tag, "COZcomment comment") == 0, "tag(%s)", s->seq[3].tag);

	uint64_t scnt = bseq_close(b);
	ut_assert(scnt == 4);
	remove(filename);
}
#endif

/**
 * end of bseq.h
 */
