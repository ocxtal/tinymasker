
/**
 * @file dbg.c
 * @brief de Bruijn hash structure for reference sequences
 *
 * @author Hajime Suzuki
 * @license MIT
 */


#ifndef UNITTEST
#  define UNITTEST				( 1 )
#endif
#define UNITTEST_UNIQUE_ID		5


#include "utils/utils.h"		/* include all */
#include "tinymasker.h"
unittest_config( .name = "dbg" );


/* forward decls */
#include "dbg.h"


/* reference index object */
// static _force_inline
void tm_ref_destroy(tm_ref_sketch_t *sk)
{
	free(_sub_offset(sk, sk->head_margin));
	return;
}

/* for debugging */
static char *tm_kmer_to_str(uint32_t const *kmer)
{
	char buf[12] = { 0 };
	for(size_t i = 0; i < 11; i++) {
		buf[i] = "ACGT"[(*kmer>>(2 * (10 - i))) & 0x03];
	}
	return(strdup(buf));
}


/* construct reference object */
typedef union {
	struct {
		uint32_t next : 8;
		uint32_t pos  : 24;
		uint32_t kmer;
	} s;
	uint64_t all;
} tm_ref_kpos_t;
typedef struct { tm_ref_kpos_t *a; size_t n, m; } tm_ref_kpos_v;
_static_assert(sizeof(tm_ref_kpos_t) == 8);

#define tm_ref_kmer(x)			( (x).all )
KRADIX_SORT_INIT(kmer, tm_ref_kpos_t, tm_ref_kmer, 8);	/* sort all */

typedef struct {
	uint32_t f, r;
} tm_ref_kmer_pair_t;

/* branch stack; max ambiguity = 4^5 */
#define TM_REF_SHIFT_OFFSET		( 6 )
#define TM_REF_KMER_STACK_SIZE	( 0x01ULL<<TM_MAX_AMB_COLUMNS )

typedef struct {
	/* constants */
	uint32_t amask, shift;

	/* dst array */
	tm_ref_kpos_t *q;

	/* kmer stack */
	uint32_t size, max_size;
	uint32_t branches, unused;
	tm_ref_kmer_pair_t kmer[TM_REF_KMER_STACK_SIZE];
} tm_ref_work_t;


static _force_inline
size_t tm_ref_dup_stack(tm_ref_work_t *w)
{
	size_t const prev_size = w->size;

	/* update size and branch pos */
	w->size      = w->size<<1;
	w->branches |= 0x01<<w->shift;

	for(size_t i = 0; i < prev_size; i++) {
		w->kmer[prev_size + i] = w->kmer[i];
	}
	return(prev_size);
}

static _force_inline
uint64_t tm_ref_update_stack(tm_ref_work_t *w)
{
	uint64_t const shrink = (w->branches & (0x01<<TM_REF_SHIFT_OFFSET)) != 0;
	w->branches = w->branches>>2;
	return(shrink);
}

static _force_inline
void tm_ref_shrink_stack(tm_ref_work_t *w)
{
	debug("shrink");
	w->size = w->size>>1;
	for(size_t i = 0; i < w->size; i++) {
		w->kmer[i] = w->kmer[2 * i];
	}
	return;
}

static _force_inline
void tm_ref_update_kmer(tm_ref_kmer_pair_t *kmer, size_t size, uint32_t amask, uint32_t shift, uint8_t base)
{
	uint32_t const f = base<<TM_REF_SHIFT_OFFSET;
	uint32_t const r = (base ^ 0x03)<<shift;
	for(size_t i = 0; i < size; i++) {
		kmer[i].f = ((kmer[i].f<<2) | f) & amask;
		kmer[i].r = ((kmer[i].r>>2) | r) & amask;
	}
	return;
}

static _force_inline
v4i32_t tm_ref_compose_kpos(v2i32_t pos, v2i32_t kmer)
{
	static uint8_t const shuf[16] __attribute__(( aligned(16) )) = {
		0, 4, 5, 6, 1, 2, 3, 0x80, 8, 12, 13, 14, 9, 10, 11, 0x80
	};

	/* fw in lower, rv in upper */
	v16i8_t const sv = _load_v16i8(shuf);
	v16i8_t const v = (v16i8_t){
		_mm_unpacklo_epi32(kmer.v1, pos.v1)
	};
	v16i8_t const w = _shuf_v16i8(v, sv);

	// _print_v16i8(v);
	// _print_v16i8(w);
	return((v4i32_t){ w.v1 });
}

static _force_inline
void tm_ref_push_pos(tm_ref_work_t *w, v2i32_t pos)
{
	for(size_t i = 0; i < w->size; i++) {
		v2i32_t const kmer = _loadu_v2i32(&w->kmer[i]);
		v4i32_t const kpos = tm_ref_compose_kpos(pos, kmer);
		_storeu_v4i32(w->q, kpos);
		w->q += 2;
	}
	return;
}

static _force_inline
void tm_ref_push_base(tm_ref_work_t *w, uint8_t base)
{
	size_t const size = w->size;		/* save current size before duplicating stack */
	uint8_t const b2 = tm_rch_to_2bit(base);

	// debug("base(%x, %c), ambiguous(%lu), size(%zu), max_size(%zu)", base, tm_rch_to_ascii(base), tm_rch_is_amb(base), w->size, w->max_size);
	if(_unlikely(tm_rch_is_amb(base) & (w->size < w->max_size))) {
		tm_ref_dup_stack(w);
		tm_ref_update_kmer(&w->kmer[size], size, w->amask, w->shift, b2>>2);
	}

	/* push base */
	tm_ref_update_kmer(&w->kmer[0], size, w->amask, w->shift, b2 & 0x03);	
	// debug("branches(%lx), b2(%x)", w->branches, b2);
	return;
}

static _force_inline
tm_ref_kpos_t *tm_ref_reserve(tm_ref_kpos_v *buf, tm_ref_kpos_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);
	// debug("cnt(%zu)", kv_cnt(*buf));

	size_t const kcnt = kv_cnt(*buf);
	tm_ref_kpos_t *p = kv_reserve(tm_ref_kpos_t, *buf, kcnt + 32 * TM_REF_KMER_STACK_SIZE);
	return(&p[kcnt]);
}

static _force_inline
size_t tm_ref_enumerate_kmer(size_t kbits, size_t max_amb_col, uint8_t const *seq, size_t slen, tm_ref_kpos_v *buf)
{
	tm_ref_work_t w = {
		/* place kmers at 6..kbits + 6 */
		.amask = ((0x01<<kbits) - 0x01)<<TM_REF_SHIFT_OFFSET,
		.shift = kbits - 2 + TM_REF_SHIFT_OFFSET,			/* k - 1 + offset */

		/* make room in dst array */
		.q = tm_ref_reserve(buf, kv_ptr(*buf)),

		/* clear stack; keeps all the combination of ambiguous bases */
		.size     = 1,
		.max_size = 1ULL<<max_amb_col,
		.branches = 0,
		.kmer = { { 0 } }
	};
	// debug("kbits(%zu), seq(%p), slen(%zu)", kbits, seq, slen);

	size_t const k = kbits>>1;
	v2i32_t const inc = _seta_v2i32(-1, 1);		/* rv, fw */
	v2i32_t pos = _seta_v2i32(k - TM_SEED_POS_OFS, TM_SEED_POS_OFS);

	/* push bases at the head */
	for(size_t i = 0; i < k; i++) {
		tm_ref_push_base(&w, seq[i]);
		tm_ref_update_stack(&w);
		pos = _add_v2i32(pos, inc);
	}

	/* body */
	for(size_t i = k; i < slen; i++) {
		tm_ref_push_base(&w, seq[i]);			/* stack will be expanded if needed */

		/* update offset; position of the next base to compose [,) range */
		pos = _add_v2i32(pos, inc);
		tm_ref_push_pos(&w, pos);

		/* update branch stack */
 		if(_unlikely(tm_ref_update_stack(&w))) {
			tm_ref_shrink_stack(&w);			/* shrink if needed */
		}

		/*
		 * check need for expansion every 32times
		 * (Skylake's branch predictor successfully learns once-every-32times event, but fails for more)
		 */
		if((i & 0x1f) != 0) { continue; }
		w.q = tm_ref_reserve(buf, w.q);
	}

	/* tail */
	tm_ref_push_base(&w, tm_rch_pack(tA));

	pos = _add_v2i32(pos, inc);
	tm_ref_push_pos(&w, pos);

	/* done; record count */
	kv_cnt(*buf) = w.q - kv_ptr(*buf);
	return(kv_cnt(*buf));
}

static _force_inline
size_t tm_ref_dedup_kmer(tm_ref_kpos_t *kptr, size_t kcnt)
{
	tm_ref_kpos_t const *p = kptr, *t = &kptr[kcnt];	/* src */
	tm_ref_kpos_t *q = kptr + 1;						/* dst */

	uint64_t const mask = ~0xffULL;		/* clear out next field */
	uint64_t prev = _loadu_u64(p) & mask;
	while(++p < t) {
		uint64_t const kmer = _loadu_u64(p);
		uint64_t const curr = kmer & mask;

		if(curr == prev) {				/* skip if both pos and kmer are the same as the previous */
			debug("dedup, i(%zu), prev(%lx), curr(%lx, %lx)", p - kptr, prev, curr, kmer);
			continue;
		}

		/* pos or kmer changed; save */
		_storeu_u64(q, kmer);
		q++;

		/* update for next */
		prev = curr;
	}
	return(q - kptr);
}

static _force_inline
size_t tm_ref_collect_kmer(tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen, tm_ref_kpos_v *buf)
{
	kv_clear(*buf);

	/* slice k-mers (with additional succeeding base) */
	if(tm_ref_enumerate_kmer(conf->kbits + 2, conf->max_amb_col, seq, slen, buf) == 0) {
		return(0);				/* abort if no kmer found */
	}

	/* sort kmers by position */
	tm_ref_kpos_t *kptr = kv_ptr(*buf);
	size_t const kcnt = kv_cnt(*buf);
	radix_sort_kmer(kptr, kcnt);

	/* dedup kmers around ambiguous bases */
	debug("kcnt(%zu)", kcnt);
	size_t const dcnt = tm_ref_dedup_kmer(kptr, kcnt);
	kv_cnt(*buf) = dcnt;
	return(dcnt);
}


/* pack kmer-position array */
typedef struct {
	uint8_t covered, exist;		/* covered flag for feeder */
	uint16_t cnt;
	uint32_t ofs;
} tm_ref_cnt_t;
typedef struct { tm_ref_cnt_t *a; size_t n, m; } tm_ref_cnt_v;


static _force_inline
size_t tm_ref_count_kmer(tm_ref_cnt_t *p, size_t ksize, tm_ref_kpos_t const *kpos, size_t kcnt)
{
	_unused(ksize);

	// debug("ksize(%zx), kcnt(%zu)", ksize, kcnt);
	for(tm_ref_kpos_t const *k = kpos, *t = &kpos[kcnt]; k < t; k++) {
		tm_ref_cnt_t *q = &p[k->s.kmer];

		q->covered |= 1;
		q->exist   |= 1;
		q->cnt++;
	}
	return(kcnt);
}

static _force_inline
void tm_ref_patch_feeder(tm_ref_cnt_t *p, size_t ksize)
{
	for(size_t i = 0; i < ksize; i++) {
		// debug("i(%r), ksize(%zx)", tm_kmer_to_str, &i, ksize - 1);
		if(p[i].exist == 0) { continue; }

		for(uint64_t mask = (ksize - 1)>>2; mask != 0; mask >>= 2) {
			uint64_t const krem = i & mask;
			p[krem].covered |= 1;

			// debug("test mask(%lx), krem(%r), covered(%u), exist(%u)", mask, tm_kmer_to_str, &krem, p[krem].covered, p[krem].exist);
		}
	}

	size_t patched = 0;
	for(size_t i = 0; i < (ksize>>2); i++) {
		if(p[i].covered == 0) {
			p[i].exist |= 1;
			patched++;
		}
	}
	debug("ksize(%zu), patched(%zu), ratio(%f)", ksize, patched, 4.0 * (double)patched / (double)ksize);
	return;
}

static _force_inline
size_t tm_ref_calc_size(tm_ref_cnt_t const *p, size_t ksize)
{
	size_t const hdr = sizeof(tm_ref_bin_t);

	size_t acc = TM_REF_BASE_OFS;
	for(size_t i = 0; i < ksize; i++) {
		if(p[i].exist == 0) { continue; }
		acc += _roundup(hdr + sizeof(uint16_t) * p[i].cnt, TM_REF_ALIGN_SIZE);

		// debug("i(%zx), cnt(%zu), acc(%zu)", i, p[i].cnt, acc);
	}
	return(acc);
}

static _force_inline
tm_ref_kpos_t const *tm_ref_pack_kpos_core(tm_ref_bin_t *bin, uint64_t kmer, tm_ref_kpos_t const *k, tm_ref_kpos_t const *t)
{
	/* extract kmer and next fields */
	uint64_t const kmask = 0xffffffff000000ff;

	/* first clear header; sizeof(tm_ref_bin_t) == sizeof(v4i32) */
	_storeu_v4i32(bin, _zero_v4i32());

	/* skip missing k-mers (won't be executed) */
	while(k < t && k->s.kmer < kmer) { trap(); k++; }

	/* pack pos array */
	size_t x = 0, y = 0;
	while(k < t && k->s.kmer == kmer) {
		uint64_t const next = k->s.next>>6;
		uint64_t const kall = k->all & kmask;
		while(y < next) {
			bin->link[y++].tail = x;
		}
		while(k < t && (k->all & kmask) == kall) {
			bin->pos[x++] = k++->s.pos;
		}
	}
	/* fill missing tail */
	while(y < 4) { bin->link[y++].tail = x; }

	/* clear tail margin */
	_storeu_v16i8(&bin->pos[x], _zero_v16i8());		/* suppress use-of-uninitialized-value error in valgrind */
	return(k);
}

static _force_inline
size_t tm_ref_pack_kpos(tm_ref_cnt_t *p, size_t ksize, tm_ref_kpos_t const *kpos, size_t kcnt, tm_ref_sketch_t *sk)
{
	_unused(ksize);

	/* first bin offset at base_ofs */
	size_t ofs = TM_REF_BASE_OFS;

	tm_ref_kpos_t const *k = kpos, *t = &kpos[kcnt];
	for(size_t i = 0; i < ksize && k < t; i++) {
		if(p[i].exist == 0) { continue; }

		/* save current offset; slice bin from the offset */
		p[i].ofs = ofs;
		tm_ref_bin_t *bin = _add_offset(sk, ofs);

		/* pack */
		k = tm_ref_pack_kpos_core(bin, i, k, t);

		/* update offset */
		ofs += _roundup(
			sizeof(tm_ref_bin_t) + sizeof(uint16_t) * p[i].cnt,
			TM_REF_ALIGN_SIZE
		);
	}
	return(ofs);
}

typedef struct {
	uint32_t len, ofs;

	/* for debugging */
	#ifdef DEBUG_KMER
		uint32_t kmer;
	#endif
} tm_ref_prefix_t;

static _force_inline
tm_ref_prefix_t tm_ref_find_link(tm_ref_cnt_t const *p, uint32_t kmer, size_t max_shift)
{
	/* breadth-first search on kmer branching tree over profile table */
	size_t const ksize = 0x01<<max_shift;
	size_t shift = max_shift;		/* where (max_shift - shift) / 2 is substitution length (zero for the first iteration) */
	do {
		size_t const kpitch = 0x01<<shift;
		uint32_t const krem = kmer & (kpitch - 0x01);

		/* for each substituted kmer */
		for(size_t i = 0; i < ksize; i += kpitch) {
			uint32_t const ksub = krem + i;
			if(p[ksub].exist == 0) { continue; }

			/* hit */
			// debug("max_shift(%zu), shift(%zu), kmer(%x), ksub(%x), len(%zu)", max_shift, shift, kmer, ksub, (max_shift - shift) / 2);
			return((tm_ref_prefix_t){
				#ifdef DEBUG_KMER
					.kmer = ksub,
				#endif

				.len = (max_shift - shift) / 2,
				.ofs = p[ksub].ofs
			});
		}
	} while((shift -= 2) != 0);

	/* never reach here?? */
	return((tm_ref_prefix_t){
		.len = max_shift / 2,
		.ofs = TM_REF_BASE_OFS
	});
}

static _force_inline
uint64_t tm_ref_build_link(tm_ref_cnt_t const *p, size_t ksize, tm_ref_sketch_t *sk)
{
	size_t const max_shift = _tzc_u64(ksize);
	uint64_t const mask = ksize - 1;

	for(size_t i = 0; i < ksize; i++) {
		if(p[i].exist == 0) { continue; }

		tm_ref_bin_t *bin = _add_offset(sk, p[i].ofs);
		#ifdef DEBUG_KMER
			bin->kmer = i;
		#endif

		/* find next link for each base */
		for(size_t j = 0; j < 4; j++) {
			tm_ref_prefix_t n = tm_ref_find_link(p, ((i<<2) + j) & mask, max_shift);
			int32_t const next = (int32_t)(n.ofs - p[i].ofs) / (int32_t)TM_REF_ALIGN_SIZE;	/* overflows; keep sign bit */
			if(_unlikely(next > INT16_MAX || next < INT16_MIN)) {
				// error("link overflow at k-mer(%zx) and next base(%c), try smaller k-mer size for this sequence or shorten the sequence.", i, "ACGT"[j]);
				return(1);
			}

			/* save link; preserve tail */
			bin->link[j].plen = n.len;
			bin->link[j].next = next;
		}
	}
	return(0);
}

static _force_inline
tm_ref_sketch_t *tm_ref_build_index(tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen, tm_ref_kpos_t const *kpos, size_t kcnt, tm_ref_cnt_v *buf)
{
	_unused(seq);
	kv_clear(*buf);

	/* reserve profile buffer */
	size_t const ksize = 0x01ULL<<conf->kbits;
	tm_ref_cnt_t *p = kv_reserve(tm_ref_cnt_t, *buf, ksize);
	memset(p, 0, ksize * sizeof(tm_ref_cnt_t));

	/* count kmers */
	tm_ref_count_kmer(p, ksize, kpos, kcnt);
	_scan_memory(p, sizeof(tm_ref_cnt_t) * ksize);		/* valgrind use-of-uninitialized-value */

	/* patch lead */
	tm_ref_patch_feeder(p, ksize);

	/* accumulate block size */
	size_t size = tm_ref_calc_size(p, ksize);
	if(size > UINT32_MAX) { trap(); }

	/* malloc; use entire page */
	size_t const rounded_size = _roundup(size + conf->margin.head + MAX2(16, conf->margin.tail), 4096);
	void *base = malloc(rounded_size);
	if(base == NULL) { return(NULL); }		/* wrapped malloc traps allocation failure */

	/* save metadata */
	tm_ref_sketch_t *sk = _add_offset(base, conf->margin.head);
	*sk = (tm_ref_sketch_t){
		.size  = size,
		.head_margin = conf->margin.head,
		.kbits = conf->kbits,
		.slen  = slen
	};

	/* pack pos array */
	tm_ref_pack_kpos(p, ksize, kpos, kcnt, sk);
	if(tm_ref_build_link(p, ksize, sk)) {
		/* failed */
		free(base);
		return(NULL);
	}
	_scan_memory(sk, sk->size);				/* valgrind use-of-uninitialized-value */

	/* anything else? */
	debug("done, size(%u), head_margin(%u), kbits(%u), slen(%u)", sk->size, sk->head_margin, sk->kbits, sk->slen);
	return(sk);
}


/* working buffer for index construction */
struct tm_ref_tbuf_s {
	/* we expect the following two members are cleard when the tm_ref_sketch is first called */
	tm_ref_kpos_v kpos;
	tm_ref_cnt_v prof;
};

// static _force_inline
tm_ref_tbuf_t *tm_ref_tbuf_init(void)
{
	tm_ref_tbuf_t *tbuf = malloc(sizeof(tm_ref_tbuf_t));
	memset(tbuf, 0, sizeof(tm_ref_tbuf_t));
	return(tbuf);
}

// static _force_inline
void tm_ref_tbuf_destroy(tm_ref_tbuf_t *tbuf)
{
	kv_destroy(tbuf->kpos);
	kv_destroy(tbuf->prof);
	return;
}

// static _force_inline
tm_ref_sketch_t *tm_ref_sketch(tm_ref_tbuf_t *self, tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen)
{
	if(slen > INT16_MAX) { return(NULL); }

	/* slice kmers from sequence (ambiguous bases are expanded here) */
	if(tm_ref_collect_kmer(conf, seq, slen, &self->kpos) == 0) {		/* with additional base */
		debug("failed");
		return(NULL);
	}

	/* build index */
	tm_ref_kpos_t *kpos = kv_ptr(self->kpos);
	size_t const kcnt = kv_cnt(self->kpos);
	return(tm_ref_build_index(conf, seq, slen, kpos, kcnt, &self->prof));
}


/**
 * end of dbg.c
 */
