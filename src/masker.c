
/**
 * @file masker.c
 */

#include <stdlib.h>


/* alphabets */
enum alphabet_iupac {
	A = 0x01,
	C = 0x02,
	G = 0x04,
	T = 0x08,
	/* IUPAC ambiguous bases */
	R = A | G,
	Y = C | T,
	S = G | C,
	W = A | T,
	K = G | T,
	M = A | C,
	/* invalid */
	N = 0
};

enum alphabet_2bit {
	nA = 0x00,
	nC = 0x01,
	nG = 0x02,
	nT = 0x03
};


/* index object, allowing matching of IUPAC ambiguous base */
typedef struct {
	uint32_t bcnt;				/* first-stage hash table size (we expect smaller than 1024) */
	uint32_t const *bofs;

	uint32_t const *pos;
} mk_idx_t;


typedef struct {
	uint32_t kmer, pos;
} mk_query_pos_t;



/* index data structure; we expect each repeat sequence is shorter than 32kbp */
typedef struct {
	uint8_t const *seq;
	uint32_t slen, plen;		/* sequence length and pos array length */
	char const *name;			/* structured name? */

	uint32_t bofs[];			/* bin offsets */
} mk_ref_sketch_t;


/*
 * reference bin layout:
 *
 * typedef struct {
 *     uint16_t krem[n];
 *     uint32_t jmp[n + 1];
 *     uint32_t pos[];
 * } mk_ref_bin_t;
 */
typedef struct {
	uint16_t cnt;
	uint16_t krem[1];		/* unaligned */
} mk_ref_bin_t;
_static_assert(sizeof(mk_ref_bin_t) == sizeof(uint32_t));

typedef struct {
	uint16_t cnt;
	uint16_t pos[1];		/* unaligned */
} mk_ref_subarr_t;
_static_assert(sizeof(mk_ref_subarr_t) == sizeof(uint32_t));

typedef struct {
	size_t idx;
	uint64_t found;
} mk_ref_match_rem_t;


/* match object */
typedef struct {
	uint16_t const *ptr;
	size_t cnt;
} mk_match_t;


static _force_inline
mk_ref_bin_t const *mk_ref_get_bin(mk_ctx_t const *self, mk_ref_sketch_t const *ref, uint64_t kbase)
{
	_unused(self);

	/* we expect bbase will be reused when this function is called multiple times; because bbase is constant among these calls */
	mk_ref_bin_t const *bbase = (mk_ref_bin_t const *)ref;
	return(&bbase[ref->bofs[kbase]]);		/* just add offset */
}

/* match k-mer remainder */
static _force_inline
mk_ref_match_rem_t mk_ref_match_rem(mk_ctx_t const *self, mk_ref_bin_t const *bin, uint64_t krem)
{
	#define _match_v16i16(a, b) ( ((v16_masku_t){ .mask = _mask_v16i16(_eq_v16i16(a, b)) }).all )

	/* broadcast kmer remainder */
	v16i16_t const vr = _set_v16i16(krem);

	/* bulk match; ~12clks / loop */
	size_t const bcnt = bin->cnt;
	for(size_t i = 0; i < bcnt; i += 16) {
		/* load and compare 16 elements */
		v16i16_t const v = _loadu_v16i16(&bin->krem[0]);
		uint32_t const mask = _match_v16i16(v, vr);

		/* return index if found */
		if(mask != 0) {
			ZCNT_RESULT size_t cnt = _tzcnt_u32(mask);
			return((mk_ref_match_rem_t){
				.idx = i + (cnt>>1),
				.found = 1
			});
		}
	}

	/* not found */
	return((mk_ref_match_rem_t){
		.idx = 0,			/* unused */
		.found = 0
	});

	#undef _match_v16i16
}

static _force_inline
mk_match_t mk_ref_slice_arr(mk_ctx_t const *self, mk_ref_bin_t const *bin, size_t idx)
{
	size_t bsize = sizeof(uint16_t) * (1 + bin->cnt);	/* including bcnt at the head */
	uint32_t const *jmp = _add_offset(bin, bsize);		/* base pointer */
	mk_ref_subarr_t const *arr = &jmp[jmp[idx]];		/* add offset (jmp[idx]) to base (jmp); target is 4byte aligned */

	/* size is saved at the head of the subarray; extract it */
	return((mk_match_t){
		.ptr = arr->pos,
		.cnt = arr->cnt
	});
}

static _force_inline
mk_match_t mk_ref_match(mk_ctx_t const *self, mk_ref_sketch_t const *ref, uint64_t kmer)
{
	/* divide kmer into base and remainder (for first and second stage matchings, respectively) */
	uint64_t const kbase = kmer & self->kmer.bmask;
	uint64_t const krem  = kmer>>self->kmer.shift;

	/* first-stage matching: naive hash table */
	mk_ref_bin_t const *bin = mk_ref_get_bin(self, ref, kbase);

	/* match second stage: linear probing */
	mk_ref_match_rem_t m = mk_ref_match_rem(self, bin, krem);
	if(m.found == 0) {
		/* not found */
		return((mk_match_t){ .ptr = NULL, .cnt = 0 });		/* we expect control paths are kept isolated */
	}

	/* hit! get subarray pointer and size */
	mk_match_t arr = mk_ref_slice_arr(self, bin, m.idx);	/* arr.ptr won't be NULL */

	/* anything else to do? */
	return(arr);
}

/* construct reference object */
mk_ref_sketch_t *mk_ref_sketch(mk_ctx_t const *self, char const *name, uint8_t const *seq, size_t slen)
{
	return(sk);
}


/* query sequence be shorter than 2Gbp */
typedef struct {
	char const *name;
	uint8_t const *seq;
	size_t slen;
	uint32_t const *kmer;
} mk_query_sketch_t;


/* seed */
typedef struct {
	uint32_t u, v;
} mk_seed_t;
_static_assert(sizeof(mk_seed_t) == 8);
typedef struct { mk_seed_t *a; size_t n, m; } mk_seed_v;


/* chain */
typedef struct {
	/* base seed position; converted to u-v coordinate */
	mk_seed_t pos;

	/* reference and query side spans */
	struct {
		uint16_t u, v;
	} span;

	/* attributes */
	union {
		struct {
			uint32_t rid : 24;		/* reference id */
			uint32_t weight : 8;	/* weight; calculated from seed count and sequence length */
		} stat;
		uint32_t scnt;				/* chained seed count */
	} attr;
} mk_chain_t;
_static_assert(sizeof(mk_chain_t) == 16);
typedef struct { mk_chain_t *a; size_t n, m; } mk_chain_v;


/* build k-mer array; we don't duplicate name and seq arrays; just copy pointers (shallow copy) */
static _force_inline
mk_query_sketch_t mk_query_sketch(mk_ctx_t *self, char const *name, uint8_t const *seq, size_t slen)
{
	/* allocate k-mer array */
	uint32_t *karr = kv_reserve(self->kmer.arr, _roundup(slen, 32) / 16 + 2);

	/* 4bit -> 2bit conversion table */
	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		[N] = 0x00,
		[A] = (nA<<6) | nA,
		[C] = (nC<<6) | nC,
		[G] = (nG<<6) | nG,
		[T] = (nT<<6) | nT
		/* FIXME: add iupac ambiguous bases */
	};
	v32i8_t const cv = _cvt_v16i8_v32i8(_load_v16i8(conv));		/* rip relative load with vbroadcasti128 */

	/* gather mask */
	static uint8_t const gather[16] __attribute__(( aligned(16) )) = {
		3, 7, 11, 15
	};
	v32i8_t const gv = _cvt_v16i8_v32i8(_load_v16i8(gather));

	/* process 32bytes at once; 5~9clks / loop */
	for(size_t i = 0; i < slen; i += 32) {
		v32i8_t const v = _loadu_v32i8(&seq[i]);
		v32i8_t const w = _shuf_v32i8(cv, v);

		/* gather 2bit cols */
		v32i8_t const g0 = _or_v32i8(
			_shl_4i32(w, 6),
			_shl_4i32(w, 18)
		);

		/* gather bytes */
		v32i8_t const g1 = _shuf_v32i8(g0, gv);
		karr[(i / 16)    ] = (uint32_t)_ext_v4i32(g1, 0);
		karr[(i / 16) + 1] = (uint32_t)_ext_v4i32(g1, 4);
	}
	return((mk_query_sketch_t){
		.name = name,
		.seq  = seq,
		.slen = slen,
		.kmer = karr
	});
}

/* copy seed chunk from reference index object converting coordinates */
static _force_inline
size_t mk_expand_seed(mk_ctx_t const *self, uint16_t const *ptr, size_t cnt, uint32_t pos, mk_seed_v *seed)
{
	/* reserve space for this batch */
	size_t scnt = kv_size(seed);
	mk_seed_t *q = kv_reserve(mk_seed_t, seed, scnt + cnt + 8);

	/* compose x-coordinate vector */
	v8i32_t const x1 = _set_v8i32(-pos);
	v8i32_t const x2 = _set_v8i32(2 * pos);

	/* bulk */
	for(size_t i = 0; i < cnt; i += 8) {
		v8i16_t const w = _loadu_v8i16(&ptr[i]);
		v8i32_t const y = _cvt_v8i16_v8i32(w);

		/* (x, y) -> (2x - y, 2y - x) */
		v8i32_t const u = _sub_v8i32(x2, y);
		v8i32_t const v = _add_v8i32(_add_v8i32(y, y), x1);

		/* fold in */
		v8i32_t const l = _lo_v8i32(u, v);
		v8i32_t const h = _hi_v8i32(u, v);

		/* save */
		_storeu_v8i32(&q[i],     l);
		_storeu_v8i32(&q[i + 4], h);
	}

	kv_cnt(seed) += cnt;
	return(cnt);
}

static _force_inline
size_t mk_collect_seed(mk_ctx_t const *self, mk_ref_sketch_t const *ref, mk_query_sketch_t const *query, mk_seed_v *seed)
{
	uint64_t const kmask = self->kmer.amask;

	for(size_t i = 0; i < query->kcnt; i += 16) {
		uint64_t kmer = _loadu_u64(&query->kmer[i / 16]);

		for(size_t j = 0; j < MIN2(16, query->kcnt - i); j++) {
			mk_match_t m = mk_ref_match(self, ref, kmer & kmask);
			kmer >>= 2;

			/* not found? */
			if(m.ptr == NULL) { continue; }

			/* found; convert matches to seeds */
			mk_expand_seed(self, m.ptr, m.cnt, i + j, seed);
		}
	}
	return(kv_cnt(seed));
}

static _force_inline
int64_t mk_chain_test_ptr(mk_seed_t const *p, mk_seed_t const *t)
{
	return((int64_t)((ptrdiff_t)(t - p - 1)));
}

static _force_inline
mk_seed_t *mk_chain_find_first(mk_ctx_t const *self, mk_seed_t *p, mk_seed_t *t, uint64_t lb)
{
	/* constants */
	uint64_t const window = self->chain.window.all;	/* window sizes */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const tmask = 0x8000000080000000;	/* extract two sign bit pair */

	/* load bounds */
	uint64_t const uv = _loadu_u64(p);			/* (ulb, vlb) */
	uint64_t ub = (uv & umask) + window;		/* (uub, vub - vlb) */

	int64_t cont = 0;		/* continuous flag */
	while(1) {
		/* check if reached tail */
		if((cont | mk_chain_test_ptr(++p, t)) < 0) { return(NULL); }

		uint64_t v = _loadu_u64(p) - lb;		/* 1: (u, v - vlb) for inclusion test */
		if(((ub - v) & tmask) == 0) { break; }	/* 2,3: break if chainable (first chainable seed found) */

		/* unchainable; test if out of uub */
		cont = ub - v;		/* save diff for testing MSb */
	}
	return(p);
}

static _force_inline
mk_seed_t *mk_chain_find_alt(mk_ctx_t const *self, mk_seed_t *p, mk_seed_t *t, uint64_t lb)
{
	/* constants */
	uint64_t const tmask = 0x8000000080000000;	/* extract two sign bit pair */

	/* calculate bounds */
	uint64_t rv = _loadu_u64(p) - lb;
	uint64_t ub = rv + (rv<<32);

	/* keep nearest seed */
	mk_seed_t *n = p;

	int64_t cont = 0;
	while((cont | mk_chain_test_ptr(++p, t)) >= 0) {
		uint64_t v = _loadu_u64(p) - lb;		/* 1: (u, v - vlb) for inclusion test */
		cont = ub - v;
		if((ub - v) & tmask) { continue; }		/* skip if unchainable */

		/* chainable; test if the seed is nearer than the previous */
		uint64_t w = v + (v<<32);				/* 2,3: (u + v - vlb, v - vlb) */
		if((ub - w) & tmask) { break; }			/* 4,5: further than previous */

		/* nearer seed found */
		ub = w;			/* update bounds */
		n = p;			/* save pointer */
	}
	return(n);
}

static _force_inline
mk_seed_t *mk_chain_find_nearest(mk_ctx_t const *self, mk_seed_t *p, mk_seed_t *t)
{
	/* load root positions */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const uv = _loadu_u64(p);			/* (ulb, vlb) */
	uint64_t lb = (uv & ~umask);				/* (0, vlb) */

	/* find first chainable seed */
	mk_seed_t *n = mk_chain_find_first(p, t, lb);
	if(n == NULL) { return(NULL); }				/* we expect control paths are kept isolated */

	/* then inspect alternative chainable seeds */
	return(mk_chain_find_alt(self, n, t, lb));
}

static _force_inline
size_t mk_chain_record(mk_ctx_t const *self, v2i32_t root, v2i32_t tail, v2i32_t scnt, mk_seed_t const *t, mk_chain_t *q)
{
	/* avoid use of general-purpose register to minimize spills */
	v2i32_t const span = _sub_v2i32(tail, root);

	/* FIXME: raw SSE operations */
	v4i32_t const chain = ({
		__m128i const x = _mm_shufflelo_epi16(span.v1, 0xd4);	/* (uint32_t [4]){ -, -, -, (vspan, uspan) } */
		__m128i const y = _mm_blend_epi16(x, scnt.v1, 0x0c);	/* (uint32_t [4]){ -, -, scnt, (vspan, uspan) } */
		__m128i const z = _mm_unpacklo_epi64(root.v1, y);		/* (uint32_t [4]){ scnt, (vspan, uspan), vpos, upos } */
		((v4i32_t){ .v1 = z });
	});

	_storeu_v4i32(q, chain);
	return(1);
}

static _force_inline
size_t mk_chain_seed(mk_ctx_t const *self, mk_seed_t const *seed, size_t scnt, mk_chain_v *chain)
{
	/*
	 * greedy extension; relatively small window size (~32) is preferred for good performance.
	 * do all of the chainablity test on general-purpose register for smaller latency of core loop.
	 * keep seed and counter on xmm regsters and do everything of the chain composing on xmm registers
	 * to minimize spill of general-purpose registers.
	 */
	uint64_t const chained = self->chain.chained;	/* offset */
	v2i32_t const inc = _set_v2i32(1, 0);

	/* src pointers */
	mk_seed_t *p = seed - 1, *t = &seed[scnt];

	/* reserve mem for chain */
	mk_chain_t *q = kv_reserve(mk_chain_t, *chain, kv_cnt(*chain) + scnt);

	/* for each seed */
	while(++p < t) {
		/* skip chained seed to find next root */
		if(p->u >= chained) { continue; }

		/* root found; keep seed on xmm register */
		v2i32_t const root = _loadu_v2i32(p);		/* we expect this won't be spilled */
		v2i32_t scnt = _zero_v2i32();

		/* iteratively link nearest seed */
		mk_seed_t *s = p;
		while(1) {
			mk_seed_t *n = mk_chain_find_nearest(self, s, t);
			if(n == NULL) { break; }

			/* increment seed count */
			scnt = _add_v2i32(scnt, inc);

			n->u += chained;	/* mark chained */
			s = n;				/* save last chained seed */
		}

		q += mk_chain_record(self, root, _loadu_v2i32(s), scnt, q);
	}

	/* update chain count */
	size_t const ccnt = q - &kv_ptr(*chain)[kv_cnt(*chain)];
	kv_cnt(*chain) = ccnt;
	return(ccnt);
}

static _force_inline
v16i8_t mk_load_forward(mk_ctx_t const *self, uint8_t const *ptr)
{
	return(_loadu_v16i8(ptr));
}

static _force_inline
v16i8_t mk_load_reverse(mk_ctx_t const *self, uint8_t const *ptr)
{
	v16i8_t const mask = _load_v16i8(self->seq.smask);
	v16i8_t const v = _loadu_v16i8(ptr);

	/* flip */
	return(_shuf_v16i8(mask, v));
}

typedef struct {
	uint8_t const *ptr;
	v16i8_t (*load)(mk_ctx_t const *, uint8_t const *);
	int64_t fwd;
} mk_filter_t;

static _force_inline
int64_t mk_filter_extend(mk_ctx_t const *self, mk_filter_t r, mk_filter_t q)
{
	/* load sequence vectors */
	v16i8_t rv = r.load(self, r.ptr), rn = r.load(self, r.ptr + r.fwd);
	v16i8_t qv = q.load(self, q.ptr), qn = q.load(self, q.ptr + q.fwd);

	/* load constants */
	v16i8_t const score_matrix = _load_v16i8(&self->filter.score_matrix);
	v16i8_t dh = _load_v16i8(&self->filter.init);
	v16i8_t dv = _load_v16i8(&self->filter.init);
	v16i8_t sv = _zero_v16i8(), mv = _zero_v16i8();

	#define _calc_a(_rv, _qv, _dh, _dv) ({ \
		/* calc score profile vector with shuffling score matrix vector */ \
		v16i8_t const match = _and_v16i8(_rv, _qv);		/* cmpeq */ \
		v16i8_t const a = _max_v16i8( \
			_shuf_v16i8(score_matrix, match), \
			_max_v16i8(_dv, _dh) \
		); \
		a; \
	})

	/* 32 vector updates */
	for(size_t j = 0; j < 16; j++) {
		/* #0 (even) */
		v16i8_t const a0 = _calc_a(rv, qv, dh, dv);
		v16i8_t const tdv = _sub_v16i8(a0, dh);
		v16i8_t const tdh = _sub_v16i8(a0, dv);
		sv = _add_v16i8(sv, tdv);
		mv = _max_v16i8(mv, sv);

		/* move rightward */
		rv = _shrd_v16i8(rv, rn, 1); rn = _shr_v16i8(rn, 1);

		/* #1 (odd) */
		v16i8_t const a1 = _calc_a(rv, qv, tdh, tdv);
		dv = _sub_v16i8(a1, tdh);
		dh = _sub_v16i8(a1, tdv);
		sv = _add_v16i8(sv, dh);
		mv = _max_v16i8(mv, sv);

		/* move downward */
		qv = _shrd_v16i8(qv, qn, 1); qn = _shr_v16i8(qn, 1);
	}

	/* fold scores */
	return(score);

	#undef _calc_a
}

static _force_inline
int64_t mk_filter_extend_inward(mk_ctx_t const *self, mk_ref_sketch_t const *ref, mk_query_sketch_t const *query, mk_chain_t const *chain)
{
	return(0);
}

static _force_inline
int64_t mk_filter_extend_outward(mk_ctx_t const *self, mk_ref_sketch_t const *ref, mk_query_sketch_t const *query, mk_chain_t const *chain)
{
	return(0);
}

static _force_inline
size_t mk_filter_save_chain(mk_ctx_t const *self, uint32_t rid, mk_chain_t *q, mk_chain_t const *p)
{
	/* load */
	v16i8_t v = _loadu_v16i8(p);

	/* calc weight */
	ZCNT_RESULT size_t zc = _lzcnt_u32(p->attr.ccnt);
	uint32_t const weight = 32 - zc;

	/* save */
	_storeu_v16i8(q, v);
	q->attr.stat.rid = rid;
	q->attr.stat.weight = weight;
	return(1);
}

static _force_inline
size_t mk_filter_chain(mk_ctx_t const *self, mk_ref_sketch_t const *ref, mk_query_sketch_t const *query, uint32_t rid, mk_chain_t *chain, size_t ccnt)
{
	/* alias pointer */
	mk_chain_t const *src = chain;
	mk_chain_t *dst = chain;

	/* load constants */
	uint64_t const min_cnt = self->filter.min_cnt;
	uint64_t const test_cnt = self->filter.test_cnt;
	int64_t const min_score = self->filter.min_score;
	size_t const ulim = self->filter.ulim;

	for(size_t i = 0; i < ccnt; i++) {
		mk_chain_t const *p = &src[i];

		/* skip if too small */
		if(p->cnt < min_cnt) { continue; }

		/* try short extension if not heavy enough */
		if(p->cnt < test_cnt) {
			int64_t score = (p->uspan > ulim
				? mk_filter_extend_inward(self, ref, query, p)		/* inward extension for long-spanning chain */
				: mk_filter_extend_outward(self, ref, query, p)
			);
			if(score < st) { continue; }
		}

		/* copy and fold in reference id */
		dst += mk_filter_save_chain(self, rid, dst, p);
	}
	return(ccnt);
}

static _force_inline
size_t mk_collect_candidate(mk_ctx_t *self, mk_ref_sketch_t const *ref, mk_query_sketch_t const *query, mk_chain_v *chain)
{
	/* enumerate seeds */
	kv_clear(self->seed);
	if(mk_collect_seed(self, ref, query, &self->seed) == 0) {
		return(0);
	}

	/* sort seeds by u-coordinates for chaining */
	radix_sort_seed(kv_ptr(self->seed), kv_cnt(self->seed));

	/* chaining; return if no chain obtained */
	size_t const cbase = kv_cnt(chain);		/* save chain count */
	size_t const ccnt  = mk_chain_seed(self, kv_ptr(self->seed), kv_cnt(self->seed), chain);
	if(ccnt == 0) { return(0); }

	/* filter (try small extension with simple score matrix) */
	size_t const cfilt = mk_filter_chain(self, ref, query, &kv_ptr(chain)[cbase], ccnt);
	if(cfilt == 0) { return(0); }

	/* finalize filtered chains */
	kv_cnt(chain) = cbase + cfilt; 
	return(kv_cnt(chain));
}

/* X-drop DP extension */
static _force_inline
mk_aln_t mk_extend_seed(mk_ctx_t const *self, mk_ref_sketch_t const *ref, mk_query_sketch_t const *query, mk_chain_t const *chain)
{
	/* upward then downward */

	return;
}

/* record-test pair */
static _force_inline
uint64_t mk_bin_is_covered(mk_bin_t *bin, mk_chain_t const *chain)
{
	/* for each overlapping alignment */

	return;
}

static _force_inline
void mk_bin_record(mk_bin_t *bin, mk_chain_t const *chain, mk_aln_t const *aln)
{
	return;
}

static _force_inline
void mk_bin_clear(mk_bin_t *bin)
{
	return;
}

static _force_inline
size_t mk_bin_cnt(mk_bin_t *bin)
{
	return(0);
}

static _force_inline
size_t mk_extend(mk_ctx_t *self, mk_ref_sketch_t const *ref, mk_query_sketch_t const *query, mk_chain_t const *chain, size_t ccnt)
{
	mk_bin_clear(&self->bin);

	dz_query_t *q = dz_pack_query(self->extend.dz, query->seq, query->slen);

	/* for each chain try extension */
	for(size_t i = 0; i < ccnt; i++) {
		mk_chain_t const *p = &chain[i];

		/* get reference sequence */
		mk_ref_sketch_t const *r = &ref[p->rid];

		/* skip if already covered (and there is no possibility that this chain surpasses the previous ones) */
		if(mk_bin_is_covered(&self->bin, p)) { continue; }

		mk_aln_t aln = mk_extend_seed(self, r, query, p);
		if(aln.score <= 0) { continue; }

		mk_record(&self->bin, p, &aln);
	}


	;
	return(mk_bin_cnt(&self->bin));
}

static _force_inline
size_t mk_scan(mk_ctx_t *self, mk_query_sketch_t const *query)
{
	/* clear chains */
	kv_clear(self->chain);

	/* collect chain for each reference sequence */
	mk_ref_sketch_t const *ref = self->ref.ptr;
	for(size_t i = 0; i < self->ref.cnt; i++) {
		mk_collect_candidate(self, &ref[i], query, &self->chain);
	}

	if(kv_cnt(self->chain) == 0) { return(0); }

	/* sort by estimated identity */
	radix_sort_chain(kv_ptr(self->chain), kv_cnt(self->chain));

	/* convert chain to alignment */
	size_t acnt = mk_extend(self, ref, query, kv_ptr(self->chain), kv_cnt(self->chain));
	return(acnt);
}





/**
 * end of masker.c
 */
