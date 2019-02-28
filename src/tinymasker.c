
/**
 * @file tinymasker.c
 */

#include <stdlib.h>

#ifndef UNITTEST
#  define UNITTEST				( 0 )
#endif

#define UNITTEST_UNIQUE_ID		1
#include "utils/utils.h"		/* include all */

#define DZ_PROTEIN
#define DZ_MAT_SIZE				( 16 )
#include "dozeu.h"


/* misc */
#ifndef MK_VERSION
#  define MK_VERSION			"tinymasker-0.0.1"
#endif
#define MK_ARCH_NAME			ARCH_NAME			/* SSE4.1 or AVX2 */

/* version string */
static _force_inline
char const *mk_version(void)
{
	char const *prefix = "tinymasker-";

	/* remove prefix */
	uint64_t spos = mm_startswith(MK_VERSION, prefix) ? strlen(prefix) : 0;
	return(&MK_VERSION[spos]);
}


unittest_config( .name = "tinymasker" );


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


/* index data structure; we expect each repeat sequence is shorter than 32kbp */
typedef struct {
	size_t size;		/* object size and sequence length */
	uint32_t kbits, slen;
	uint8_t const *seq;
	char const *name;			/* structured name? */
} mk_ref_sketch_t;

/* construct reference object */
typedef struct {
	uint32_t kmer, pos;
} mk_ref_kpos_t;
_static_assert(sizeof(mk_ref_kpos_t) == 8);

typedef struct { mk_ref_kpos_t *a; size_t n, m; } mk_ref_kpos_v;

#define mk_ref_kmer(x)			( (x).kmer )
KRADIX_SORT_INIT(kmer, mk_ref_kpos_t, mk_ref_kmer, 4);

typedef struct {
	uint32_t f, r;
} mk_ref_kmer_pair_t;

typedef struct {
	/* constants */
	uint32_t amask, shift;

	/* dst array */
	mk_ref_kpos_t *q;

	/* kmer stack */
	uint32_t size;
	uint32_t branches;
	mk_ref_kmer_pair_t kmer[256];	/* branch stack; max ambiguity = 4^4 */
} mk_ref_work_t;

static _force_inline
uint64_t mk_ref_base_to_2bit(uint8_t base)
{
	/* convert iupac nucleotide to pair of 2-bit encoded bases */
	uint64_t const magic = 0x498e4dc399824100;
	return((magic>>(4 * base)) & 0x0f);
}

static _force_inline
uint64_t mk_ref_is_ambiguous(uint8_t base)
{
	uint64_t const magic = 0x1111111011101000;
	return((magic>>(4 * base)) & 0x01);
}

static _force_inline
size_t mk_ref_dup_stack(mk_ref_work_t *w)
{
	size_t prev_size = w->size;

	/* update size and branch pos */
	w->size = w->size<<1;
	w->branches = (w->branches>>2) | (0x01<<w->shift);

	for(size_t i = 0; i < prev_size; i++) {
		w->kmer[w->size + i] = w->kmer[i];
	}
	return(prev_size);
}

static _force_inline
void mk_ref_shrink_stack(mk_ref_work_t *w)
{
	w->size = w->size>>1;
	for(size_t i = 0; i < w->size; i++) {
		w->kmer[i] = w->kmer[2 * i];
	}
	return;
}

static _force_inline
void mk_ref_update_kmer(mk_ref_kmer_pair_t *kmer, size_t size, uint32_t amask, uint32_t shift, uint8_t base)
{
	for(size_t i = 0; i < size; i++) {
		kmer[i].f = ((kmer[i].f<<2) | base) & amask;
		kmer[i].r =  (kmer[i].r>>2) | ((0x03 ^ base)<<shift);
	}
	return;
}

static _force_inline
void mk_ref_push_pos(mk_ref_work_t *w, uint32_t pos)
{
	for(size_t i = 0; i < w->size; i++) {
		*w->q++ = (mk_ref_kpos_t){
			.kmer = w->kmer[i].f,
			.pos = pos
		};
		*w->q++ = (mk_ref_kpos_t){
			.kmer = w->kmer[i].r,
			.pos = pos			/* FIXME: direction flag needed? */
		};
	}
	return;
}

static _force_inline
void mk_ref_push_base(mk_ref_work_t *w, uint8_t base)
{
	size_t const size = w->size;		/* save current size before duplicating stack */
	uint8_t const b2 = mk_ref_base_to_2bit(base);

	if(_unlikely(mk_ref_is_ambiguous(base))) {
		mk_ref_dup_stack(w);
		mk_ref_update_kmer(&w->kmer[size], size, w->amask, w->shift, b2>>2);
	}

	/* push base */
	mk_ref_update_kmer(&w->kmer[size], size, w->amask, w->shift, b2 & 0x03);
	return;
}

static _force_inline
mk_ref_kpos_t *mk_ref_reserve(mk_ref_kpos_v *buf, mk_ref_kpos_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	mk_ref_kpos_t *p = kv_reserve(mk_ref_kpos_t, *buf, 1024);
	return(&p[kv_cnt(*buf)]);
}

static _force_inline
size_t mk_ref_collect_kmer(size_t kbits, uint8_t const *seq, size_t slen, mk_ref_kpos_v *buf)
{
	mk_ref_work_t w = {
		.amask = (0x01<<kbits) - 0x01,
		.shift = kbits - 2,

		/* make room in dst array */
		.q = mk_ref_reserve(buf, kv_ptr(*buf)),

		/* clear stack; keeps all the combination of ambiguous bases */
		.size = 1,
		.branches = 0,
		.kmer = { { 0 } }
	};

	/* push bases at the head */
	size_t const k = kbits>>1;
	for(size_t i = 0; i < k - 1; i++) {
		mk_ref_push_base(&w, seq[i]);
	}

	/* body */
	for(size_t i = k - 1; i < slen; i++) {
		mk_ref_push_base(&w, seq[i]);			/* stack will be expanded if needed */
		mk_ref_push_pos(&w, i);

		if(_unlikely(w.branches & 0x01)) {
			mk_ref_shrink_stack(&w);			/* shrink if needed */
		}

		/*
		 * check need for expansion every 32times:
		 * Skylake's branch predictor successfully learns once-every-32times event, but fails for more
		 */
		if((i & 0x1f) != 0) { continue; }
		w.q = mk_ref_reserve(buf, w.q);
	}

	/* done; record count */
	kv_cnt(*buf) = w.q - kv_ptr(*buf);
	return(kv_cnt(*buf));
}


/* rotate */
static _force_inline
size_t mk_ref_rotate_kmer(size_t kbits, mk_ref_kpos_t *kpos, size_t kcnt)
{
	uint32_t const mask = (0x01<<kbits) - 0x01;

	/* automatically vectorized? */
	for(size_t i = 0; i < kcnt; i++) {
		uint32_t const kmer = kpos[i].kmer;
		kpos[i].kmer = ((kmer & mask)<<2) | (kmer>>kbits);
	}
	return(kcnt);
}


/* pack kmer-position array */
typedef struct {
	uint32_t tail : 12;			/* pos count */
	uint32_t plen : 4;			/* unmatching prefix length */
	int32_t next  : 16;			/* diff to next bin */
} mk_ref_link_t;

typedef struct {
	mk_ref_link_t link[4];
	uint16_t pos[];
} mk_ref_bin_t;


_static_assert(sizeof(mk_ref_bin_t) == 16);		/* smallest bin size == 16 */

#define MK_REF_ALIGN_SIZE		( sizeof(uint64_t) )
#define MK_REF_BASE_OFS			( sizeof(mk_ref_sketch_t) )


typedef struct {
	uint32_t cnt, ofs;
} mk_ref_profile_t;
typedef struct { mk_ref_profile_t *a; size_t n, m; } mk_ref_profile_v;

static _force_inline
size_t mk_ref_build_profile(mk_ref_profile_t *p, size_t ksize, mk_ref_kpos_t const *kpos, size_t kcnt)
{
	_unused(ksize);

	for(mk_ref_kpos_t const *k = kpos, *t = &kpos[kcnt]; k < t; k++) {
		p[k->kmer>>2].cnt++;
	}
	return(kcnt);
}

static _force_inline
size_t mk_ref_calc_size(mk_ref_profile_t const *p, size_t ksize)
{
	size_t const hdr = sizeof(mk_ref_bin_t);

	size_t acc = 0;
	for(size_t i = 0; i < ksize; i++) {
		if(p[i].cnt == 0) { continue; }
		acc += _roundup(hdr + sizeof(uint16_t) * p[i].cnt, MK_REF_ALIGN_SIZE);
	}
	return(acc);
}

static _force_inline
size_t mk_ref_pack_kpos(mk_ref_profile_t *p, size_t ksize, mk_ref_kpos_t const *kpos, size_t kcnt, mk_ref_sketch_t *sk)
{
	_unused(ksize);

	/* first bin */
	size_t const hdr = sizeof(mk_ref_bin_t);
	size_t ofs = MK_REF_BASE_OFS;

	mk_ref_kpos_t const *k = kpos, *t = &kpos[kcnt];
	while(k < t) {
		/* save current offset */
		uint32_t const kbase = k->kmer>>2;
		p[kbase].ofs = ofs;

		/* slice bin from current offset */
		mk_ref_bin_t *bin = _add_offset(sk, ofs);
		_storeu_v4i32(bin, _zero_v4i32());		/* clear link */

		/* pack pos array */
		size_t i = 0, j = 0;
		do {
			uint32_t const kall = k->kmer;
			do {
				bin->pos[i++] = k++->pos;
			} while(k < t && k->kmer == kall);

			do {
				bin->link[j++].tail = i;
			} while(j < (kall & 0x03));
		} while(k < t && (k->kmer>>2) == kbase);

		/* update offset */
		ofs += _roundup(hdr + sizeof(uint16_t) * p[kbase].cnt, MK_REF_ALIGN_SIZE);
	}
	return(ofs);
}

typedef struct {
	uint32_t len, ofs;
} mk_ref_prefix_t;

static _force_inline
mk_ref_prefix_t mk_ref_find_link(mk_ref_profile_t const *p, uint32_t kmer, size_t max_depth)
{
	/* breadth-first search on kmer branching tree over profile table */
	size_t depth = 0;
	do {
		size_t const slen = 2 * depth;		/* substitution length (zero for the first iteration) */
		uint32_t const krem = kmer & (0 - (0x01<<slen));

		/* for each substituted kmer */
		for(size_t i = 0; i < 0x01ULL<<slen; i++) {
			uint32_t const ksub = krem + i;
			if(p[ksub].cnt == 0) { continue; }

			/* hit */
			return((mk_ref_prefix_t){
				.len = depth,
				.ofs = p[ksub].ofs
			});
		}
	} while(++depth < max_depth);

	/* never reach here?? */
	return((mk_ref_prefix_t){
		.len = max_depth,
		.ofs = MK_REF_BASE_OFS
	});
}

static _force_inline
void mk_ref_build_link(mk_ref_profile_t const *p, size_t ksize, mk_ref_sketch_t *sk)
{
	size_t const shift = _tzcnt_u64(ksize);
	for(size_t i = 0; i < ksize; i++) {
		if(p[i].cnt == 0) { continue; }

		mk_ref_bin_t *bin = _add_offset(sk, p[i].ofs);

		/* find next link for each base */
		for(size_t j = 0; j < 4; j++) {
			mk_ref_prefix_t n = mk_ref_find_link(p, (i + (j<<shift))>>2, shift>>1);

			/* save link; preserve tail */
			bin->link[j].plen = n.len;
			bin->link[j].next = n.ofs - p[i].ofs;	/* overflows; keep sign bit */
		}
	}
	return;
}

static _force_inline
mk_ref_sketch_t *mk_ref_build_index(size_t kbits, mk_ref_kpos_t const *kpos, size_t kcnt, mk_ref_profile_v *buf)
{
	kv_clear(*buf);

	/* reserve profile buffer */
	size_t const ksize = 0x01ULL<<kbits;
	mk_ref_profile_t *p = kv_reserve(mk_ref_profile_t, *buf, ksize);

	/* count kmers */
	mk_ref_build_profile(p, ksize, kpos, kcnt);

	/* accumulate block size */
	size_t size = mk_ref_calc_size(p, ksize);
	if(size > UINT32_MAX) { trap(); }

	mk_ref_sketch_t *sk = (mk_ref_sketch_t *)malloc(_roundup(size, 4096));	/* use entire page */
	if(sk == NULL) { return(NULL); }
	sk->kbits = kbits;
	sk->size = size;

	/* pack pos array */
	mk_ref_pack_kpos(p, ksize, kpos, kcnt, sk);
	mk_ref_build_link(p, ksize, sk);

	/* anything else? */
	return(sk);
}


/* working buffer */
typedef struct {
	mk_ref_kpos_v kpos;
	mk_ref_profile_v prof;
} mk_ref_t;

static _force_inline
mk_ref_sketch_t *mk_ref_sketch(mk_ref_t *self, size_t kbits, char const *name, uint8_t const *seq, size_t slen)
{
	/* slice kmers from sequence (ambiguous bases are expanded here) */
	if(mk_ref_collect_kmer(kbits + 2, seq, slen, &self->kpos) == 0) {		/* with additional base */
		return(NULL);
	}

	mk_ref_kpos_t *kpos = kv_ptr(self->kpos);
	size_t const kcnt = kv_cnt(self->kpos);

	/* rotate by a base before sort */
	mk_ref_rotate_kmer(kbits, kpos, kcnt);

	/* sort kmers by position */
	radix_sort_kmer(kpos, kcnt);

	/* build index */
	mk_ref_sketch_t *sk = mk_ref_build_index(kbits, kpos, kcnt, &self->prof);
	if(sk == NULL) { return(NULL); }

	/* save metadata (we do not copy them) */
	sk->seq  = seq;
	sk->slen = slen;
	sk->name = name;
	return(sk);
}


/* matcher */
typedef struct {
	size_t src, dst, cnt;
} mk_ref_squash_t;

typedef struct {
	size_t prefix, unmatching;
	mk_ref_bin_t const *bin;
} mk_ref_state_t;

typedef struct {
	mk_ref_state_t state;
	mk_ref_squash_t squash;
} mk_ref_next_t;

typedef struct {
	uint16_t const *ptr;
	size_t cnt;
} mk_ref_match_t;


static _force_inline
mk_ref_state_t mk_ref_match_init(mk_ref_sketch_t const *ref)
{
	return((mk_ref_state_t){
		.prefix     = ref->kbits>>1,
		.unmatching = 1,
		.bin        = _add_offset(ref, MK_REF_BASE_OFS)	/* bin for AAA... */
	});
}

static _force_inline
mk_ref_next_t mk_ref_match_next(mk_ref_state_t s, uint8_t next)
{
	int64_t const pitch = MK_REF_ALIGN_SIZE;

	/* get link to the next bin */
	mk_ref_link_t const link = s.bin->link[next];
	mk_ref_bin_t const *bin = _add_offset(s.bin, pitch * link.next);

	/* update matching state */
	size_t const prefix = s.prefix - s.unmatching + link.plen;
	size_t const sidx = bin->link[next - 1ULL].tail;
	return((mk_ref_next_t){
		.state = {
			.prefix     = prefix,
			.unmatching = prefix != 0,
			.bin        = bin
		},

		/* extract squashable subbin */
		.squash = {
			.src = next == 0 ? 0 : sidx,
			.dst = bin->link[next].tail,
			.cnt = bin->link[3].tail
		}
	});
}

static _force_inline
mk_ref_match_t mk_ref_get_arr(mk_ref_state_t s)
{
	return((mk_ref_match_t){
		.ptr = s.bin->pos,
		.cnt = s.bin->link[3].tail
	});
}


/* index builder */

typedef struct {
	uint32_t dummy;;
} mk_idx_conf_t;

typedef struct {
	uint32_t dummy;
} mk_idx_t;

static _force_inline
void mk_idx_destroy(mk_idx_t *mi)
{
	_unused(mi);

	return;
}

static _force_inline
mk_idx_t *mk_idx_gen(mk_idx_conf_t const *conf, pt_t *pt, char const *filename, FILE *log)
{
	_unused(conf);
	_unused(pt);
	_unused(filename);
	_unused(log);

	return(NULL);
}


/* index I/O */
static _force_inline
size_t mk_idx_dump(mk_idx_t const *idx, pg_t *pg, write_t write)
{
	_unused(idx);
	_unused(pg);
	_unused(write);

	size_t const size = 0;
	return(size);
}

static _force_inline
mk_idx_t *mk_idx_load(pg_t *pg, read_t read)
{
	_unused(pg);
	_unused(read);

	mk_idx_t *mi = NULL;
	return(mi);
}


/* seed */
typedef struct {
	uint32_t u, v;
} mk_seed_t;
_static_assert(sizeof(mk_seed_t) == 8);
typedef struct { mk_seed_t *a; size_t n, m; } mk_seed_v;

#define mk_seed_upos(x)			( (x).u )
KRADIX_SORT_INIT(seed, mk_seed_t, mk_seed_upos, 4);


typedef struct {
	uint16_t src, dst, cnt;
} mk_intv_t;
_static_assert(sizeof(mk_intv_t) == 6);
typedef struct { mk_intv_t *a; size_t n, m; } mk_intv_v;

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

#define mk_chain_attr(x)		( (x).attr.scnt )
KRADIX_SORT_INIT(chain, mk_chain_t, mk_chain_attr, 4);


/* alignment result */
typedef struct {
	uint32_t r, q;
} mk_pair_t;

typedef struct {
	/* rbt */
	rbt_header_t h;
	uint32_t qmax;

	/* stats */
	uint32_t min_weight;
	int64_t score;
	dz_alignment_t *aln;

	/* positions */
	mk_pair_t pos, span;

	/* alignment path */
	struct {
		uint8_t const *ptr;
		size_t len;
	} path;
} mk_aln_t;
_static_assert(sizeof(mk_aln_t) == 64);
typedef struct { mk_aln_t *a; size_t n, m; } mk_aln_v;

#define mk_aln_rbt_header(a)		( &(a)->h )
#define mk_aln_rbt_cmp(a, b)		( (a)->pos.q < (b)->pos.q )
#define mk_aln_ivt_cmp_head(a, b)	( (a)->qmax                > (b)->pos.q )
#define mk_aln_ivt_cmp_tail(a, b)	( (a)->pos.q               < (b)->pos.q + (b)->span.q )
#define mk_aln_ivt_cmp_iter(a, b)	( (a)->pos.q + (a)->span.q > (b)->pos.q + (b)->span.q )

static _force_inline
uint64_t mk_aln_ivt_update(mk_aln_t *parent, mk_aln_t *child)
{
	if(child == NULL) { return(1); }

	uint32_t cepos = child->pos.q + child->span.q;
	uint32_t qmax = parent->qmax;

	if(cepos > qmax) {
		parent->qmax = cepos;
		return(1);
	}
	return(0);
}

RBT_INIT_IVT(aln, mk_aln_t, mk_aln_rbt_header,
	mk_aln_rbt_cmp,
	mk_aln_ivt_update
);
RBT_INIT_ITER(aln, mk_aln_t, mk_aln_rbt_header,
	mk_aln_ivt_cmp_head,
	mk_aln_ivt_cmp_tail,
	mk_aln_ivt_cmp_iter
);


/* working buffer */
typedef struct {
	uint32_t dummy;
} mk_scan_conf_t;

typedef struct {
	struct {
		mk_ref_sketch_t const *ptr;
		size_t cnt;
	} ref;

	/* working buffers */
	struct {
		mk_seed_v arr;		/* for sorting seeds */
		mk_intv_v intv;
	} seed;

	/* chaining */
	struct {
		union {
			struct { uint32_t u, v; } sep;
			uint64_t all;
		} window;
		uint32_t min_scnt, chained;

		/* working buffer */
		mk_chain_v arr;
	} chain;

	/* chain filter */
	struct {
		uint8_t smask[16];
		uint8_t score_matrix[16];
		struct {
			uint8_t cv[16], pv[16];
		} init;

		/* extension test if < test_cnt, do inward extension if longer than uspan_thresh */
		uint32_t test_cnt, uspan_thresh;
		int32_t min_score;		/* discard if sum of extension scores is smaller than min_score */
	} filter;

	/* alignment bin */
	struct {
		uint32_t min_weight;
		dz_t *dz;				/* everything contained */

		/* result bin */
		mk_aln_v arr;
	} extend;
} mk_scan_t;


static _force_inline
void mk_scan_init_static(mk_scan_t *self, size_t cnt, mk_scan_conf_t const *conf)
{
	_unused(self);
	_unused(cnt);
	_unused(conf);

	return;
}

static _force_inline
void mk_scan_destroy_static(mk_scan_t *self, size_t cnt)
{
	_unused(self);
	_unused(cnt);

	return;
}

static _force_inline
void mk_scan_clear(mk_scan_t *self)
{
	kv_clear(self->seed.arr);
	kv_clear(self->seed.intv);
	kv_clear(self->chain.arr);
	kv_clear(self->extend.arr);
	return;
}


static _force_inline
mk_seed_t *mk_seed_reserve(mk_seed_v *buf, mk_seed_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	mk_seed_t *p = kv_reserve(mk_seed_t, *buf, 65536);
	return(&p[kv_cnt(*buf)]);
}

static _force_inline
mk_intv_t *mk_intv_reserve(mk_intv_v *buf, mk_intv_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	mk_intv_t *p = kv_reserve(mk_intv_t, *buf, 64);
	return(&p[kv_cnt(*buf)]);
}

static _force_inline
size_t mk_save_intv(mk_ref_state_t s, mk_ref_squash_t sq, mk_intv_t *r)
{
	/* just copy */
	*r = (mk_intv_t){
		.src = sq.src,
		.dst = sq.dst,
		.cnt = sq.cnt
	};
	return(1 - s.unmatching);
}

static _force_inline
size_t mk_expand_seed(mk_ref_state_t s, v4i32_t uofs, v4i32_t vofs, mk_seed_t *q)
{
	v4i32_t const wmask = _set_v4i32(0x20000000);

	mk_ref_match_t m = mk_ref_get_arr(s);
	for(size_t i = 0; i < m.cnt; i += 4) {
		v4i32_t const z = ({
			__m128i const x = _mm_loadl_epi64((__m128i const *)&m.ptr[i]);	/* movq (mem), xmm */
			__m128i const y = _mm_cvtepi16_epi32(x);	/* pmovsxwd; sign expansion */
			((v4i32_t){ .v1 = y });
		});
		v4i32_t const w = _and_v4i32(z, wmask);			/* disjoin forward and reverse */

		v4i32_t const u = _sub_v4i32(uofs, w);
		v4i32_t const v = _add_v4i32(vofs, _add_v4i32(w, w));

		/* save */
		_storeu_v4i32(&q[i],     _lo_v4i32(u, v));
		_storeu_v4i32(&q[i + 2], _hi_v4i32(u, v));
	}
	return(m.cnt);
}

static _force_inline
size_t mk_collect_seed(mk_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, mk_seed_v *seed, mk_intv_v *intv)
{
	/* load coordinate constants */
	v4i32_t const uinc = _set_v4i32(2), vinc = _set_v4i32(-1);
	v4i32_t uofs = _set_v4i32(0x20000), vofs = _set_v4i32(0x10000);

	mk_seed_t *q = mk_seed_reserve(seed, kv_ptr(*seed));
	mk_intv_t *r = mk_intv_reserve(intv, kv_ptr(*intv));

	/* initial state (of three general-purpose registers) */
	mk_ref_state_t s = mk_ref_match_init(ref);

	/* for each base... */
	for(size_t i = 0; i < qlen; i++) {
		/* update coordinates */
		uofs = _add_v4i32(uofs, uinc);
		vofs = _add_v4i32(vofs, vinc);

		/* update matching status for the next bin (prefetch) */
		mk_ref_next_t n = mk_ref_match_next(s, query[i]);
		r += mk_save_intv(n.state, n.squash, r);

		/* skip current bin if not matching */
		if(!s.unmatching) {
			q += mk_expand_seed(s, uofs, vofs, q);
		}

		/* alias new state (we expect the next bin has already arrived) */
		s = n.state;

		/* test array expansion every 32 times */
		if((i & 0x1f) != 0) { continue; }
		q = mk_seed_reserve(seed, q);
		r = mk_intv_reserve(intv, r);
	}
	return(kv_cnt(*seed));
}

static _force_inline
size_t mk_squash_seed(mk_intv_t const *intv, size_t icnt, mk_seed_v *seed)
{
	mk_seed_t *q = kv_ptr(*seed);

	/* squash succeeding match */
	size_t idx = 0;
	for(size_t i = 0; i < icnt; i++) {
		mk_intv_t const *p = &intv[i];

		/* determine source and destination indices */
		size_t const src = idx + p->src;
		size_t const dst = idx + p->dst;
		size_t const tail = idx + p->cnt;

		/* copy with ymm */
		for(size_t i = src, j = dst; j < tail; i += 8, j += 8) {
			v32i8_t const v = _loadu_v32i8(&q[i]);
			_storeu_v32i8(&q[j], v);
		}

		/* update index for next iteration */
		idx += tail - dst + src;
	}

	kv_cnt(*seed) = idx;
	return(idx);
}

static _force_inline
int64_t mk_chain_test_ptr(mk_seed_t const *p, mk_seed_t const *t)
{
	return((int64_t)((ptrdiff_t)(t - p - 1)));
}

static _force_inline
mk_seed_t *mk_chain_find_first(mk_scan_t const *self, mk_seed_t *p, mk_seed_t *t, uint64_t lb)
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
mk_seed_t *mk_chain_find_alt(mk_scan_t const *self, mk_seed_t *p, mk_seed_t *t, uint64_t lb)
{
	_unused(self);

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
mk_seed_t *mk_chain_find_nearest(mk_scan_t const *self, mk_seed_t *p, mk_seed_t *t)
{
	/* load root positions */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const uv = _loadu_u64(p);			/* (ulb, vlb) */
	uint64_t lb = (uv & ~umask);				/* (0, vlb) */

	/* find first chainable seed */
	mk_seed_t *n = mk_chain_find_first(self, p, t, lb);
	if(n == NULL) { return(NULL); }				/* we expect control paths are kept isolated */

	/* then inspect alternative chainable seeds */
	return(mk_chain_find_alt(self, n, t, lb));
}

static _force_inline
size_t mk_chain_record(mk_scan_t const *self, v2i32_t root, v2i32_t tail, v2i32_t scnt, mk_chain_t *q)
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

	/* forward if scnt > min_scnt */
	v2i32_t const min_scnt = _load_v2i32(&self->chain.min_scnt);
	v2i32_t const inc = _set_v2i32(1);		/* kept on register */
	v2i32_t const fwd = _and_v2i32(inc, _gt_v2i32(scnt, min_scnt));
	return(_ext_v2i32(fwd, 0));				/* movq */
}

static _force_inline
size_t mk_chain_seed(mk_scan_t const *self, mk_seed_t *seed, size_t scnt, mk_chain_v *chain)
{
	/*
	 * greedy extension; relatively smaller window size (~32) is preferred for good performance.
	 * do all of the chainablity test on general-purpose register for smaller latency of core loop.
	 * keep seed and counter on xmm regsters and do everything of the chain composing on xmm registers
	 * to minimize spill of general-purpose registers.
	 */
	uint64_t const chained = self->chain.chained;	/* offset */
	v2i32_t const inc = _set_v2i32(1);

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
v16i8_t mk_load_forward(mk_scan_t const *self, uint8_t const *ptr)
{
	_unused(self);
	return(_loadu_v16i8(ptr));
}

static _force_inline
v16i8_t mk_load_reverse(mk_scan_t const *self, uint8_t const *ptr)
{
	v16i8_t const mask = _load_v16i8(self->filter.smask);
	v16i8_t const v = _loadu_v16i8(ptr);

	/* flip */
	return(_shuf_v16i8(mask, v));
}

typedef struct {
	uint8_t const *ptr;
	v16i8_t (*load)(mk_scan_t const *, uint8_t const *);
	ptrdiff_t fwd;
} mk_filter_t;

static _force_inline
int64_t mk_filter_extend(mk_scan_t const *self, mk_filter_t r, mk_filter_t q)
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
		rv = _bsrd_v16i8(rv, rn, 1); rn = _shr_v16i8(rn, 1);

		/* #1 (odd) */
		v16i8_t const a1 = _calc_a(rv, qv, tdh, tdv);
		dv = _sub_v16i8(a1, tdh);
		dh = _sub_v16i8(a1, tdv);
		sv = _add_v16i8(sv, dh);
		mv = _max_v16i8(mv, sv);

		/* move downward */
		qv = _bsrd_v16i8(qv, qn, 1); qn = _shr_v16i8(qn, 1);
	}

	/* fold scores */
	int64_t const score = 0;
	return(score);

	#undef _calc_a
}

static _force_inline
int64_t mk_filter_extend_inward(mk_scan_t const *self, mk_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, mk_chain_t const *chain)
{
	_unused(self);
	_unused(ref);
	_unused(query);
	_unused(qlen);
	_unused(chain);
	return(0);
}

static _force_inline
int64_t mk_filter_extend_outward(mk_scan_t const *self, mk_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, mk_chain_t const *chain)
{
	_unused(self);
	_unused(ref);
	_unused(query);
	_unused(qlen);
	_unused(chain);
	return(0);
}

static _force_inline
size_t mk_filter_save_chain(mk_scan_t const *self, uint32_t rid, mk_chain_t *q, mk_chain_t const *p)
{
	_unused(self);

	/* load */
	v16i8_t v = _loadu_v16i8(p);

	/* calc weight */
	ZCNT_RESULT size_t zc = _lzcnt_u32(p->attr.scnt);
	uint32_t const weight = 32 - zc;

	/* save */
	_storeu_v16i8(q, v);
	q->attr.stat.rid = rid;
	q->attr.stat.weight = weight;
	return(1);
}

static _force_inline
size_t mk_filter_chain(mk_scan_t const *self, mk_ref_sketch_t const *ref, uint32_t rid, uint8_t const *query, size_t qlen, mk_chain_t *chain, size_t ccnt)
{
	/* alias pointer */
	mk_chain_t const *src = chain;
	mk_chain_t *dst = chain;

	/* load constants */
	uint64_t const test_cnt = self->filter.test_cnt;
	int64_t const min_score = self->filter.min_score;
	size_t const uspan_thresh = self->filter.uspan_thresh;

	for(size_t i = 0; i < ccnt; i++) {
		mk_chain_t const *p = &src[i];

		/* try short extension if not heavy enough */
		if(p->attr.scnt < test_cnt) {
			int64_t score = (p->span.u > uspan_thresh
				? mk_filter_extend_inward(self, ref, query, qlen, p)		/* inward extension for long-spanning chain */
				: mk_filter_extend_outward(self, ref, query, qlen, p)
			);
			if(score < min_score) { continue; }
		}

		/* copy and fold in reference id */
		dst += mk_filter_save_chain(self, rid, dst, p);
	}
	return(ccnt);
}

static _force_inline
size_t mk_collect_candidate(mk_scan_t *self, mk_ref_sketch_t const *ref, uint32_t rid, uint8_t const *query, size_t qlen, mk_chain_v *chain)
{
	/* enumerate seeds */
	kv_clear(self->seed.arr);
	if(mk_collect_seed(ref, query, qlen, &self->seed.arr, &self->seed.intv) == 0) {
		return(0);
	}

	mk_intv_t const *intv = kv_ptr(self->seed.intv);
	size_t icnt = kv_cnt(self->seed.intv);

	/* squash overlapping seeds */
	if(mk_squash_seed(intv, icnt, &self->seed.arr) == 0) {
		return(0);
	}

	mk_seed_t *sptr = kv_ptr(self->seed.arr);
	size_t scnt = kv_cnt(self->seed.arr);

	/* sort seeds by u-coordinates for chaining */
	radix_sort_seed(sptr, scnt);

	/* chaining; return if no chain obtained */
	size_t const cbase = kv_cnt(*chain);		/* save chain count */
	size_t const ccnt  = mk_chain_seed(self, sptr, scnt, chain);
	if(ccnt == 0) { return(0); }

	/* filter (try small extension with simple score matrix) */
	size_t const cfilt = mk_filter_chain(self, ref, rid, query, qlen, &kv_ptr(*chain)[cbase], ccnt);
	if(cfilt == 0) { return(0); }

	/* finalize filtered chains */
	kv_cnt(*chain) = cbase + cfilt; 
	return(kv_cnt(*chain));
}



/* X-drop DP extension */
static _force_inline
uint64_t mk_extend_is_complete(mk_scan_t *self)
{
	_unused(self);

	/* FIXME */
	return(0);
}

static _force_inline
uint64_t mk_extend_is_covered(mk_scan_t *self, mk_aln_t const *caln, uint32_t weight)
{
	/* convert chain position to (pseudo) alignment range */
	mk_aln_t const *v = kv_ptr(self->extend.arr);

	/* for each overlapping alignment */
	rbt_iter_t it;
	rbt_init_iter_aln(&it, v, caln);

	mk_aln_t const *p = rbt_fetch_head_aln(&it, v, caln);
	while(p != NULL) {
		if(p->min_weight > weight) {
			return(1);		/* already covered */
		}

		p = rbt_fetch_next_aln(&it, v, caln);
	}
	return(0);				/* possiblilty remains for better alignment */
}

static _force_inline
void mk_extend_record(mk_scan_t *self, mk_chain_t const *chain, mk_aln_t aln)
{
	_unused(chain);

	kv_push(mk_aln_t, self->extend.arr, aln);
	rbt_insert_aln(kv_ptr(self->extend.arr), kv_cnt(self->extend.arr) - 1);
	return;
}

static _force_inline
void mk_extend_clear(mk_scan_t *self)
{
	kv_clear(self->extend.arr);
	return;
}


static _force_inline
mk_aln_t mk_chain_as_aln(mk_scan_t const *self, mk_chain_t const *chain)
{
	_unused(self);
	_unused(chain);

	/* FIXME */
	return((mk_aln_t){
		.pos  = { .r = 0, .q = 0 },
		.span = { .r = 0, .q = 0 },
		.qmax = 0
	});
}

static _force_inline
mk_aln_t mk_extend_seed(mk_scan_t const *self, mk_ref_sketch_t const *ref, dz_query_t *const *query, mk_aln_t const *caln)
{
	/* upward then downward */
	int64_t rrpos = caln->pos.r;
	dz_forefront_t const *r = dz_extend(self->extend.dz, query[1],
		dz_root(self->extend.dz), 1,
		(char const *)&ref->seq[rrpos], -rrpos, 0
	);

	/* turnback */
	uint64_t spos = dz_calc_max_pos(self->extend.dz, r);
	uint64_t rspos = spos>>32, qspos = spos & 0xffffffff;
	dz_forefront_t const *f = dz_extend(self->extend.dz, query[0],
		dz_root(self->extend.dz), 1,
		(char const *)&ref->seq[rspos], ref->slen - rspos, 0
	);

	/* traceback */
	dz_alignment_t *aln = dz_trace(self->extend.dz, f);

	uint32_t rspan = aln->ref_length, qspan = aln->query_length;
	return((mk_aln_t){
		/* coordinates */
		.qmax = qspos + qspan,
		.pos  = { .r = rspos, .q = qspos },
		.span = { .r = rspan, .q = qspan },

		/* path */
		.path = { .ptr = aln->path, .len = aln->path_length },

		/* stats */
		.min_weight = 0,
		.score = aln->score,
		.aln = aln				/* save original */
	});
}

static _force_inline
size_t mk_extend(mk_scan_t *self, mk_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, mk_chain_t const *chain, size_t ccnt)
{
	/* clear bin */
	mk_extend_clear(self);
	dz_flush(self->extend.dz);

	/* build reference-side profile */
	dz_query_t *q[2] = {
		dz_pack_query_forward(self->extend.dz, (char const *)query, qlen),
		dz_pack_query_reverse(self->extend.dz, (char const *)query, qlen)
	};

	/* for each chain try extension */
	for(size_t i = 0; i < ccnt; i++) {
		mk_chain_t const *p = &chain[i];
		if(p->attr.stat.weight < self->extend.min_weight) { break; }

		mk_aln_t caln = mk_chain_as_aln(self, p);

		/* get reference sequence */
		mk_ref_sketch_t const *r = &ref[p->attr.stat.rid];

		/* skip if already covered (and there is no possibility that this chain surpasses the previous ones) */
		if(mk_extend_is_covered(self, &caln, p->attr.stat.weight)) { continue; }

		mk_aln_t aln = mk_extend_seed(self, r, &q[0], &caln);
		if(aln.score <= 0) { continue; }

		mk_extend_record(self, p, aln);
		if(mk_extend_is_complete(self)) { break; }
	}

	/* anything else? */
	return(kv_cnt(self->extend.arr));
}


/* evaluate all; query sequence be shorter than 2Gbp */
static _force_inline
size_t mk_scan_mask(mk_scan_t *self, uint8_t const *seq, size_t slen)
{
	mk_scan_clear(self);

	/* collect chain for each reference sequence */
	mk_ref_sketch_t const *ref = self->ref.ptr;
	for(size_t i = 0; i < self->ref.cnt; i++) {
		mk_collect_candidate(self, &ref[i], i, seq, slen, &self->chain.arr);
	}

	mk_chain_t *cptr = kv_ptr(self->chain.arr);
	size_t const ccnt = kv_cnt(self->chain.arr);
	if(ccnt == 0) { return(0); }

	/* sort by estimated identity */
	radix_sort_chain(cptr, ccnt);

	/* convert chain to alignment */
	size_t acnt = mk_extend(self, ref, seq, slen, cptr, ccnt);
	return(acnt);
}




/* printer */
typedef struct {
	uint32_t dummy;
} mk_print_conf_t;

typedef struct {
	uint32_t dummy;
} mk_print_t;


static _force_inline
void mk_print_destory_static(mk_print_t *self)
{
	_unused(self);

	return;
}

static _force_inline
void mk_print_init_static(mk_print_t *self, mk_print_conf_t const *conf, char const *args)
{
	_unused(self);
	_unused(conf);
	_unused(args);

	return;
}





/* multithreaded scan-and-mask */
typedef struct {
	uint32_t dummy;
} mk_mtscan_t;

static _force_inline
mk_mtscan_t *mk_mtscan_init(mk_scan_conf_t const *conf, mk_idx_t const *mi, mk_print_t *printer, pt_t *pt)
{
	_unused(conf);
	_unused(mi);
	_unused(printer);
	_unused(pt);

	return(NULL);
}

static _force_inline
void mk_mtscan_destroy(mk_mtscan_t *self)
{
	free(self);
	return;
}

static _force_inline
int mk_mtscan_file(mk_mtscan_t *self, char const *fn)
{
	_unused(self);
	_unused(fn);

	return(0);
}


/* return codes */
enum main_error_codes {
	ERROR_INTERNAL = 255,
	ERROR_NO_ARG = 1,

	ERROR_OPEN_IDX = 2,
	ERROR_LOAD_IDX = 3,

	ERROR_OPEN_RSEQ = 4,
	ERROR_OPEN_QSEQ = 5
};

typedef struct {
	/* index dump mode if not NULL */
	char const *idxdump;

	/* global args */
	char const *args;
	uint32_t verbose, help;
	size_t nth;

	/* scan-and-mask params */
	mk_idx_conf_t index;
	mk_scan_conf_t scan;
	mk_print_conf_t print;

	/* option parser */
	FILE *log;
	opt_t opt;
} mk_conf_t;

static _force_inline
void mk_conf_init_static(mk_conf_t *conf, char const *const *argv, FILE *fp)
{
	_unused(conf);
	_unused(argv);
	_unused(fp);
	return;
}

static _force_inline
void mk_conf_destroy_static(mk_conf_t *conf)
{
	_unused(conf);
	return;
}


/* determine help and verbose level */
typedef struct {
	FILE *fp;
	uint64_t help, quit;
} mk_conf_outfp_t;

static _force_inline
mk_conf_outfp_t mk_conf_get_outfp(mk_conf_t *conf)
{
	/* use stdout for explicit version (-v) and help (-h) options */
	if(conf->verbose == 1 || conf->help > 0) {
		return((mk_conf_outfp_t){
			.fp = stdout,
			.help = conf->help,
			.quit = 1
		});
	}

	/* restore default verbose level */
	if(conf->verbose == 0) {
		conf->verbose = 1;
	}

	/* implicit help invoked when no input file found; always redirected to stderr */
	return((mk_conf_outfp_t){
		.fp = conf->log,
		.help = opt_parg_cnt(&conf->opt) == 0,
		.quit = 0
	});
}

/* print help */
static _force_inline
void mk_conf_print_help(mk_conf_t const *conf, FILE *log)
{
	if(conf->verbose == 0) { return; }

	#define _msg_impl(_level, _fmt, ...) { logger_printf(log, _fmt "%s\n", __VA_ARGS__); }
	#define _msg(_level, ...) { if((_level) <= conf->verbose + 1) { _msg_impl(_level, __VA_ARGS__, ""); } }

	_msg(2, "\n"
			"  tinymasker - fast repeat masking tool\n"
			"");
	_msg(2, "Usage:\n"
			"  index construction:\n"
			"    $ tinymasker -t4 -d index.tmi repeats.fa\n"
			"  mask:\n"
			"    $ tinymasker -t4 index.tmi contigs.fa\n"
			"")
	_msg(2, "General options:");
	_msg(2, "");
	if(conf->verbose < 2) {
		_msg(2, "  Pass -hh to show all the options.");
		_msg(2, "");
	}

	#undef _msg_impl
	#undef _msg

	return;
}


/* index construction */

static _force_inline
int main_index_error(mk_conf_t *conf, int error_code, char const *filename)
{
	_unused(conf);
	switch(error_code) {
	/* argument missing */
	case ERROR_NO_ARG: error("argument is not enough. at least one reference file is required."); break;

	/* opening files */
	case ERROR_OPEN_IDX: error("failed to open index file `%s' in write mode. Please check file path and its permission.", filename); break;
	case ERROR_OPEN_RSEQ: error("failed to open sequence file `%s' for building index. Please check file path and its format.", filename); break;
	}
	return(error_code);
}

static _force_inline
int main_index_intl(mk_conf_t *conf, pg_t *pg, pt_t *pt)
{
	if(opt_parg_cnt(&conf->opt) == 0) {
		return(main_index_error(conf, ERROR_NO_ARG, NULL));
	}

	/* iterate over index blocks */
	kv_foreach(void *, opt_pargv(&conf->opt), {

		/* generate index */
		mk_idx_t *mi = mk_idx_gen(&conf->index, pt, *p, stderr);
		if(mi == NULL) {
			/* failed to open file */
			return(main_index_error(conf, ERROR_OPEN_RSEQ, conf->idxdump));
		}

		/* dump index */
		size_t size = mk_idx_dump(mi, pg, (write_t)pgwrite);

		/* flush output for next batch */
		pg_flush(pg);

		mk_idx_destroy(mi);
		message(conf->log, "built and dumped index for `%s', on-memory size of this chunk: %.1f MB", (char const *)*p, (double)size / (1024ULL * 1024));
	});
	return(0);
}

static _force_inline
int main_index(mk_conf_t *conf, pt_t *pt)
{
	/* add suffix if missing and if /dev/xxx */
	if(!mm_startswith(conf->idxdump, "/dev") && !mm_endswith(conf->idxdump, ".mai")) {
		message(conf->log, "index filename does not end with `.mai' (added).");
		conf->idxdump = opt_append(&conf->opt, conf->idxdump, ".mai");
	}

	/* open file in write mode */
	FILE *fp = fopen(conf->idxdump, "wb");
	if(fp == NULL) {
		/* failed open file (locked?) */
		return(main_index_error(conf, ERROR_OPEN_IDX, conf->idxdump));
	}

	/* initialize compressor */
	pg_t pg;
	pg_init_static(&pg, fp, pt);

	/* index_intl does everything */
	int error_code = main_index_intl(conf, &pg, pt);

	/* done */
	pg_stat_t stat = pg_destroy_static(&pg);
	message(conf->log, "done. total index size (compressed) on disk: %.1f MB.", (double)stat.out / (1024ULL * 1024));

	fclose(fp);
	return(error_code);
}


/* scan-and-mask */

static _force_inline
int main_scan_error(mk_conf_t *conf, int error_code, char const *file)
{
	_unused(conf);
	switch(error_code) {
	/* unknown */
	case ERROR_INTERNAL: error("failed to instanciate alignment context."); break;

	/* in mapping */
	case ERROR_OPEN_QSEQ: error("failed to open sequence file `%s'. Please check file path and format.", file); break;

	/* argument missing */
	case ERROR_NO_ARG: error("argument is not enough. a reference and at least one query files are required."); break;

	/* index loading */
	case ERROR_OPEN_IDX: error("failed to open index file `%s'. Please check file path and permission.", file); break;
	case ERROR_LOAD_IDX: error("failed to load index block from `%s'. Please check file path and version, or rebuild the index.", file); break;

	/* index construction */
	case ERROR_OPEN_RSEQ: error("failed to open sequence file `%s' for building index. Please check file path and format.", file); break;
	}
	return(error_code);
}

/* working buffer */
typedef struct {
	/* always available */
	mk_conf_t *conf;
	size_t fcnt;					/* fetched count */

	pt_t *pt;
	mk_print_t *printer;

	char const *ref;				/* reference filename */
	char const *const *query;		/* query filename */
	size_t qcnt;


	/* for prebuilt index */
	FILE *fp;				/* != NULL if prebuilt index is available */
	pg_t pg;
} main_scan_tbuf_t;

static _force_inline
int main_scan_tbuf_destroy_static(main_scan_tbuf_t *w)
{
	if(w->fp != NULL) {
		fclose(w->fp);
		pg_destroy_static(&w->pg);
	}
	return(0);
}

static _force_inline
int main_scan_tbuf_init_static(main_scan_tbuf_t *w, mk_conf_t *conf, char const *const *parg, size_t pcnt, mk_print_t *printer, pt_t *pt)
{
	/* save args */
	*w = (main_scan_tbuf_t){
		.conf = conf,
		.fcnt = 0,

		.pt = pt,
		.printer = printer,

		/* inputs */
		.ref = parg[0],
		.query = &parg[1],
		.qcnt = pcnt - 1
	};

	/* check suffix to determine if the file is prebuilt index */
	if(!mm_endswith(parg[0], ".mai")) {
		return(0);			/* not a prebuilt index */
	}

	/* open; if fails, it would be invalid path or permission */
	if((w->fp = fopen(parg[0], "rb")) == NULL) {
		return(main_scan_error(conf, ERROR_OPEN_IDX, parg[0]));
	}
	pg_init_static(&w->pg, w->fp, pt);
	return(0);
}


static _force_inline
int main_scan_idx_gen(main_scan_tbuf_t *w, mk_idx_t **pmi)
{
	if(++w->fcnt > 1) { return(0); }

	/* first call */
	mk_idx_t *mi = mk_idx_gen(&w->conf->index, w->pt, w->ref, stderr);
	if(mi == NULL) {
		/* failed to open file */
		return(main_scan_error(w->conf, ERROR_OPEN_RSEQ, w->ref));
	}
	message(stderr, "built index for `%s'.", w->ref);
	*pmi = mi;
	return(0);
}

static _force_inline
int main_scan_idx_load(main_scan_tbuf_t *w, mk_idx_t **pmi)
{
	if(pg_eof(&w->pg)) { return(0); }

	/* prebuilt index available; try to fetch next block */
	mk_idx_t *mi = mk_idx_load(&w->pg, (read_t)pgread);
	if(++w->fcnt == 1 && mi == NULL) {
		/* would be broken */
		return(main_scan_error(w->conf, ERROR_LOAD_IDX, w->ref));
	}

	/* NULL for tail */
	if(mi != NULL) { message(stderr, "loaded index block from `%s'.", w->ref); }
	*pmi = mi;
	return(0);
}

/* for each query file */
static _force_inline
int main_scan_foreach_qfile(main_scan_tbuf_t *w, mk_mtscan_t *mt)
{
	for(size_t i = 0; i < w->qcnt; i++) {
		if(mk_mtscan_file(mt, w->query[i])) {
			return(main_scan_error(w->conf, ERROR_OPEN_QSEQ, NULL));
		}
	}
	return(0);
}

/* for each index chunk */
static _force_inline
int main_scan_foreach_idx(main_scan_tbuf_t *w)
{
	while(1) {
		mk_idx_t *mi = NULL;

		int fetcher_error_code = (w->fp == NULL
			? main_scan_idx_gen(w, &mi)
			: main_scan_idx_load(w, &mi)
		);
		if(fetcher_error_code != 0 || mi == NULL) {
			return(fetcher_error_code);		/* error should be handled inside */
		}

		/* instanciate multithreading context for this index chunk */
		mk_mtscan_t *mt = mk_mtscan_init(&w->conf->scan, mi, w->printer, w->pt);
		if(mt == NULL) {
			return(main_scan_error(w->conf, ERROR_INTERNAL, NULL));
		}

		/* do the task; for each query file */
		int scan_error_code = main_scan_foreach_qfile(w, mt);

		/* done for this index chunk */
		mk_mtscan_destroy(mt);
		mk_idx_destroy(mi);
		if(scan_error_code != 0) { return(scan_error_code); }
	}
	return(0);
}

/* for each index filename */
static _force_inline
int main_scan_intl(mk_conf_t *conf, mk_print_t *printer, pt_t *pt)
{
	char const *const *parg = (char const *const *)opt_parg(&conf->opt);
	size_t pcnt = opt_parg_cnt(&conf->opt);

	/* error if no argument is given */
	if(pcnt < 2) {
		return(main_scan_error(conf, ERROR_NO_ARG, NULL));
	}

	/* instanciate thread-local working buffers */
	main_scan_tbuf_t w;
	int init_error_code = main_scan_tbuf_init_static(&w, conf, parg, pcnt, printer, pt);
	if(init_error_code != 0) {
		return(init_error_code);
	}

	/* we consider the first argument as reference */
	int scan_error_code = main_scan_foreach_idx(&w);

	/* done */
	main_scan_tbuf_destroy_static(&w);
	return(scan_error_code);
}

/* printer */
static _force_inline
int main_scan(mk_conf_t *conf, pt_t *pt)
{
	/* instanciate alignment formatter */
	mk_print_t printer;
	mk_print_init_static(&printer, &conf->print, conf->args);

	/* dispatch */
	int error_code = main_scan_intl(conf, &printer, pt);

	mk_print_destory_static(&printer);
	return(error_code);
}

/* create worker threads for indexing and mapping */
static _force_inline
int main_dispatch(mk_conf_t *conf)
{
	pt_t *pt = pt_init(conf->nth);
	if(pt == NULL) {
		return(main_scan_error(conf, ERROR_INTERNAL, NULL));
	}

	/* dispatch either index or scan */
	int error_code = (conf->idxdump ? main_index : main_scan)(conf, pt);

	/* done */
	pt_destroy(pt);
	return(error_code);
}


/* entry */
int main(int argc, char *argv[])
{
	int error_code = ERROR_NO_ARG;

	/* unittest dispatcher */
	#if defined(UNITTEST) && UNITTEST != 0
		if(argc > 1 && strcmp(argv[1], "unittest") == 0) {
			return(unittest_main(argc, argv));
		}
	#else
		_unused(argc);
	#endif

	/* instanciate option object */
	logger_init();
	mk_conf_t conf;
	mk_conf_init_static(&conf, (char const *const *)argv, stderr);

	/* always print version */
	mk_conf_outfp_t out = mk_conf_get_outfp(&conf);
	message(out.fp, "Version: %s, Build: %s", mk_version(), MK_ARCH_NAME);

	/* when -h is passed or no input file is given, print help message */
	if(out.help) {
		if(conf.help > 0) { error_code = 0; }	/* also exit status is 0 (not an error) */
		mk_conf_print_help(&conf, out.fp);		/* we use stdout when invoked by -h option */
	}
	if(out.quit) { goto _main_final; }

	/* dispatch tasks and get return code */
	if((error_code = main_dispatch(&conf)) == 0) {
		message(conf.log, "Command: %s", conf.args);	/* print log when succeeded */
	}
	logger_destroy();

_main_final:;
	mk_conf_destroy_static(&conf);
	return(error_code);
}

/**
 * end of tinymasker.c
 */
