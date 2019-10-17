
/**
 * @file align.c
 * @brief seeding, chaining, filtering and alignment generation (core implementations)
 *
 * @author Hajime Suzuki
 * @license MIT
 */


#ifndef UNITTEST
#  define UNITTEST				( 1 )
#endif
#define UNITTEST_UNIQUE_ID		3


#include "utils/utils.h"		/* include all */
#include "tinymasker.h"
unittest_config( .name = "align" );


/* others */
#include "dozeu.h"
#include "baln.h"
#include "dbg.h"
#include "index.h"


/* forward decl */
#include "align.h"


/* general position */
typedef struct {
	uint32_t r, q;
} tm_pair_t;

static _force_inline
tm_pair_t tm_add_pair(tm_pair_t a, tm_pair_t b)
{
	return((tm_pair_t){
		.r = a.r + b.r,
		.q = a.q + b.q
	});
}


/* position with direction */
typedef struct {
	uint16_t r;
	uint8_t dir, unused;		/* dir == 1 for reverse, 0 for forward */
	uint32_t q;
} tm_pos_t;

static _force_inline
char *tm_pos_to_str(tm_pos_t const *p)
{
	return(xbprintf("dir(%u), pos(%u, %u), unused(%u)", p->dir, p->r, p->q, p->unused));
}

static _force_inline
tm_pair_t tm_pos_to_pair(tm_pos_t pos)
{
	/* drop dir and unused */
	return((tm_pair_t){
		.r = pos.r,
		.q = pos.q
	});
}

static _force_inline
tm_pos_t tm_add_pos(tm_pos_t a, tm_pair_t b)
{
	return((tm_pos_t){
		.r   = a.r + (a.dir ? -b.r : b.r),
		.q   = a.q + b.q,
		.dir = a.dir,
		.unused = a.unused
	});
}

static _force_inline
tm_pos_t tm_sub_pos(tm_pos_t a, tm_pair_t b)
{
	return((tm_pos_t){
		.r   = a.r + (a.dir ? b.r : -b.r),
		.q   = a.q - b.q,
		.dir = a.dir,
		.unused = a.unused
	});
}

static _force_inline
tm_pos_t tm_canon_pos(tm_pos_t pos, tm_pair_t span)
{
	return((tm_pos_t){
		.r   = pos.r - (pos.dir ? span.r : 0),	/* convert to head position */
		.q   = pos.q,
		.dir = pos.dir,
		.unused = pos.unused
	});
}



/* alignment scores */
typedef struct {
	int32_t raw;			/* raw SW score + end bonus */
	int32_t patched;		/* raw SW score + end bonus + adjustment */
} tm_score_t;


/* alignment result */
typedef struct {
	/* rbt */
	rbt_header_t h;
	uint32_t qmax;

	/* positions */
	struct {
		uint32_t rid : 24;
		uint32_t max_weight : 8;
	} attr;

	tm_pos_t pos;
	tm_pair_t span;

	/* scores */
	tm_score_t score;

	/* everything else in dz_alignment_t */
	dz_alignment_t const *aln;
} tm_aln_t;
_static_assert(sizeof(tm_aln_t) == 48);



/* seed; see comment below for the definition */
typedef struct {
	uint32_t v, u;
} tm_seed_t;
_static_assert(sizeof(tm_seed_t) == 8);
typedef struct { tm_seed_t *a; size_t n, m; } tm_seed_v;

#define tm_seed_upos(x)			( (x).u )
KRADIX_SORT_INIT(seed, tm_seed_t, tm_seed_upos, 4);

#if 0
static _force_inline
tm_pair_t tm_seed_decode(tm_seed_t const *p)
{
	/* remove chained flag */
	uint32_t const u = p->u, v = p->v;

	/*
	 * convert to r-q corrdinates;
	 * u = 2 * q - r + 0x10000 + d
	 * v = 2 * r - q - 0x20000 - 2 * d
	 *     where d = 0 for forward and 0x40000000 for reverse
	 *
	 * 2 * u + v = 3 * q
	 * 2 * v + u = 3 * (r - 0x10000 - d)
	 */
	uint32_t const rt = (2 * v + u) & 0x7fffffff, r = (rt + 0x80030000) / 3U;
	uint32_t const q  = (2 * u + v) / 3U;

	uint32_t const vv = v & 0x3fffffff;
	uint32_t const uu = u & 0x3fffffff;
	uint32_t const rt2 = 2 * vv + uu, rr = rt2 + 0x30000, r2 = (rr - 0xc0000000) / 3U;
	uint32_t const qt2 = 2 * uu + vv, qq = qt2 & 0x3fffffff, q2 = qq / 3U;
	debug("r(%x), q(%x), v(%x), u(%x), r2(%x, %x, %x), q2(%x, %x, %x)", r, q, vv, uu, rt2, rr, r2, qt2, qq, q2);

	/* subtract offset */
	return((tm_pair_t){
		.r = r2,		/* leaves direction flag at bit 30 */
		.q = q2
	});
}
#endif


static _force_inline
uint64_t tm_seed_dir(tm_seed_t const *s)
{
	return(1 - (s->v>>31));
}

static char *tm_seed_to_str(tm_seed_t const *s)
{
	// tm_pair_t p = tm_seed_decode(s);
	uint32_t const u = s->u & 0x7fffffff, v = s->v;		/* clear chained flag */
	uint32_t const rt = (2 * v + u), r = (rt + 0xc0030000) / 3U;
	uint32_t const qt = (2 * u + v), q = (qt & 0x7fffffff) / 3U;
	uint32_t const dir = (u & 0x80000000) != 0;
	return(xbprintf("seed(%x, %x), t(%x, %x), pos(%x, %x, %x)", s->v, s->u, qt, rt, dir, q, r));
}


typedef struct {
	uint16_t dst, src, cnt, unused;	/* upper 4bits must be cleared before use (overlaps with plen) */
} tm_sqiv_t;
_static_assert(sizeof(tm_sqiv_t) == 8);
typedef struct { tm_sqiv_t *a; size_t n, m; } tm_sqiv_v;



/* chain */

typedef struct {
	/*
	 * base seed position;
	 * rpos comes first because vpos comes first in seed_t.
	 * vpos is placed lower of (upos, vpos) tuple to make chaining simpler.
	 */
	tm_pos_t pos;				/* dir == 1 for reverse, 0 for forward */

	/* reference and query side spans */
	tm_pair_t span;
} tm_chain_raw_t;

static _force_inline
tm_pair_t tm_chain_raw_pos(tm_chain_raw_t const *chain)
{
	return((tm_pair_t){
		.q = chain->pos.q,
		.r = chain->pos.r		/* upper 16bits are masked out */
	});
}

static _force_inline
tm_pair_t tm_chain_raw_span(tm_chain_raw_t const *chain)
{
	return((tm_pair_t){
		.q = chain->span.q,
		.r = chain->span.r
	});
}


typedef struct {
	/* see tm_chain_compose */
	tm_pos_t pos;					/* chain end position; rspos + rspan for forward and rspos for reverse */

	union {
		struct {
			uint32_t rid : 24;		/* reference id */
			uint32_t weight : 8;	/* weight; calculated from seed count and sequence length */
		} sep;
		uint32_t all;
	} attr;

	/* half remaining of span */
	struct {
		uint32_t q;
	} span;
} tm_chain_t;
typedef struct { tm_chain_t *a; size_t n, m; } tm_chain_v;

_static_assert(sizeof(tm_chain_t) == 16);
_static_assert(sizeof(tm_chain_t) == sizeof(tm_chain_raw_t));

#define tm_chain_attr(x)		( (x).attr.all )
KRADIX_SORT_INIT(chain, tm_chain_t, tm_chain_attr, 4);


static char *tm_chain_raw_to_str(tm_chain_raw_t const *c)
{
	uint32_t const rspos = c->pos.dir ? c->pos.r + c->span.r : c->pos.r;
	uint32_t const repos = c->pos.dir ? c->pos.r : c->pos.r + c->span.r;
	return(xbprintf("dir(%u), (%u, %u) -- (%u, %u) --> (%u, %u)", c->pos.dir, c->pos.q, rspos, c->span.q, c->span.r, c->pos.q + c->span.q, repos));
}
static char *tm_chain_to_str(tm_chain_t const *c)
{
	return(xbprintf("dir(%u), pos(%u, %u), qspan(%u), rid(%u), weight(%u)", c->pos.dir & 0x01, c->pos.q, c->pos.r, c->span.q, c->attr.sep.rid, c->attr.sep.weight));
}

static _force_inline
tm_pair_t tm_chain_head(tm_chain_t const *chain)
{
	debug("(%x, %x)", chain->pos.q, chain->pos.r);

	return((tm_pair_t){
		.q = chain->pos.q,
		.r = chain->pos.r		/* upper 16bits are masked out */
	});
}

static _force_inline
tm_pos_t tm_chain_pos(tm_chain_t const *chain)
{
	return(chain->pos);
}


/* for debugging */
static
char *tm_aln_to_str(tm_aln_t const *p)
{
	return(xbprintf("p(%p), score(%d, %d), dir(%u), pos(%u, %u), span(%u, %u)", p, p->score.raw, p->score.patched, p->pos.dir, p->pos.q, p->pos.r, p->span.q, p->span.r));
}


/* alignment utils */

typedef struct {
	double identity;
} tm_aln_stat_t;

static _force_inline
tm_aln_stat_t tm_aln_calc_stat(tm_aln_t const *aln, uint64_t flip)
{
	/* match / ref_length */
	dz_alignment_t const *a = aln->aln;
	uint32_t const span = flip ? aln->span.q : aln->span.r;


	// fprintf(stderr, "count(%u, %u, %u, %u), span(%u)\n", a->match_count, a->mismatch_count, a->ins_count, a->del_count, span);
	return((tm_aln_stat_t){
		.identity = (double)a->match_count / (double)span
	});
}


/* alignment bin and interval tree */
typedef struct {
	tm_aln_t *a;
	size_t n, m;
} tm_aln_v;

#define tm_aln_rbt_header(a)		( &(a)->h )
#define tm_aln_rbt_cmp(a, b)		( (a)->pos.q <  (b)->pos.q )

/* matcher */
#define tm_aln_rbt_cmp_head(a, b)	( (a)->pos.q >= (b)->pos.q )
#define tm_aln_rbt_cmp_tail(a, b)	( (a)->pos.q <= (b)->pos.q )
#define tm_aln_rbt_cmp_iter(a, b)	( 1 )

/* intersection */
#if 1
#define tm_aln_ivt_cmp_head(a, b)	({ (a)->qmax                > (b)->pos.q; })
#define tm_aln_ivt_cmp_tail(a, b)	({ (a)->pos.q               < (b)->pos.q + (b)->span.q; })
#define tm_aln_ivt_cmp_iter(a, b)	({ (a)->pos.q + (a)->span.q > (b)->pos.q; })
#else
#define tm_aln_ivt_cmp_head(a, b)	({ \
	debug("compare head: a(%u, %u, %u), b(%u, %u, %u), cmp(%u > %u)", (a)->pos.q, (a)->pos.q + (a)->span.q, (a)->qmax, (b)->pos.q, (b)->pos.q + (b)->span.q, (b)->qmax, (a)->qmax, (b)->pos.q); \
	(a)->qmax                > (b)->pos.q; \
})
#define tm_aln_ivt_cmp_tail(a, b)	({ \
	debug("compare tail: a(%u, %u, %u), b(%u, %u, %u), cmp(%u < %u)", (a)->pos.q, (a)->pos.q + (a)->span.q, (a)->qmax, (b)->pos.q, (b)->pos.q + (b)->span.q, (b)->qmax, (a)->pos.q, (b)->pos.q + (b)->span.q); \
	(a)->pos.q               < (b)->pos.q + (b)->span.q; \
})
#define tm_aln_ivt_cmp_iter(a, b)	({ \
	debug("compare iter: a(%u, %u, %u), b(%u, %u, %u), cmp(%u > %u)", (a)->pos.q, (a)->pos.q + (a)->span.q, (a)->qmax, (b)->pos.q, (b)->pos.q + (b)->span.q, (b)->qmax, (a)->pos.q + (a)->span.q, (b)->pos.q); \
	(a)->pos.q + (a)->span.q > (b)->pos.q; \
})
#endif

static _force_inline
uint64_t tm_aln_ivt_update(tm_aln_t *node, tm_aln_t const *left, tm_aln_t const *right)
{
	debug("update, node(%p), left(%p), right(%p)", node, left, right);

	if(left == NULL && right == NULL) {
		debug("initialize, pos(%u), span(%u), qmax(%u)", node->pos.q, node->span.q, node->pos.q + node->span.q);
		node->qmax = node->pos.q + node->span.q;	/* initialize */
		return(1);		/* further update needed */
	}

	/* calculate max of children */
	uint32_t const lmax = left  == NULL ? 0 : left->qmax;
	uint32_t const rmax = right == NULL ? 0 : right->qmax;
	uint32_t const cmax = MAX2(lmax, rmax);

	debug("lmax(%u), rmax(%u), cmax(%u), qmax(%u)", lmax, rmax, cmax, node->qmax);

	/* update if needed */
	if(cmax > node->qmax) {
		node->qmax = cmax;
		return(1);		/* further update needed */
	}
	return(0);			/* not updated; end update */
}

RBT_INIT_IVT(aln, tm_aln_t, tm_aln_rbt_header,
	tm_aln_rbt_cmp,
	tm_aln_ivt_update
);
RBT_INIT_ITER_PATCH(match_aln, tm_aln_t, tm_aln_rbt_header,
	tm_aln_rbt_cmp_head,
	tm_aln_rbt_cmp_tail,
	tm_aln_rbt_cmp_iter,
	tm_aln_ivt_update
);
RBT_INIT_ITER(isct_aln, tm_aln_t, tm_aln_rbt_header,
	tm_aln_ivt_cmp_head,
	tm_aln_ivt_cmp_tail,
	tm_aln_ivt_cmp_iter
);


/* alignment dedup hash; we put alignment end position in this bin */
typedef struct {
	uint64_t key;
	union {
		struct {
			uint32_t untouched;
			uint32_t rid;			/* record reference-side id */
		} s;
		uint64_t val;
	} u;
} tm_dedup_t;

enum tm_dedup_state_e {
	NOT_EVALD        = 0,
	EVALD_BUT_FAILED = 1,
	EVALD_AND_SUCCD  = 2,
};

#define tm_dedup_key(_p)				(_p)->key
#define tm_dedup_val(_p)				(_p)->u.val
#define MM_HASH_INIT_VAL				( RH_INIT_VAL )
RH_INIT(dedup, tm_dedup_t,
	uint64_t, tm_dedup_key,
	uint64_t, tm_dedup_val
);


/* working buffers */
struct tm_scan_s {
	struct {
		tm_seed_v arr;		/* for sorting seeds */
		tm_sqiv_v sqiv;
	} seed;

	struct {
		tm_chain_v arr;
	} chain;

	struct {
		uint32_t max_weight;	/* working variable */
		dz_arena_t *fill, *trace;

		/* dedup hash */
		rh_dedup_t pos;

		/* result bin */
		tm_aln_v arr;
	} extend;
};


// static _force_inline
tm_scan_t *tm_scan_init(void)
{
	tm_scan_t *self = malloc(sizeof(tm_scan_t));
	memset(self, 0, sizeof(tm_scan_t));

	/* init hash */
	rh_init_static_dedup(&self->extend.pos, 0);

	/* clear DP working space */
	size_t const size = DZ_MEM_INIT_SIZE - _roundup(sizeof(dz_arena_t), DZ_MEM_ALIGN_SIZE);
	self->extend.fill = dz_arena_create(size);
	self->extend.trace = NULL;
	return(self);
}

// static _force_inline
void tm_scan_destroy(tm_scan_t *self)
{
	kv_destroy(self->seed.arr);
	kv_destroy(self->seed.sqiv);
	kv_destroy(self->chain.arr);
	kv_destroy(self->extend.arr);

	/* extend */
	rh_destroy_static_dedup(&self->extend.pos);
	dz_arena_destroy(self->extend.fill);

	free(self);
	return;
}

static _force_inline
size_t tm_scan_clip_size(size_t size)
{
	if(size == 0) { return(0); }
	size_t const clipped = 0x8000000000000000>>_lzc_u64(size - 1);
	return(clipped);
}

// static _force_inline
void tm_scan_start_batch(tm_scan_t *self, dz_arena_t *mem)
{
	/* update mem arena */
	self->extend.trace = mem;
	dz_arena_pop_stack(self->extend.fill);

	/* skip if buffer is empty */
	if(kv_ptr(self->extend.arr) == NULL) { return; }

	/* sort of garbage collection */
	size_t const seed = tm_scan_clip_size(kv_max(self->seed.arr));
	size_t const sqiv = tm_scan_clip_size(kv_max(self->seed.sqiv));
	size_t const chain = tm_scan_clip_size(kv_max(self->chain.arr));
	size_t const aln = tm_scan_clip_size(kv_max(self->extend.arr));
	size_t const dedup = tm_scan_clip_size(rh_max_dedup(&self->extend.pos));

	kv_resize(tm_seed_t, self->seed.arr, seed);
	kv_resize(tm_sqiv_t, self->seed.arr, sqiv);
	kv_resize(tm_chain_t, self->chain.arr, chain);
	kv_resize(tm_aln_t, self->extend.arr, aln);
	rh_resize_dedup(&self->extend.pos, dedup);
	return;
}

static _force_inline
void tm_scan_clear(tm_scan_t *self)
{
	kv_clear(self->seed.arr);
	kv_clear(self->seed.sqiv);
	kv_clear(self->chain.arr);
	kv_clear(self->extend.arr);
	return;
}


static _force_inline
tm_seed_t *tm_seed_reserve(tm_seed_v *buf, tm_seed_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	size_t const scnt = kv_cnt(*buf);
	tm_seed_t *p = kv_reserve(tm_seed_t, *buf, scnt + 65536);

	// debug("scnt(%zu), q(%p), p(%p, %p)", scnt, q, p, &p[scnt]);
	return(&p[scnt]);
}

static _force_inline
size_t tm_expand_seed(tm_ref_state_t s, v4i32_t uofs, v4i32_t vofs, tm_seed_t *q)
{
	uint32_t const mask = 0xc000ffff;
	tm_ref_match_t m = tm_ref_get_arr(s);

	#if 1
	v4i32_t const wmask = _set_v4i32(mask);		/* leave lower 15bits and direction flag */
	for(size_t i = 0; i < m.cnt; i += 4) {
		/*
		 * u = 2 * r  - q'
 		 * v = 2 * q' - r
 		 */
		v4i32_t const z = ({
			__m128i const x = _mm_loadl_epi64((__m128i const *)&m.ptr[i]);	/* movq (mem), xmm */
			__m128i const y = _mm_cvtepi16_epi32(x);	/* pmovsxwd; sign expansion (negative for reverse) */
			((v4i32_t){ .v1 = y });
		});
		v4i32_t const w = _and_v4i32(z, wmask);			/* disjoin forward and reverse of the reference sequence */

		v4i32_t const u = _sub_v4i32(uofs, w);					/* 2 * q - r + 0x10000 */
		v4i32_t const v = _add_v4i32(vofs, _add_v4i32(w, w));	/* 2 * r - q - 0x20000 */

		/* save: u in upper and v in lower */
		_storeu_v4i32(&q[i],     _lo_v4i32(u, v));
		_storeu_v4i32(&q[i + 2], _hi_v4i32(u, v));
	}
	#else
	/* serial; has a bug */
	for(size_t i = 0; i < m.cnt; i++) {
		int16_t const z = m.ptr[i];
		uint32_t const w = ((int32_t)z) & mask;			/* signed expansion */
		uint32_t const u = _ext_v4i32(uofs, 0) - w;
		uint32_t const v = _ext_v4i32(vofs, 0) + 2 * w;
		q[i].v = v;
		q[i].u = u;

		_scan_memory(q, sizeof(tm_seed_t));				/* valgrind use-of-uninitialized-value */
	}
	#endif
	return(m.cnt);
}


static _force_inline
tm_sqiv_t *tm_sqiv_reserve(tm_sqiv_v *buf, tm_sqiv_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	size_t const vcnt = kv_cnt(*buf);
	tm_sqiv_t *p = kv_reserve(tm_sqiv_t, *buf, vcnt + 64);
	return(&p[vcnt]);
}

static _force_inline
size_t tm_save_sqiv(tm_ref_squash_t sq, tm_sqiv_t *r)
{
	/* just copy */
	_storeu_v2i32(r, sq.v);
	return(1);
}

static _force_inline
size_t tm_save_last_sqiv(tm_ref_state_t s, tm_sqiv_t *r)
{
	tm_ref_squash_t sq = tm_ref_load_squash(s.bin, -1LL, 0);	/* unmatching = -1 to insert zeros to lower 32bit */
	_storeu_v2i32(r, sq.v);
	return(1 + s.unmatching);
}

static _force_inline
size_t tm_collect_seed(tm_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, size_t intv, tm_seed_v *seed, tm_sqiv_v *sqiv)
{
	/* load coordinate constants */
	v4i32_t const uinc = _set_v4i32(2), vinc = _set_v4i32(-1);
	v4i32_t uofs = _set_v4i32(0x10000 + 2 * (TM_SEED_POS_OFS + 1));
	v4i32_t vofs = _set_v4i32(0x80000000 - 0x20000 - (TM_SEED_POS_OFS + 1));	/* query offsets */

	tm_seed_t *q = tm_seed_reserve(seed, kv_ptr(*seed));
	tm_sqiv_t *r = tm_sqiv_reserve(sqiv, kv_ptr(*sqiv));

	size_t const mask = intv - 1;

	/* initial state (of three general-purpose registers) */
	tm_ref_state_t s = tm_ref_match_init(ref);

	/* for each base... */
	uint8_t const *t = &query[qlen];
	size_t rem = qlen + 1;
	while(--rem != 0) {
		/* update coordinates */
		uofs = _add_v4i32(uofs, uinc);
		vofs = _add_v4i32(vofs, vinc);

		/* update matching status for the next bin (prefetch) */
		tm_ref_next_t n = tm_ref_match_next(s, (rem & mask) == 0, t[-rem]);

		/* save matching state of the previous bin */
		r += tm_save_sqiv(n.squash, r) + s.unmatching;	/* do not forward pointer if previous bin did not match */

		/* skip current bin if not matching */
		s = n.state;
		if(s.unmatching == 0) {
			q += tm_expand_seed(s, uofs, vofs, q);
		}

		/* test array expansion every 32 times */
		if((rem & 0x1f) != 0) { continue; }
		q = tm_seed_reserve(seed, q);
		r = tm_sqiv_reserve(sqiv, r);
	}

	/* save the last sqiv */
	r += tm_save_last_sqiv(s, r);

	/* update seed count */
	size_t const scnt = q - kv_ptr(*seed);
	size_t const vcnt = r - kv_ptr(*sqiv);
	kv_cnt(*seed) = scnt;
	kv_cnt(*sqiv) = vcnt;
	return(scnt);
}

static _force_inline
size_t tm_squash_seed(tm_sqiv_t const *sqiv, size_t icnt, tm_seed_v *seed)
{
	tm_seed_t *dp = kv_ptr(*seed);
	tm_seed_t const *sp = dp;
	// debug("sqiv(%p), icnt(%zu)", sqiv, icnt);

	/* squash succeeding match */
	for(size_t i = 0; i < icnt; i++) {
		_scan_memory(&sqiv[i], sizeof(tm_sqiv_t));		/* valgrind use-of-uninitialized-value */
		tm_sqiv_t sq = {
			.dst = sqiv[i].dst & 0xfff,
			.src = sqiv[i].src & 0xfff,
			.cnt = sqiv[i].cnt & 0xfff
		};
		// debug("spos(%zu), dpos(%zu), cnt(%u, %u, %u)", sp - kv_ptr(*seed), dp - kv_ptr(*seed), sq.dst, sq.src, sq.cnt);

		/* former half */
		for(size_t j = 0; j < sq.dst; j += 4) {
			v32i8_t const v = _loadu_v32i8(&sp[j]);
			_storeu_v32i8(&dp[j], v);
		}

		/* latter half */
		dp -= sq.src - sq.dst;
		for(size_t j = sq.src; j < sq.cnt; j += 4) {
			v32i8_t const v = _loadu_v32i8(&sp[j]);
			_storeu_v32i8(&dp[j], v);
		}

		/* forward pointers */
		sp += sq.cnt;
		dp += sq.cnt;
		// debug("(%u, %u, %u), idx(%zu, %zu), cnt(%zu), squashed(%zu)", sq.d, sq.s, sq.c, dp - kv_ptr(*seed), sp - kv_ptr(*seed), sq.c, sq.s - sq.d);
	}

	// debug("sp(%zu), scnt(%zu)", sp - kv_ptr(*seed), kv_cnt(*seed));
	debugblock({
		if((size_t)(sp - kv_ptr(*seed)) != kv_cnt(*seed)) { trap(); }
	});

	size_t const cnt = dp - kv_ptr(*seed);
	kv_cnt(*seed) = cnt;
	return(cnt);
}

static _force_inline
int64_t tm_chain_test_ptr(tm_seed_t const *p, tm_seed_t const *t)
{
	// debug("check tail, test(%ld)", (int64_t)((ptrdiff_t)(t - p - 1)));
	return((int64_t)((ptrdiff_t)(t - p - 1)));
}

static _force_inline
tm_seed_t *tm_chain_find_first(uint64_t ruv, tm_seed_t *p, tm_seed_t *t, uint64_t lb, uint64_t window)
{
	/* constants */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const tmask = 0x8000000080000000;	/* extract two sign bit pair */

	/* load bounds */
	uint64_t const uv = ruv;					/* (ulb, vlb) */
	uint64_t ub = (uv & umask) + window;		/* (uub, vub - vlb) */
	// debug("uub(%lx), vwlen(%lx)", ub>>32, ub & 0xffffffff);

	int64_t cont = 0;		/* continuous flag */
	while(1) {
		/* check if reached tail */
		if((cont | tm_chain_test_ptr(++p, t)) < 0) { return(NULL); }

		uint64_t const v = _loadu_u64(p) - lb;	/* 1: (u, v - vlb) for inclusion test */
		// debug("test inclusion, pv(%lx), pu(%lx), u(%lx), {v - vlb}(%lx), test(%lx)", p->v, p->u, v>>32, v & 0xffffffff, (ub - v) & tmask);

		cont = ub - v;		/* save diff for testing MSb */
		if(((cont | v) & tmask) == 0) {			/* 2,3: break if chainable (first chainable seed found) */
			// debug("chainable");
			break;
		}

		/* unchainable */
		// debug("unchainable");
	}
	return(p);
}

static _force_inline
tm_seed_t *tm_chain_find_alt(tm_seed_t *p, tm_seed_t *t, uint64_t lb)
{
	/* constants */
	uint64_t const tmask = 0x8000000080000000;	/* extract two sign bit pair */

	/* calculate bounds */
	uint64_t rv = _loadu_u64(p) - lb;
	uint64_t ub = rv + (rv<<32);
	// debug("u(%lx), {v - vlb}(%lx), uub(%lx), vub(%lx)", rv>>32, rv & 0xffffffff, ub>>32, ub & 0xffffffff);

	/* keep nearest seed */
	tm_seed_t *n = p;

	int64_t cont = 0;
	while((cont | tm_chain_test_ptr(++p, t)) >= 0) {
		uint64_t const v = _loadu_u64(p) - lb;	/* 1: (u, v - vlb) for inclusion test */
		// debug("u(%lx), {v - vlb}(%lx), uub(%lx), vub(%lx), test(%lx), term(%lx)", p->u, v & 0xffffffff, ub>>32, ub & 0xffffffff, (ub - v) & tmask, (int64_t)(ub - v) < 0);
		cont = ub - v;
		if((cont | v) & tmask) {		/* skip if unchainable */
			// debug("unchainable");
			continue;
		}

		/* chainable; test if the seed is nearer than the previous */
		uint64_t const w = v + (v<<32);			/* 2,3: (u + v - vlb, v - vlb) */
		// debug("{u + v - vlb}(%lx), {v - vlb}(%lx)", w>>32, w & 0xffffffff);
		if((ub - w) & tmask) {			/* 4,5: further than previous */
			// debug("further; terminate");
			break;
		}

		/* nearer seed found */
		ub = w;			/* update bounds */
		n  = p;			/* save pointer */
	}
	return(n);
}

static _force_inline
tm_seed_t *tm_chain_find_nearest(uint64_t ruv, tm_seed_t *p, tm_seed_t *t, uint64_t window)
{
	/* load root positions */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const uv = ruv;					/* (ulb, vlb) */
	uint64_t lb = (uv & ~umask);				/* (0, vlb) */
	// debug("ulb(%lx), vlb(%lx)", uv>>32, lb);

	/* find first chainable seed */
	tm_seed_t *n = tm_chain_find_first(ruv, p, t, lb, window);
	if(n == NULL) { return(NULL); }				/* we expect control paths are kept isolated */

	/* then inspect alternative chainable seeds */
	return(tm_chain_find_alt(n, t, lb));
}

#define _print_v2i32x(a) { \
	debug("(v2i32_t) %s(%x, %x)", #a, _ext_v2i32(a, 1), _ext_v2i32(a, 0)); \
}
#define _print_v4i32x(a) { \
	debug("(v4i32_t) %s(%x, %x, %x, %x)", #a, _ext_v4i32(a, 3), _ext_v4i32(a, 2), _ext_v4i32(a, 1), _ext_v4i32(a, 0)); \
}

/*
 * uv -> xy coordinate conversion: epos[2] must be zero
 */
static _force_inline
v4i32_t tm_chain_compose(v2i32_t spos, v2i32_t epos)
{
	v2i64_t const coef = _seta_v2i64(0x55555556, 0x55555556);	/* 0.3333... + roundup alpha */
	v4i32_t const ofs  = _seta_v4i32(
		0, 0x80000000 - 3 * TM_SEED_POS_OFS,	/* -3 * base_offset */
		0, 0x00030000 - 3 * TM_SEED_POS_OFS		/* 3 * direction flag - 3 * base_offset */
	);

	v4i32_t const a = _lo_v4i32(epos, spos);	/* 1; (ue, us, ve, vs) */
	v4i32_t const b = _add_v4i32(				/* 3; (ye, ys, xe, xs) = (2*ue + ve, 2*us + vs, 2*ve + ue, 2*vs + us) */
		_add_v4i32(a, a),						/* 2; (2*ue, 2*us, 2*ve, 2*vs) */
		_shuf_v4i32(a, _km_v4i32(1, 0, 3, 2))	/* 2; (ve, vs, ue, us) */
	);

	/* decode span */
	v4i32_t const d = (v4i32_t){
		_mm_srli_epi64(b.v1, 32)				/* 4; (0, ye, 0, xe) */
	};
	v4i32_t const c = _sub_v4i32(d, b);			/* 5; (-ye, ye - ys, -xe, xe - xs) */
	v4i32_t const span = (v4i32_t){
		_mm_mul_epu32(coef.v1, c.v1)			/* 6-; (qspan, -, rspan, -) */
	};
	// _print_v4i32x(span);

	/* decode pos */
	v4i32_t const e = _bsr_v4i32(c, 1);			/* 6; (-, -ye, -, -xe) */
	v4i32_t const f = (v4i32_t){
		_mm_blendv_epi32(b.v1, e.v1, epos.v1)	/* 7; (-, ys, -, dir ? -xe : xs) */
	};
	v4i32_t const g = _add_v4i32(f, ofs);		/* 8; (-, ys + 0x80000000, -, (dir ? -xe : xs) + 0xc0030000) */
	v4i32_t const pos = (v4i32_t){
		_mm_mul_epu32(coef.v1, g.v1)			/* 9-; (qpos, 0, rpos, -) */
	};
	// _print_v4i32x(pos);

	/* fold pos and span */
	v4i32_t const h = (v4i32_t){
		_mm_shuffle2_epi32(pos.v1, span.v1, 0xdd)		/* 15; (qspan, rspan, qpos, rpos) */
	};
	// _print_v4i32x(h);
	return(h);		/* (qspan, rspan, qpos, rpos) */
}

static _force_inline
tm_chain_t *tm_chain_reserve(tm_chain_v *chain, size_t scnt)
{
	size_t const used = kv_cnt(*chain);
	tm_chain_t *base = kv_reserve(tm_chain_t, *chain, used + scnt);
	return(&base[used]);
}

static _force_inline
size_t tm_chain_finalize(tm_chain_v *chain, tm_chain_t *q)
{
	tm_chain_t *base = kv_ptr(*chain);
	size_t const used = kv_cnt(*chain);
	size_t const ccnt = q - &base[used];
	// debug("cnt(%zu, %zu), q(%p, %p)", ccnt, used, q, base);
	kv_cnt(*chain) = used + ccnt;
	return(ccnt);
}

static _force_inline
size_t tm_chain_seed(tm_idx_profile_t const *profile, tm_seed_t *seed, size_t scnt, tm_chain_v *chain)
{
	/*
	 * greedy extension; relatively smaller window size (~32) is preferred for good performance.
	 * do all of the chainablity test on general-purpose register for smaller latency of core loop.
	 * keep seed and counter on xmm regsters and do everything of the chain composing on xmm registers
	 * to minimize spill of general-purpose registers.
	 */
	uint64_t const chained = 0x80000000;				/* offset */
	uint64_t const window  = profile->chain.window.all;	/* window sizes */

	v4i32_t const inc = _seta_v4i32(0, 1, 0, 0);		/* count #seeds */
	v4i32_t const adj = _load_v4i32(&profile->chain.kadj[0]);	/* (-, min_scnt, k, k) */

	/* src pointers */
	tm_seed_t *p = seed - 1, *t = &seed[scnt];

	/* reserve mem for chain */
	tm_chain_t *q = tm_chain_reserve(chain, scnt);
	// debug("scnt(%zu)", scnt);

	/* for each seed */
	while(++p < t) {
		/* skip chained seed to find next root */
		// debug("u(%u), chained(%lu)", p->u, p->u >= chained);
		if(p->u >= chained) { continue; }

		/* root found; keep seed on xmm register */
		v2i32_t const root = _load_v2i32(p);	/* we expect this won't be spilled */
		v4i32_t spos = _sub_v4i32((v4i32_t){ root.v1 }, adj);	/* (-, -min_scnt, vspos - k, uspos - k) */

		/* iteratively link nearest seed */
		tm_seed_t *s = p;
		uint64_t ruv = _loadu_u64(s);
		while(1) {
			tm_seed_t *n = tm_chain_find_nearest(ruv, s, t, window);
			if(n == NULL) { break; }

			/* increment seed count */
			spos = _add_v4i32(spos, inc);

			s = n;					/* save last chained seed */
			ruv = _loadu_u64(s);	/* before flag set */
			s->u += chained;		/* mark chained */
		}

		/*
		debugblock({
			if((int32_t)_ext_v2i32(spos, 2) >= 0) {
				debug("%r ---> %r, cnt(%d)", tm_seed_to_str, p, tm_seed_to_str, s, _ext_v2i32(spos, 2));
			}
		});
		*/

		v2i32_t const epos = (v2i32_t){ .v1 = _mm_cvtsi64_si128(ruv) };
		if((int32_t)_ext_v2i32(spos, 2) < 0) { continue; }		/* skip if too short */

		/* save */
		v4i32_t const c = tm_chain_compose((v2i32_t){ spos.v1 }, epos);
		_store_v4i32(q, c);
		q++;
	}

	/* update chain count */
	return(tm_chain_finalize(chain, q));
}



/* filter */

typedef struct {
	v16i8_t scv, gv, cv, pv;
	uint8_t fw[64], rv[64];
} tm_filter_work_t;


static _force_inline
void tm_filter_work_init(tm_filter_work_t *w, tm_idx_profile_t const *profile)
{
	w->scv = _load_v16i8(&profile->filter.score_matrix); 
	w->gv  = _load_v16i8(&profile->filter.gap);
	w->pv  = _load_v16i8(&profile->filter.init);

	/* derive p = 1 vector */
	w->cv  = _add_v16i8(w->pv, w->gv);
	w->cv  = _maxu_v16i8(w->cv, _bsr_v16i8(w->cv, 1));
	return;
}


/* 4-bit encoding for reference side; return in reversed orientation */
static _force_inline
v32i8_t tm_filter_load_rf(uint8_t const *p) {
	v32i8_t const x = _loadu_v32i8(p);
	v32i8_t const y = _swap_v32i8(x);

	// _print_v32i8x(y);
	return(y);
}

static _force_inline
v32i8_t tm_filter_load_rr(uint8_t const *p) {
	v32i8_t const x = _loadu_v32i8(p - 32);

	// _print_v32i8x(x);
	return(x);
}


/* 2-bit encoding at [3:2] for query side */
static _force_inline
v32i8_t tm_filter_load_qf(uint8_t const *p) {
	return(_loadu_v32i8(p));
}

static _force_inline
v32i8_t tm_filter_load_qr(uint8_t const *p) {
	v32i8_t const v = _loadu_v32i8(p - 32);
	return(_swap_v32i8(v));
}

static _force_inline
v32i8_t tm_filter_conv_r(v32i8_t v)
{
	static uint8_t const rconv[16] __attribute__(( aligned(16) )) = {
		[tA] = A, [tC] = C, [tG] = G, [tT] = T,
		[tR] = R, [tY] = Y, [tS] = S,	/* tS ^ 0x0f and tW ^ 0x0f never happen since the input sequence is always canonical */
		[tK] = K, [tM] = M, [tW] = W,
		[tB] = B, [tD] = D, [tH] = H, [tV] = V
	};
	v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(rconv));
	v32i8_t const x = _shr_v32i8(v, 1);
	return(_shuf_v32i8(cv, x));
}

static _force_inline
v32i8_t tm_filter_conv_q(v32i8_t v, uint64_t dir)
{
	/* convert to 4bit */
	static uint8_t const qconv[16] __attribute__(( aligned(16) )) = {
		[nA]     = A, [nC]     = C, [nG]     = G, [nT]     = T,	/* forward */
		[nA + 1] = T, [nC + 1] = G, [nG + 1] = C, [nT + 1] = A	/* reverse-complemented */
	};
	v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(qconv));
	v32i8_t const x = _set_v32i8(dir);
	v32i8_t const y = _add_v32i8(v, x);
	return(_shuf_v32i8(cv, y));
}

static _force_inline
void tm_filter_load_seq(tm_filter_work_t *w, uint8_t const *r, uint8_t const *q, tm_chain_raw_t const *c)
{
	/* load positions */
	uint64_t const dir = c->pos.dir;
	tm_pair_t const pos  = tm_chain_raw_pos(c);
	tm_pair_t const span = tm_chain_raw_span(c);
	// debug("rr(%u), rf(%u), qr(%u), qf(%u), dir(%x)", pos.r, pos.r + span.r, pos.q, pos.q + span.q, dir);

	/* ref; convert before save */
	v32i8_t const rf = tm_filter_conv_r(tm_filter_load_rf(&r[pos.r + span.r]));
	v32i8_t const rr = tm_filter_conv_r(tm_filter_load_rr(&r[pos.r]));
	_store_v32i8(&w->fw[0], rf);
	_store_v32i8(&w->rv[0], rr);

	/* query; already converted */
	v32i8_t const qf = tm_filter_conv_q(tm_filter_load_qf(&q[pos.q + span.q]), dir);
	v32i8_t const qr = tm_filter_conv_q(tm_filter_load_qr(&q[pos.q]), dir);

	uint8_t *qfp = (dir & 0x01) ? &w->rv[32] : &w->fw[32];
	uint8_t *qrp = (dir & 0x01) ? &w->fw[32] : &w->rv[32];
	_store_v32i8(qfp, qf);
	_store_v32i8(qrp, qr);

	/*
	debugblock({
		v32i8_t const sf = _swap_v32i8(rf);
		v32i8_t const sr = _swap_v32i8(rr);
		v32i8_t const tf = _load_v32i8(qfp);
		v32i8_t const tr = _load_v32i8(qrp);
		_print_v32i8(sf);
		_print_v32i8(tf);
		_print_v32i8(sr);
		_print_v32i8(tr);
	});
	*/
	return;
}

static _force_inline
int64_t tm_filter_extend_core(tm_filter_work_t *w, uint8_t const *buf)
{
	uint8_t const *r = &buf[32 - 8], *q = &buf[32 - 8];
	v16i8_t const scv = w->scv, gv = w->gv;
	v16i8_t pv = w->pv, cv = w->cv, mv = w->pv;

	#define _calc_next(_rv, _qv, _pv, _cv, _shift) ({ \
		/* calc score profile vector with shuffling score matrix vector */ \
		v16i8_t const match = _and_v16i8(_rv, _qv);		/* cmpeq */ \
		v16i8_t const sv = _shuf_v16i8(scv, match); \
		v16i8_t const x = _add_v16i8(_pv, sv); \
		v16i8_t const y = _add_v16i8(_cv, gv); \
		v16i8_t const z = _maxu_v16i8(y, _shift(y, 1)); \
		v16i8_t const n = _maxu_v16i8(x, z); \
		n; \
	})

	/* 32 vector updates */
	v16i8_t qv = _loadu_v16i8(q);		/* load query side initial vector */
	for(size_t i = 0; i < 16; i++) {
		/* #0 (even) */
		r--;
		v16i8_t const rv = _loadu_v16i8(r);	/* move rightward */
		v16i8_t const ev = _calc_next(rv, qv, pv, cv, _bsl_v16i8);
		mv = _maxu_v16i8(mv, ev);
		// _print_v16u8(ev);
		// _print_v16u8(mv);

		/* #1 (odd) */
		q++;
		qv = _loadu_v16i8(q);	/* move downward */
		v16i8_t const ov = _calc_next(rv, qv, cv, ev, _bsr_v16i8);
		mv = _maxu_v16i8(mv, ov);
		// _print_v16u8(ov);
		// _print_v16u8(mv);

		/* alias registers for the next loop */
		pv = ev;
		cv = ov;
	}

	/* fold scores */
	return(_hmaxu_v16i8(mv) - 128LL);

	#undef _calc_next
}

static _force_inline
int64_t tm_filter_extend(tm_filter_work_t *w, tm_idx_profile_t const *profile, uint8_t const *r, uint8_t const *q, tm_chain_raw_t const *c)
{
	/* skip if long enough */
	uint32_t const qth = profile->filter.span_thresh + 1;
	if(c->span.q >= qth) { return(1); }

	/* load sequences */
	tm_filter_load_seq(w, r, q, c);

	/* extend */
	int64_t const fw = tm_filter_extend_core(w, w->fw);
	int64_t const rv = tm_filter_extend_core(w, w->rv);
	// debug("fw(%ld), rv(%ld), score(%ld)", fw, rv, fw + rv);
	return(fw + rv >= profile->filter.min_score);
}

#if 1
static _force_inline
v2i32_t tm_filter_calc_pos(tm_chain_raw_t const *p)
{
	/* load */
	v2i32_t const spos = _loadu_v2i32(&p->pos.r);
	v2i32_t const span = _loadu_v2i32(&p->span.r);

	/* create mask: (0, dir ? -1 : 0) */
	v2i32_t const thresh = _seta_v2i32(0x3fffffff, 0xffff);
	v2i32_t const mask   = _gt_v2i32(spos, thresh);

	/* (qpos, rpos) for forward, (qpos, rpos + rspan) for reverse */
	v2i32_t const epos = _add_v2i32(spos, _and_v2i32(mask, span));

	// _print_v2i32x(spos);
	// _print_v2i32x(span);
	// _print_v2i32x(mask);
	// _print_v2i32x(epos);
	return(epos);
}
#else
static _force_inline
v2i32_t tm_filter_calc_pos(tm_chain_raw_t const *p)
{
	/* load */
	v2i32_t const spos = _loadu_v2i32(&p->pos.r);
	v2i32_t const span = _loadu_v2i32(&p->span.r);

	/* create mask: (0, dir ? -1 : 0) */
	v2i32_t const thresh = _seta_v2i32(0x3fffffff, 0xffff);
	v2i32_t const mask   = _gt_v2i32(spos, thresh);

	/* (qpos + qspan, rpos + rspan) for forward, (qpos + qspan, rpos) for reverse */
	v2i32_t const epos = _add_v2i32(spos, _andn_v2i32(mask, span));

	_print_v2i32x(spos);
	_print_v2i32x(span);
	_print_v2i32x(mask);
	_print_v2i32x(epos);
	return(epos);
}
#endif

static _force_inline
size_t tm_filter_save_chain(uint32_t rid, tm_chain_t *q, tm_chain_raw_t const *p)
{
	// debug("p(%p), %r", p, tm_chain_raw_to_str, p);

	/* calc weight */
	ZCNT_RESULT size_t weight = _lzc_u32(p->span.r);

	/* adjust span to obtain end position */
	v2i32_t const v = tm_filter_calc_pos(p);
	_storeu_v2i32(q, v);

	q->attr.sep.rid    = rid;		/* overwrite */
	q->attr.sep.weight = weight;	/* negated and offsetted by 32; the larger for the longer chain */

	/* copy span */
	q->span.q = p->span.q;
	// debug("q(%p), %r, rid(%u), weight(%u)", q, tm_chain_to_str, q, rid, weight);
	return(1);
}

static _force_inline
size_t tm_filter_chain(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, tm_chain_t *chain, size_t ccnt)
{
	/* alias pointer */
	tm_chain_raw_t const *src = (tm_chain_raw_t const *)chain;
	tm_chain_t *dst = chain;

	uint32_t const rid = tm_idx_ref_rid(si);		/* FIXME: kept on xmm? */
	uint8_t const *ref = tm_idx_ref_seq_ptr(si);

	/* init working buffer (load score vectors on xmm) */
	tm_filter_work_t w __attribute__(( aligned(32) ));
	tm_filter_work_init(&w, profile);

	// debug("min_score(%d), span_thresh(%u), rid(%u)", profile->filter.min_score, profile->filter.span_thresh, rid);
	for(size_t i = 0; i < ccnt; i++) {
		tm_chain_raw_t const *p = &src[i];
		// debug("i(%zu), ccnt(%zu), %r", i, ccnt, tm_chain_raw_to_str, p);

		/* try short extension if not heavy enough */
		if(!tm_filter_extend(&w, profile, ref, query, p)) { continue; }

		/* copy and fold in reference id */
		dst += tm_filter_save_chain(rid, dst, p);
	}
	// debug("cfilt(%zu)", dst - chain);
	return(dst - chain);
}

static _force_inline
size_t tm_seed_and_sort(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, size_t qlen, tm_seed_v *seed, tm_sqiv_v *sqiv)
{
	/* enumerate seeds */
	kv_clear(*seed);
	kv_clear(*sqiv);
	if(tm_collect_seed(&si->ref, query, qlen, profile->chain.squash_intv, seed, sqiv) == 0) {
		return(0);
	}

	/* squash overlapping seeds */
	if(tm_squash_seed(kv_ptr(*sqiv), kv_cnt(*sqiv), seed) == 0) {
		return(0);
	}

	/* sort seeds by u-coordinates for chaining */
	tm_seed_t *sptr = kv_ptr(*seed);
	size_t const scnt = kv_cnt(*seed);

	radix_sort_seed(sptr, scnt);

	/*
	debug("scnt(%zu)", scnt);
	for(size_t i = 0; i < scnt; i++) {
		debug("i(%zu), %r", i, tm_seed_to_str, &sptr[i]);
	}
	*/
	return(scnt);
}

static _force_inline
size_t tm_chain_and_filter(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, tm_seed_t *seed, size_t scnt, tm_chain_v *chain)
{
	/* chaining; return if no chain obtained */
	size_t const cbase = kv_cnt(*chain);		/* save chain count */
	size_t const ccnt  = tm_chain_seed(profile, seed, scnt, chain);
	if(ccnt == 0) { return(0); }
	// debug("cbase(%zu), ccnt(%zu)", cbase, ccnt);

	/* filter (try small extension with simple score matrix) */
	size_t const cfilt = tm_filter_chain(si, profile, query, &kv_ptr(*chain)[cbase], ccnt);
	
	/* update chain count */
	kv_cnt(*chain) = cbase + cfilt;
	return(cfilt);
}

static _force_inline
size_t tm_collect_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *seq, size_t slen)
{
	tm_idx_sketch_t const **si = (tm_idx_sketch_t const **)idx->sketch.arr;
	tm_idx_profile_t const **pf = (tm_idx_profile_t const **)idx->profile.arr;

	/* clear result bin */
	kv_clear(self->chain.arr);

	/* collect chain for each reference sequence */
	for(size_t i = 0; i < idx->sketch.cnt; i++) {
		/* collect seeds */
		if(tm_seed_and_sort(si[i], pf[si[i]->h.pid], seq, slen, &self->seed.arr, &self->seed.sqiv) == 0) { continue; }

		/* chain and filter */
		tm_seed_t *seed = kv_ptr(self->seed.arr);
		size_t scnt = kv_cnt(self->seed.arr);
		tm_chain_and_filter(si[i], pf[si[i]->h.pid], seq, seed, scnt, &self->chain.arr);

		// fprintf(stderr, "sketch(%zu)\n", i);
	}

	debug("ccnt(%zu)", kv_cnt(self->chain.arr));
	return(kv_cnt(self->chain.arr));
}



/* X-drop DP extension */
#define TM_WEIGHT_MARGIN			( 3U )

/* clear everything; called once every query */
static _force_inline
void tm_extend_clear(tm_scan_t *self)
{
	/* clear test threshold */
	self->extend.max_weight = 32;

	/* clear hash */
	rh_clear_dedup(&self->extend.pos);

	/* clear rbt */
	kv_clear(self->extend.arr);
	kv_pushp(tm_aln_t, self->extend.arr);
	rbt_init_static_aln(kv_ptr(self->extend.arr));

	/* we must not clear trace stack because we have saved alignments and metadata in it */
	// dz_arena_flush(self->extend.trace);
	return;
}


static _force_inline
tm_aln_t tm_chain_as_aln(tm_chain_t const *c)
{
	return((tm_aln_t){
		.pos  = tm_chain_pos(c),
		.span = { .q = c->span.q, .r = c->span.q },	/* duplicate qspan because rspan is missing in chain object */
		.qmax = 0
	});
}


/* start position hash */
static _force_inline
uint64_t tm_extend_hash_pos(uint32_t rid, tm_pos_t spos)
{
	/* FIXME: better performance */
	tm_pair_t const p = tm_pos_to_pair(spos);

	uint64_t const x = _loadu_u64(&p);			/* r in lower 32bit, q in upper 32bit */
	uint64_t const y = (rid<<1) + spos.dir;		/* only the lowest bit matters for dir */

	uint64_t const magic = 0xf324a24a1111ULL;
	return((0x10001001 * x) ^ x ^ (x>>31) ^ (x>>18) ^ (magic * y));
}

static _force_inline
uint64_t tm_extend_mark_pos(tm_scan_t *self, uint32_t rid, tm_pos_t spos)
{
	/* duplicated if state is not zero */
	uint64_t const h = tm_extend_hash_pos(rid, spos);
	tm_dedup_t *bin = rh_put_ptr_dedup(&self->extend.pos, h);

	/* we don't expect bin be NULL but sometimes happen, or already evaluated (we don't mind the last state) */
	if(bin == NULL || bin->u.s.untouched == 0ULL) {
		debug("already evaluated, bin(%p), untouched(%u), rid(%u, %u), h(%lx, %lx)", bin, bin->u.s.untouched, bin->u.s.rid, rid, bin->key, h);
		return(1);
	}

	/* not duplicated; put current pos */
	bin->u.s.untouched = 0ULL;
	bin->u.s.rid       = rid;
	// rh_put_dedup(&self->extend.pos, h, 0ULL);
	debug("found new bin(%p), h(%lx)", bin, h);
	return(0);
}


/* alignment bin (alignment collection); sorted by its span */

static _force_inline
uint64_t tm_extend_is_complete(tm_scan_t *self)
{
	_unused(self);

	/* FIXME */
	return(0);
}

static _force_inline
uint64_t tm_extend_test_range(tm_chain_t const *c, tm_aln_t const *p)
{
	/*
	 * test q-range; skip if entirely covered
	 * c->pos.q > p->pos.q && c->pos.q + c->span.q < p->pos.q + p->span.q
	 */
	tm_pair_t const pos = tm_chain_head(c);
	if((uint32_t)(pos.q - p->pos.q) < (uint32_t)(p->span.q - c->span.q)) {
		debug("covered, (%u, %u), (%u, %u)", pos.q, pos.q + c->span.q, p->pos.q, p->pos.q + p->span.q);

		if(c->attr.sep.weight + TM_WEIGHT_MARGIN > p->attr.max_weight + 1U) {
			debug("weight not enough, weight(%u), max_weight(%u)", c->attr.sep.weight, p->attr.max_weight);
			return(1);	/* still too shorter than overlapping alignment */
		}
	}
	return(0);
}

static _force_inline
uint64_t tm_extend_is_covered(tm_scan_t *self, tm_chain_t const *c)
{
	return(0);			/* FIXME: removed */
	debug("%r", tm_chain_to_str, c);

	/* for each overlapping alignment */
	tm_aln_t const *v = kv_ptr(self->extend.arr);
	if(v == NULL) { return(0); }	/* no result found */

	/* create anchor; convert chain to alignment object so that it can be compared to existing alignments */
	tm_aln_t const caln = tm_chain_as_aln(c);		/* convert chain position to (pseudo) alignment range */

	/* init iterator */
	rbt_iter_t it;
	rbt_init_iter_isct_aln(&it);

	tm_aln_t const *p = rbt_fetch_isct_aln(&it, v, &caln);
	while(p != NULL) {
		debug("p(%p), weight(%u), (%u) --- %u --> (%u)", p, p->attr.max_weight, p->pos.q, p->span.q, p->pos.q + p->span.q);

		/* overlapping and heavy enough */
		if(c->attr.sep.weight > p->attr.max_weight) {	/* smaller weight is better */
			debug("weight too large, weight(%u), max_weight(%u)", c->attr.sep.weight, p->attr.max_weight);
			return(1);		/* already covered */
		}

		if(tm_extend_test_range(c, p)) {
			return(1);		/* still too shorter than overlapping alignment */
		}

		p = rbt_fetch_isct_aln(&it, v, &caln);
	}
	return(0);				/* possiblilty remains for better alignment */
}


typedef struct {
	uint64_t collided;
	uint64_t replaced;
} tm_extend_replace_t;


static _force_inline
int64_t tm_extend_compare_aln(tm_aln_t const *x, tm_aln_t const *y)
{
	/* strcmp equivalent */
	return((int64_t)(x->score.raw - y->score.raw));		/* expand sign */
}

static _force_inline
uint64_t tm_extend_patch_bin(tm_scan_t const *self, tm_aln_t *v, rbt_iter_t *it, tm_aln_t *old, tm_aln_t const *new)
{
	_unused(self);

	/* if the old one is larger than the new one, do nothing */
	if(tm_extend_compare_aln(old, new) >= 0) {
		return(0);		/* not replaced (discarded) */
	}

	/* new one is larger; overwrite content except for the header */
	old->qmax = new->pos.q + new->span.q;				/* initialize qmax */
	memcpy(&old->attr, &new->attr, sizeof(tm_aln_t) - offsetof(tm_aln_t, attr));

	/* and update interval tree */
	rbt_patch_match_aln(it, v);
	return(1);			/* replaced */
}

static _force_inline
tm_extend_replace_t tm_extend_slice_bin(tm_scan_t *self, tm_aln_t const *new)
{
	tm_extend_replace_t const notfound = {
		.collided = 0,
		.replaced = 0
	};

	/* apparently no element matches when the result array is empty */
	tm_aln_t *v = kv_ptr(self->extend.arr);
	if(v == NULL) { return(notfound); }

	/* position is compared at once on GP register */
	uint64_t const x = _loadu_u64(&new->pos);	/* uncanonized spos so that it points the tail end of the second extension */

	/* init iterator */
	rbt_iter_t it;
	rbt_init_iter_match_aln(&it);

	while(1) {
		tm_aln_t const *p = rbt_fetch_match_aln(&it, v, new);
		if(p == NULL) { break; }

		/* compare end pos, return iterator if matched */
		uint64_t const y = _loadu_u64(&p->pos);

		// xfprintf(stderr, "compare aln, bid(%zu), eq(%lu), %r, %r\n", p - v, x == y, tm_aln_to_str, p, tm_aln_to_str, new);
		if(x == y) {		/* direction contained */
			/* alignment end position collides; take better one */
			// xfprintf(stderr, "duplication found, patch bin, bid(%zu), replace(%lu), %r, %r\n", p - v, tm_extend_compare_aln(p, new) < 0, tm_aln_to_str, p, tm_aln_to_str, new);
			return((tm_extend_replace_t){
				.collided = 1,
				.replaced = tm_extend_patch_bin(self, v, &it, (tm_aln_t *)p, new)
			});
		}
	}

	/* not found; we need to allocate new bin for the alignment */
	return(notfound);
}

static _force_inline
void tm_extend_push_bin(tm_scan_t *self, tm_aln_t const *aln)
{
	debug("allocate new bin, bid(%zu)", kv_cnt(self->extend.arr));

	/* copy */
	tm_aln_t *p = kv_pushp(tm_aln_t, self->extend.arr);
	memcpy(p, aln, sizeof(tm_aln_t));

	/* update tree */
	rbt_insert_aln(kv_ptr(self->extend.arr), kv_cnt(self->extend.arr) - 1);
	// rbt_print_aln(kv_ptr(self->extend.arr), kv_cnt(self->extend.arr), tm_aln_to_str);

	/* test for patch query */
	debugblock({
		debug("patch");
		tm_aln_t *v = kv_ptr(self->extend.arr);
		uint64_t const x = _loadu_u64(&aln->pos);

		rbt_iter_t it;
		rbt_init_iter_match_aln(&it);

		while(1) {
			tm_aln_t const *q = rbt_fetch_match_aln(&it, v, aln);
			if(q == NULL) { break; }

			/* compare end pos, return iterator if matched */
			uint64_t const y = _loadu_u64(&q->pos);
			if(x == y && q->pos.dir == aln->pos.dir) { debug("found"); break; }
		}
		rbt_patch_match_aln(&it, v);

		debug("done");
	});
	return;
}

static _force_inline
uint64_t tm_extend_record(tm_scan_t *self, tm_aln_t const *aln)
{
	/* returns 1 if recorded */
	tm_extend_replace_t const r = tm_extend_slice_bin(self, aln);
	if(r.collided) { return(r.replaced); }

	/* allocate new bin */
	tm_extend_push_bin(self, aln);
	return(1);		/* recorded */
}



static _force_inline
tm_pair_t tm_aln_span(dz_alignment_t const *aln)
{
	return((tm_pair_t){
		.r = aln->ref_length,
		.q = aln->query_length
	});
}

static _force_inline
tm_aln_t tm_extend_compose_aln(tm_scan_t const *self, tm_chain_t const *chain, tm_pos_t epos, dz_alignment_t const *aln, tm_score_t score)
{
	_unused(self);
	debug("compose, aln(%p), path(%s)", aln, aln->path);

	/* calc pos and span */
	tm_pair_t const span = tm_aln_span(aln);
	tm_pos_t const spos  = tm_sub_pos(epos, span);		/* convert to head position */
	// tm_pos_t const pos  = tm_canon_pos(spos, span);		/* spos.r < epos.r whichever direction is */

	// xfprintf(stderr, "compose, chain(%p), %r, aln(%p), %r\n", chain, tm_chain_to_str, chain, aln, tm_pos_to_str, &spos);

	tm_aln_t const a = {
		/* coordinates and scores */
		.qmax = epos.q,		/* spos.q + span.q, */
		.pos  = spos,
		.span = span,
		.score = score,

		/* attributes */
		.attr = {
			.rid        = chain->attr.sep.rid,
			.max_weight = chain->attr.sep.weight + TM_WEIGHT_MARGIN		/* x8 */
		},

		/* save original */
		.aln  = aln
	};
	return(a);
}




/* reference fetcher and query converter for dz */
typedef struct {
	__m128i mv;
	uint8_t const *p;
	ptrdiff_t inc;

	dz_fill_fetch_t next;
} tm_extend_fetcher_t;

static _force_inline
dz_fill_fetch_t tm_extend_fetch_core(uint8_t const *p, ptrdiff_t inc)
{
	uint32_t const x = (uint32_t)(*p);
	uint32_t const y = x ^ (uint32_t)inc;

	return((dz_fill_fetch_t){
		.is_term = x == '\0',
		.rch     = y & 0x1e					/* note: LSb is always cleared */
	});
}

static _force_inline
dz_fill_fetch_t tm_extend_fetch_patch(dz_fill_fetch_t curr, uint32_t is_term)
{
	return((dz_fill_fetch_t){
		.is_term = curr.is_term,
		.rch     = curr.rch + is_term		/* shift to bonus matrix if tail */
	});
}

static _force_inline
void tm_extend_fetcher_init(tm_extend_fetcher_t *self, uint8_t const *ref, uint32_t dir)
{
	self->p    = ref - dir;
	self->inc  = dir ? -1LL : 1LL;
	self->next = tm_extend_fetch_core(self->p, self->inc);

/*
	uint8_t const *p = ref - dir;
	while(*p != '\0') {
		fprintf(stderr, "%c", "NACMGRSVTWYHKDBN"[*p]);
		p += self->inc;
	}
	fprintf(stderr, "\n");
*/
	return;
}

static _force_inline
dz_fill_fetch_t tm_extend_fetch_next(tm_extend_fetcher_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	_unused(query);

	dz_fill_fetch_t const prev = self->next;

	/* do nothing if reached tail */
	if(prev.is_term) { return(prev); }

	/* fetch next base */
	dz_fill_fetch_t const next = tm_extend_fetch_core(self->p + self->inc, self->inc);
	dz_fill_fetch_t const curr = tm_extend_fetch_patch(prev, next.is_term);

	/* load score matrix for current char */
	self->mv = _mm_cvtsi64_si128(
		_loadu_u32(&score_matrix[curr.rch * DZ_QUERY_MAT_SIZE / 2])		/* already x2 */
	);

	/* save */
	self->p += self->inc;
	self->next = next;

	// debug("p(%p), inc(%ld), is_term(%u), c(%x, %x)", self->p, self->inc, next.is_term, curr.rch, next.rch);
	// _print_v16i8((v16i8_t){ self->mv });
	return(curr);

	#if 0
	if(*self->p == '\0') {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.rch     = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = (self->conv>>(4 * c)) & 0x0f;
	self->mv = _mm_cvtsi64_si128(_loadu_u32(&score_matrix[e * DZ_QUERY_MAT_SIZE]));	/* <<2 */
	// debug("p(%p), c(%x, %c), e(%x), inc(%ld)", self->p, c, "NACMGRSVTWYHKDBN"[c], e, self->inc);
	// _print_v16i8((v16i8_t){ self->mv });

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.rch     = e
	});
	#endif
}

static _force_inline
__m128i tm_extend_get_profile(tm_extend_fetcher_t *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	_unused(score_matrix);

	uint8_t const *packed = dz_query_packed_array(query);
	__m128i const v = ({
		__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
		__m128i const sc = _mm_shuffle_epi8(self->mv, qv);
		// _print_v16i8((v16i8_t){ qv });
		_mm_cvtepi8_epi16(sc);
	});
	// _print_v16i8((v16i8_t){ v });
	return(v);
}

static _force_inline
dz_trace_match_t tm_extend_get_match(void const *unused, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link)
{
	_unused(unused);

	uint64_t const ch = link->rch;
	uint8_t const *packed = dz_query_packed_array(query);
	int8_t const score = score_matrix[ch * DZ_QUERY_MAT_SIZE / 2 + packed[qidx]];

	// debug("rch(%x), qidx(%zx), c(%x, %x), match(%x)", ch, qidx, packed[qidx], tm_2bit_to_qch(packed[qidx]), tm_qrch_is_match(tm_2bit_to_qch(packed[qidx]), ch & 0x1e));
	return((dz_trace_match_t){
		.score = score,
		.match = tm_qrch_is_match(tm_2bit_to_qch(packed[qidx]), ch & 0x1e)
	});
}

static _force_inline
size_t tm_extend_calc_dim(size_t qlen)
{
	_unused(qlen);
	return(1);
}

static _force_inline
__m128i tm_extend_conv(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	_unused(score_matrix);

	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		[nA]     = 0x00, [nC]     = 0x01, [nG]     = 0x02, [nT]     = 0x03,	/* forward */
		[nA + 1] = 0x03, [nC + 1] = 0x02, [nG + 1] = 0x01, [nT + 1] = 0x00,	/* reverse-complemented */
		[0x0f]   = 0x04		/* invalid */
	};
	v16i8_t const cv = _load_v16i8(conv);
	v16i8_t const d = _set_v16i8(dir);
	v16i8_t const x = _or_v16i8((v16i8_t){ v }, d);
	v16i8_t const y = _shuf_v16i8(cv, x);

	// _print_v16i8((v16i8_t){ v });
	// _print_v16i8(y);
	return(y.v1);
}

static _force_inline
dz_query_t *tm_pack_query_wrap(dz_arena_t *mem, dz_profile_t const *profile, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, uint32_t qdir, tm_pos_t pos)
{
	dz_pack_query_t const pack = {
		.dir      = 0,		/* ignored */
		.invalid  = 0x0f,
		.calc_dim = tm_extend_calc_dim,
		.conv     = tm_extend_conv
	};

	/* calc reference remaining length */
	size_t const rlen = tm_idx_ref_seq_len(sk);
	size_t const rrem = pos.dir ? pos.r : rlen - pos.r;

	if(qdir) {
		size_t const qrem = MIN2(2 * rrem, pos.q);
		debug("reverse: qpos(%u) --- qspan(%u) --> qspos(%u)", pos.q, qrem, pos.q - qrem);

		dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, &query[pos.q - qrem], qrem, &pack);
		dz_pack_query_reverse_core(q, profile, &query[pos.q - qrem], qrem, &pack);
		return(q);
	} else {
		size_t const qrem = MIN2(2 * rrem, qlen - pos.q);
		debug("forward: qpos(%u) --- qspan(%u) --> qepos(%u)", pos.q, qrem, pos.q + qrem);

		dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, &query[pos.q], qrem, &pack);
		dz_pack_query_forward_core(q, profile, &query[pos.q], qrem, &pack);
		return(q);
	}
}

static _force_inline
dz_state_t const *tm_extend_wrap(dz_arena_t *mem, dz_profile_t const *profile, tm_idx_sketch_t const *sk, dz_query_t const *q, tm_pos_t pos)
{
	uint8_t const *ref = tm_idx_ref_seq_ptr(sk);

	/* use 4bit fetcher */
	tm_extend_fetcher_t w __attribute__(( aligned(16) ));
	tm_extend_fetcher_init(&w, &ref[pos.r], pos.dir);
	debug("ref(%p, %p, %zu), rlen(%zu), rdir(%u), rpos(%u)", ref, &ref[pos.r], &ref[pos.r] - ref, tm_idx_ref_seq_len(sk), pos.dir, pos.r);

	dz_fetcher_t fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fetch_next_t)tm_extend_fetch_next,
		.get_profile = (dz_get_profile_t)tm_extend_get_profile,		/* direct conversion */
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, q, &fetcher, (dz_state_t const **)&profile->root, 1));
}

static _force_inline
tm_pos_t tm_calc_max_wrap(dz_profile_t const *profile, dz_query_t const *q, dz_state_t const *r, tm_pos_t rpos)
{
	if(r == NULL || r->max.cap == NULL) {
		return(rpos);
	}

	/* get downward max */
	dz_max_pos_t const s = dz_calc_max_pos_core(profile, q, r);

	/* insert bonus */
	tm_pos_t const pos = {
		.r   = rpos.r,
		.q   = rpos.q,
		.dir = rpos.dir
	};

	/* compose span */
	tm_pair_t const span = {
		.r = s.rpos,
		.q = s.qpos
	};

	debug("abs(%d), bonus(%d)", r->max.score.abs, s.bonus);
	return(tm_add_pos(pos, span));			/* direction copied */
}

#if 1
static _force_inline
tm_pos_t tm_extend_load_rpos(tm_scan_t const *self, tm_idx_profile_t const *pf, tm_chain_t const *c)
{
	_unused(self);
	_unused(pf);

	tm_pos_t const pos = tm_chain_pos(c);		/* head pos with direction */
	return(pos);
}
#else
static _force_inline
tm_pair_t tm_extend_load_rpos(tm_scan_t const *self, tm_idx_profile_t const *pf, tm_chain_t const *c)
{
	_unused(self);

	uint32_t const kadj = pf->chain.kadj[0];
	uint32_t const rdir = c->pos.dir;
	tm_pair_t const epos = tm_chain_head(c);
	return((tm_pair_t){
		.r = epos.r - (rdir ? -kadj : kadj),
		.q = epos.q -  kadj
	});
}
#endif

static _force_inline
uint64_t tm_extend_is_tail_maximal(tm_scan_t const *self, tm_idx_sketch_t const *sk, tm_pos_t epos)
{
	_unused(self);

	uint32_t const rend = epos.dir ? 0 : tm_idx_ref_seq_len(sk);
	return(epos.r == rend);
}

typedef struct {
	tm_pos_t pos;
	uint64_t fail;
} tm_extend_first_t;

static _force_inline
tm_extend_first_t tm_extend_first(tm_scan_t *self, tm_idx_profile_t const *pf, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, tm_pos_t rpos)
{
	debug("forward: rdir(%u), rpos(%u, %u)", rpos.dir, rpos.q, rpos.r);

	/* first extension is reference forward */
	dz_query_t const *qr = tm_pack_query_wrap(self->extend.fill,
		pf->extend.dz, sk,
		query, qlen, 0,
		rpos
	);

	/* extend; reference forward */
	dz_state_t const *r = tm_extend_wrap(self->extend.fill,
		pf->extend.dz, sk, qr, rpos
	);
	tm_pos_t const epos = tm_calc_max_wrap(pf->extend.dz, qr, r, rpos);

	debug("forward: rpos(%u, %u) --- score(%d) --> epos(%u, %u)", rpos.q, rpos.r, r != NULL ? r->max.score.abs : 0, epos.q, epos.r);
	return((tm_extend_first_t){
		.pos  = epos,
		.fail = r->max.score.abs < (tm_extend_is_tail_maximal(self, sk, epos) ? pf->extend.bonus : 0)
	});
}

static _force_inline
dz_alignment_t const *tm_extend_second(tm_scan_t *self, tm_idx_profile_t const *pf, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, tm_pos_t epos)
{
	/* second extension is reference reverse: flip direction */
	epos.dir ^= 0x01;
	dz_query_t const *qf = tm_pack_query_wrap(self->extend.fill,
		pf->extend.dz, sk,
		query, qlen, 1,		/* reverse */
		epos
	);
	dz_state_t const *f = tm_extend_wrap(self->extend.fill,
		pf->extend.dz, sk, qf, epos
	);
	debug("reverse: score(%d)", (f != NULL ? f->max.score.abs : -1));

	/* we don't expect f becomes NULL but might happen. score == 0 indicates there is no significant alignment with length > 0 */
	if(f == NULL || f->max.score.abs == 0) {
		return(NULL);		/* not promising or something is wrong */
	}

	/* traceback */
	dz_getter_t getter = {
		.opaque    = NULL,
		.get_match = tm_extend_get_match
	};
	dz_alignment_t const *aln = dz_trace_core(self->extend.trace,
		pf->extend.dz, qf, &getter, f
	);
	debug("reverse: spos(%u, %u) --- score(%d) --> epos(%u, %u), f(%p), aln(%p)", epos.q - aln->query_length, epos.r + (epos.dir ? aln->ref_length : -aln->ref_length), aln->score, epos.q, epos.r, f, aln);

	// fprintf(stderr, "count(%u, %u, %u, %u)\n", aln->match_count, aln->mismatch_count, aln->ins_count, aln->del_count);
	return(aln);
}

typedef struct {
	tm_pos_t epos;
	dz_alignment_t const *aln;
} tm_extend_res_t;

static _force_inline
tm_extend_res_t tm_extend_core(tm_scan_t *self, tm_idx_profile_t const *pf, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, tm_chain_t const *c)
{
	/* flush working memory */
	dz_arena_flush(self->extend.fill);
	tm_extend_res_t const failed = { .aln  = NULL };

	/* load root position */
	tm_pos_t const rpos = tm_extend_load_rpos(self, pf, c);

	/* reference forward */
	tm_extend_first_t const e = tm_extend_first(self, pf, sk, query, qlen, rpos);
	debug("fail(%lu)", e.fail);
	if(e.fail || tm_extend_mark_pos(self, tm_idx_ref_rid(sk), e.pos)) {
		return(failed);			/* it seems the tail position is already evaluated */
	}

	/* reference reverse */
	dz_alignment_t const *aln = tm_extend_second(self, pf, sk, query, qlen, e.pos);	/* direction flipped internally */
	return((tm_extend_res_t){
		.epos = e.pos,
		.aln  = aln
	});
}


static _force_inline
void tm_extend_count_base_core(size_t *acc, v32i8_t v, uint32_t window)
{
	v32i8_t const x = _shl_v32i8(v, 5);
	v32i8_t const y = _shl_v32i8(v, 4);

	uint64_t const l = ((v32_masku_t){ .mask = _mask_v32i8(x) }).all & ~(uint64_t)window;	/* vpmovmskb, andnq */
	uint64_t const h = ((v32_masku_t){ .mask = _mask_v32i8(y) }).all & ~(uint64_t)window;

	size_t const ccnt = _popc_u64(~h & l);		/* andnq, popcntq */
	size_t const gcnt = _popc_u64(h & ~l);
	size_t const tcnt = _popc_u64(h & l);

	acc[0] += ccnt;
	acc[1] += gcnt;
	acc[2] += tcnt;
	return;
}

static _force_inline
void tm_extend_count_base(size_t *acc, uint8_t const *query, size_t qlen)
{
	debug("qlen(%zu)", qlen);

	/* body */
	size_t qpos = 0;
	while((qpos += 32) < qlen) {
		v32i8_t const v = _loadu_v32i8(&query[qpos - 32]);
		tm_extend_count_base_core(&acc[1], v, 0);
	}

	/* tail */ {
		v32i8_t const v = _loadu_v32i8(&query[qpos - 32]);
		uint32_t const window = 0xffffffff<<(qlen & 0x1f);
		tm_extend_count_base_core(&acc[1], v, window);
	}

	/* derive count of 'A' from the others and length */
	acc[0] = qlen - (acc[1] + acc[2] + acc[3]);
	return;
}

static _force_inline
int32_t tm_extend_calc_complexity(tm_scan_t const *self, tm_idx_profile_t const *pf, uint8_t const *query, size_t qlen)
{
	_unused(self);

	/* initialize accumulators */
	size_t acc[4] = { 0 };
	tm_extend_count_base(acc, query, qlen);

	/* adjustment */
	size_t const tot = acc[0] + acc[1] + acc[2] + acc[3];
	double adj = 0.0;
	for(size_t i = 0; i < 4; i++) {
		if(acc[i] == 0) { continue; }

		double const c = (double)acc[i];
		double const f = (double)(acc[i]<<2) / (double)tot;

		adj -= c * log(f);
	}

	debug("(%zu, %zu, %zu, %zu), adj(%f)", acc[0], acc[1], acc[2], acc[3], adj);
	return((int32_t)(adj * pf->stat.rlambda));
}

static _force_inline
tm_score_t tm_extend_patch_score(tm_scan_t const *self, tm_idx_profile_t const *pf, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, tm_pos_t epos, dz_alignment_t const *aln)
{
	_unused(qlen);

	/* add bonus if anchored at the head */
	uint32_t const bonus = tm_extend_is_tail_maximal(self, sk, epos) ? pf->extend.bonus : 0;
	int32_t const raw    = bonus + aln->score;
	// fprintf(stderr, "dir(%u), (%u, %u), rend(%u), maximal(%u), bonus(%u)\n", epos.dir, epos.r, tm_aln_span(aln).r, rend, epos.r == rend, bonus);

	/* calc complexity */
	uint32_t const qspan = aln->query_length;
	uint32_t const qpos  = epos.q - qspan;
	int32_t const adj = tm_extend_calc_complexity(self, pf, &query[qpos], qspan);

	// fprintf(stderr, "score(%d, %d)\n", raw, raw + adj);
	return((tm_score_t){
		.raw     = raw,
		.patched = MAX2(0, raw + adj)
	});
}

static _force_inline
uint64_t tm_extend_single(tm_scan_t *self, tm_idx_profile_t const *pf, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, tm_chain_t const *q)
{
	tm_extend_res_t const r = tm_extend_core(self, pf, sk, query, qlen, q);
	if(r.aln == NULL) { return(0); }

	/* check score */
	tm_score_t const score = tm_extend_patch_score(self, pf, sk, query, qlen, r.epos, r.aln);
	int32_t const s = pf->extend.use_raw ? score.raw : score.patched;
	if(s <= pf->extend.min_score) { return(0); }

	/* save alignment; discard traceback object if the alignment is filtered out by an existing one */
	tm_aln_t const a = tm_extend_compose_aln(self, q, r.epos, r.aln, score);

	/* interval tree not update when the new one is discarded */
	return(tm_extend_record(self, &a));		/* 1 if succeeded */
}

static _force_inline
size_t tm_extend_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *query, size_t qlen, tm_chain_t const *chain, size_t ccnt)
{
	tm_idx_sketch_t const **sk  = (tm_idx_sketch_t const **)idx->sketch.arr;
	tm_idx_profile_t const **pf = (tm_idx_profile_t const **)idx->profile.arr;

	/* clear bin */
	tm_extend_clear(self);
/*
	for(size_t i = 0; i < qlen; i++) {
		fprintf(stderr, "%c", "A   C   G   T   "[query[i]]);
	}
	fprintf(stderr, "\n");
	for(size_t j = 0; j < DZ_REF_MAT_SIZE; j++) {
		for(size_t i = 0; i < DZ_QUERY_MAT_SIZE; i++) {
			fprintf(stderr, "ref(%zu), query(%zu), score(%d)\n", j, i, pf[0]->extend.dz->matrix[j * DZ_QUERY_MAT_SIZE + i]);
		}
	}
*/

	/* for each chain try extension */
	for(size_t i = 0; i < ccnt; i++) {
		tm_chain_t const *q = &chain[i];
		if(q->attr.sep.weight >= self->extend.max_weight) { break; }

		/* skip if already covered (and there is no possibility that this chain surpasses the previous ones) */
		if(tm_extend_is_covered(self, q)) { continue; }

		/* save current traceback stack for unwinding */
		dz_freeze_t const *fz = dz_arena_freeze(self->extend.trace);

		/* extend */
		tm_idx_sketch_t const *s  = sk[q->attr.sep.rid];
		tm_idx_profile_t const *p = pf[s->h.pid];
		debug("%r, rid(%u), pid(%u), rname(%s)", tm_chain_to_str, q, q->attr.sep.rid, tm_idx_ref_pid(s), tm_idx_ref_name_ptr(s));

		if(tm_extend_single(self, p, s, query, qlen, q) == 0) {
			/* failure; unwind traceback stack */
			dz_arena_restore(self->extend.trace, fz);
		} else {
			/* succeeded */
			if(tm_extend_is_complete(self)) { break; }
		}
	}

	/* anything else? */
	return(kv_cnt(self->extend.arr) - 1);
}




/* evaluate all; query sequence must be encoded in 2bit at [3:2] and shorter than 2Gbp */
static _force_inline
void tm_scan_convert_aln(tm_idx_sketch_t const *si, char const *qname, size_t qlen, baln_aln_t *dst, tm_aln_t const *src)
{
	char const *rname = (char const *)tm_idx_ref_name_ptr(si);
	size_t const rlen = tm_idx_ref_seq_len(si);

	tm_pair_t const span = src->span;
	tm_pos_t const spos = src->pos;
	tm_pos_t const pos = tm_canon_pos(spos, span);
	tm_aln_stat_t const stat = tm_aln_calc_stat(src, 0);

	*dst = (baln_aln_t){
		.name = { .r = rname,  .q = qname },
		.len  = { .r = rlen,   .q = qlen },
		.pos  = { .r = pos.r,  .q = pos.q },
		.span = { .r = span.r, .q = span.q },

		.dir  = pos.dir,
		.identity = stat.identity,
		.score = {
			.raw     = src->score.raw,
			.patched = src->score.patched
		},
		.path = {
			.ptr = src->aln->path,
			.len = src->aln->path_length
		}
	};
	return;
}

static _force_inline
size_t tm_scan_pack_aln(tm_idx_t const *idx, char const *qname, uint8_t const *query, size_t qlen, baln_aln_t *dst, tm_aln_t const *src)
{
	_unused(query);
	tm_idx_sketch_t const **si = (tm_idx_sketch_t const **)idx->sketch.arr;

	/* sort by qpos: traverse tree */
	tm_aln_t const all = {
		.pos  = { .q = 0 },
		.span = { .q = INT32_MAX },
		// .qmax = INT32_MAX
	};

	rbt_iter_t it;
	rbt_init_iter_isct_aln(&it);

	baln_aln_t *q = dst;
	while(1) {
		tm_aln_t const *p = rbt_fetch_isct_aln(&it, src, &all);
		if(p == NULL) { break; }

		debug("bid(%zu), %r", q - dst, tm_aln_to_str, p);
		tm_scan_convert_aln(si[p->attr.rid], qname, qlen, q++, p);
	}
	size_t const saved = q - dst;
	return(saved);
}

static _force_inline
baln_alnv_t *tm_scan_finalize(tm_scan_t *self, tm_idx_t const *idx, char const *sname, uint8_t const *seq, size_t slen)
{
	/* results */
	tm_aln_t const *src = kv_ptr(self->extend.arr);
	size_t const cnt = kv_cnt(self->extend.arr) - 1;

	/* allocate dst array */
	baln_alnv_t *dst = dz_arena_malloc(self->extend.trace,
		  sizeof(baln_alnv_t)		/* header */
		+ sizeof(baln_aln_t) * cnt	/* payload */
	);
	debug("src(%p), cnt(%zu), dst(%p)", src, cnt, dst);

	debugblock({
		rbt_print_aln(src, cnt, tm_aln_to_str);
		for(size_t i = 1; i < kv_cnt(self->extend.arr); i++) {
			tm_aln_t const *p = &src[i];
			debug("bid(%zu), %r", i, tm_aln_to_str, p);
		}
	});

	/* pack; cnt returned */
	dst->cnt = tm_scan_pack_aln(idx, sname, seq, slen, dst->arr, src);
	dst->bin = NULL;
	dst->free.body = NULL;
	dst->free.bin  = NULL;
	debug("cnt(%zu, %zu)", cnt, dst->cnt);
	return(dst);
}

// static _force_inline
baln_alnv_t *tm_scan_all(tm_scan_t *self, tm_idx_t const *idx, char const *sname, uint8_t const *seq, size_t slen)
{
	tm_scan_clear(self);

	/* seed and chain */
	if(tm_collect_all(self, idx, seq, slen) == 0) {
		return(0);
	}

	/* sort by estimated identity */
	tm_chain_t *cptr = kv_ptr(self->chain.arr);
	size_t const ccnt = kv_cnt(self->chain.arr);
	radix_sort_chain(cptr, ccnt);

	/* convert chain to alignment */
	if(tm_extend_all(self, idx, seq, slen, cptr, ccnt) == 0) {
		return(NULL);		/* alignment not found */
	}

	/* copy alignment to dst buffer */
	return(tm_scan_finalize(self, idx, sname, seq, slen));
}


/**
 * end of align.c
 */
