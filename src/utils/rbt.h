
/**
 * @file rbt.h
 * @brief red-black tree and interval tree
 * FIXME: tightly dependent on x86_64 architecture
 */

#ifndef _RBT_H_INCLUDED
#define _RBT_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"
#include "arch.h"


typedef struct rbt_header_s { uint32_t children[2]; } rbt_header_t;
#define RBT_SCALING_FACTOR			( 4 )
#define RBT_IMASK					( 0xfffffffc )
#define RBT_BMASK					( 0x100000001ULL )
#define _rbt_ptr(_type, _p, _i)		(_type *)(&(_p)[(_i) * sizeof(_type) / (RBT_SCALING_FACTOR * sizeof(uint32_t))])
#define _rbt_idx(_raw_idx)			( (_raw_idx) & RBT_IMASK )
#define _rbt_nop(_x, _y)			;

/* node iterator */
typedef struct rbt_iter_s {
	uint64_t dv;
	uint32_t *b;
	uint32_t ibuf[64];
} rbt_iter_t;

static _force_inline
void rbt_iter_init_static(struct rbt_iter_s *w)
{
	/*
	 * initialize direction bitvector
	 * root is right node of the god node; offsetted by 5bits; removed afterward (see below)
	 */
	w->dv = 0x01ULL;

	/*
	 * initialize index buffer
	 * w->ibuf[0] is for "parent of root" node, w->ibuf[1] for root node
	 */
	w->b = &w->ibuf[1];
	w->ibuf[0] = 0;
	return;
}
static _force_inline
uint64_t rbt_iter_remove_left_parents(struct rbt_iter_s *w)
{
	/* remove the tail contiguous Rs */
	uint64_t rcnt = _tzc_u64(~w->dv);
	debug("remove left parent, rcnt(%lu), dv(%lx -> %lx), b(%p)", rcnt, w->dv, w->dv>>rcnt, w->b - rcnt);

	w->dv >>= rcnt;
	w->b -= rcnt;
	return(w->dv == 0);			/* return true if reached "parent of root" */
}

enum rbt_shuffle_e { rbCl = 0, rbCr, rbBl, rbBr, rbAl, rbAr, rbAi, rbR, rbB };
#define _rbt_mask(_cl, _cr, _bl, _br, _al, _ar, _ai) \
	2*(_cl), 2*(_cl)+1, 2*(_cr), 2*(_cr)+1, 2*(_bl), 2*(_bl)+1, 2*(_br), 2*(_br)+1, \
	2*(_al), 2*(_al)+1, 2*(_ar), 2*(_ar)+1, 2*(_ai), 2*(_ai)+1, 0x80, 0x80

/*
 * We assume the nodes of tree is packed in a single vector which is addressable by non-negative index i on `arr'.
 * The nodes are located by uint32_t *ptr and index multiple of 4, which forces the the sizeof(node) be multiple of 16 (== 4 * sizeof(uint32_t)).
 * The first node of the array (arr[0]) is reserved for parent node of root, which cannot be used for storing any data.
 * The "parent of root" node is initialized by calling rbt_init_static function.
 *
 * _bucket_t:              container type, must contain rbt_header_t.
 * _hdr(p):                a macro which takes pointer to _bucket_t and returns pointer to rbt_header_t.
 * _cmp(p1, p2):           a macro which takes two pointers to _bucket_t and returns comparison result of p1->val < p2->val.
 * _update(parent, child): finalizer for parent -> child relation.
 */
#define RBT_INIT_IVT(_sfx, _bucket_t, _hdr, _cmp, _update) \
	typedef uint64_t (*rbt_callback_##_sfx##_t)(_bucket_t *node, _bucket_t *parent); \
	static _force_inline \
	uint64_t rbt_print_intl_##_sfx(FILE *fp, uint32_t *depth, uint8_t *check, uint32_t const *p, uint32_t raw_idx, uint32_t tab) { \
		static char const *spaces = "                                                                                                "; \
		if(_rbt_idx(raw_idx) == 0) { \
			if(tab > depth[0]) { depth[0] = tab; } \
			if(tab < depth[1]) { depth[1] = tab; } \
			return(0); \
		} \
		uint64_t err = check[raw_idx>>2] ? 1 : 0; check[raw_idx>>2]++; \
		uint32_t left = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx)))->children[0]; \
		uint32_t right = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx)))->children[1]; \
		if((raw_idx & 0x01) == 0 && ((left & 0x01) == 0 || (right & 0x01) == 0)) { \
			err |= 0x80; check[raw_idx>>2] |= 0x80; \
		} \
		err |= rbt_print_intl_##_sfx(fp, depth, check, p, right, tab + 1); \
		fprintf(fp, "%.*s(%u, %c)\n", (int)(2 * tab), spaces, _rbt_idx(raw_idx), (raw_idx & 0x01) ? 'B' : 'R'); \
		err |= rbt_print_intl_##_sfx(fp, depth, check, p, left, tab + 1); \
		return(err); \
	} \
	static _force_inline \
	void rbt_print_##_sfx(_bucket_t const *arr, size_t len) { \
		uint8_t *check = calloc(1, len + 1); \
		uint32_t depth[2] = { 0, UINT32_MAX }; \
		uint32_t const *p = (uint32_t const *)arr; \
		uint64_t err = rbt_print_intl_##_sfx(stderr, depth, check, p, _hdr(_rbt_ptr(_bucket_t, p, 0))->children[1], 0); \
		printf("depth(%u, %u)\n", depth[0], depth[1]); \
		if(err) { \
			fprintf(stderr, "broken\n"); \
			for(size_t i = 0; i < len; i++) { \
				if(check[i] < 2) { continue; } \
				fprintf(stderr, "(%lu, %u, (%u, %c), (%u, %c)), count(%u)\n", i * RBT_SCALING_FACTOR, _loadu_u32(&p[(i * RBT_SCALING_FACTOR) + 2]), \
					p[i * RBT_SCALING_FACTOR] & 0xfffc, (p[i * RBT_SCALING_FACTOR] & 0x01) ? 'B' : 'R', \
					p[(i * RBT_SCALING_FACTOR) + 1] & 0xfffc, (p[(i * RBT_SCALING_FACTOR) + 1] & 0x01) ? 'B' : 'R', \
					check[i]); \
			} \
		} \
		free(check); \
		return; \
	} \
	static _force_inline \
	void rbt_init_static_##_sfx(_bucket_t *arr) { \
		uint32_t *p = (uint32_t *)arr; \
		_storeu_u64(_hdr(_rbt_ptr(_bucket_t, p, 0)), RBT_BMASK);					/* leaves; both black */ \
		return; \
	} \
	/* find path to a node that is appropriate for the leaf to be inserted */ \
	static _force_inline \
	void rbt_find_leaf_intl_##_sfx(struct rbt_iter_s *w, uint32_t const *p, uint32_t raw_root_idx, _bucket_t const *leaf) { \
		debug("root(%u, %u)", _rbt_idx(raw_root_idx), raw_root_idx & 0x01); \
		/* add 5-bit offset */ \
		w->dv <<= 5; \
		/* dig tree from root; "parent of root" is at [0] and root is right node of the node */ \
		for(uint32_t raw_idx = raw_root_idx; ((int64_t)w->dv) > 0 && _rbt_idx(raw_idx) != 0;) { \
			_bucket_t const *node = _rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx));		/* extract pointer for the current node */ \
			uint64_t c = _loadu_u64(_hdr(node));						/* load (right, left) index pair */ \
			uint64_t shift = _cmp(node, leaf) ? 0x20 : 0;				/* node->val < leaf->val ? RIGHT : LEFT */ \
			debug("node at raw_idx(%u), left(%lu, %lu), right(%lu, %lu), go %s", \
				_rbt_idx(raw_idx), \
				_rbt_idx(c), c & 0x01, \
				_rbt_idx((c>>32)), (c>>32) & 0x01, shift ? "right" : "left"); \
			/* go down the tree by one */ \
			*w->b++ = raw_idx;			/* save current node index */ \
			w->dv = (w->dv<<1) + shift;	/* save direction */ \
			raw_idx = c>>shift;			/* extract next node index */ \
		} \
		w->dv >>= 5;					/* remove 5-bit offset (see above) */ \
		return; \
	} \
	/* \
	 * locate the first node that node->val >= anchor->val \
	 * (note: _cmp(p1, p2) must return p1->val < p2->val, not p1->val <= p2->val so that this function work correctly) \
	 */ \
	static _force_inline \
	_bucket_t const *rbt_find_right_##_sfx(_bucket_t const *arr, struct rbt_iter_s *w, _bucket_t const *anchor) { \
		uint32_t const *p = (uint32_t const *)arr; \
		/* clear iterator */ \
		rbt_iter_init_static(w); \
		/* find leaf location */ \
		uint32_t const raw_root_idx = _hdr(_rbt_ptr(_bucket_t, p, 0))->children[1]; \
		rbt_find_leaf_intl_##_sfx(w, p, raw_root_idx, anchor); \
		/* remove the tail contiguous Rs; returns true if no node is found */ \
		if(rbt_iter_remove_left_parents(w)) { return(NULL); } \
		/* here b[0]->val < anchor->val and b[-1]->val >= anchor->val */ \
		return(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]))); \
	} \
	static _force_inline \
	void rbt_insert_load_node_##_sfx(rbt_iter_t *w, uint32_t *p, uint32_t node_idx) { \
		/* save leaf */ \
		uint32_t raw_cidx = RBT_SCALING_FACTOR * node_idx;		/* leaf node index */ \
		_bucket_t *nptr = _rbt_ptr(_bucket_t, p, raw_cidx); \
		_storeu_u64(_hdr(nptr), RBT_BMASK);						/* mark two children of the leaf black (both empty) */ \
		_update(nptr, NULL);						/* invoke finalizer if needed */ \
		/* clear iterator */ \
		rbt_iter_init_static(w); \
		/* find parent; dv[0] for B->C direction */ \
		uint32_t const raw_root_idx = _hdr(_rbt_ptr(_bucket_t, p, 0))->children[1]; \
		rbt_find_leaf_intl_##_sfx(w, p, raw_root_idx, nptr); \
		/* parent of leaf found; save leaf index at the tail */ \
		*w->b = raw_cidx; \
		return; \
	} \
	static _force_inline \
	uint64_t rbt_insert_recolor_##_sfx(rbt_iter_t *w, uint32_t *p) { \
		uint32_t cidx = _rbt_idx(w->b[0]);			/* current node (leaf) index */ \
		uint32_t pidx = _rbt_idx(w->b[-1]);			/* parent node index */ \
		uint64_t ghdr = 0; \
		/* until (parent, uncle, grandparent) == (red, red, black) breaks */ \
		while(_likely(w->b > &w->ibuf[2])) { \
			/* update parent node */ \
			rbt_header_t *pphdr = _hdr(_rbt_ptr(_bucket_t, p, pidx));			/* load parent points */ \
			pphdr->children[w->dv & 0x01] = cidx;				/* mark red (since color bit is masked out by _rbt_idx); dv[0] for B->C direction */ \
			/* load grandparent */ \
			uint32_t gidx = _rbt_idx(w->b[-2]);						/* grandparent index */ \
			ghdr = _loadu_u64(_hdr(_rbt_ptr(_bucket_t, p, gidx)));				/* load parent and uncle (index, color) tuples */ \
			/* check if we need further fixup; break the loop if either parent or uncle is black; here b points grandparent and dv[0] indicates direction of C */ \
			if(_unlikely(ghdr & RBT_BMASK)) { break; } \
			/* need; mark parent and uncle black and continue */ \
			_storeu_u64(_hdr(_rbt_ptr(_bucket_t, p, gidx)), ghdr | RBT_BMASK);/* overwrite; parent and uncle are black */ \
			/* current and parent are fixed */ \
			_update(_rbt_ptr(_bucket_t, p, pidx), _rbt_ptr(_bucket_t, p, cidx));	/* parent -> current */ \
			_update(_rbt_ptr(_bucket_t, p, gidx), _rbt_ptr(_bucket_t, p, pidx));	/* grandparent -> parent */ \
			/* move up tree by two nodes */ \
			w->b -= 2; \
			w->dv >>= 2; \
			/* load next parent */ \
			cidx = gidx;				/* color is already masked out (see above) */ \
			pidx = _rbt_idx(w->b[-1]);	/* mask out color */ \
		} \
		return(ghdr); \
	} \
	static _force_inline \
	void rbt_insert_rebalance_##_sfx(rbt_iter_t *w, uint32_t *p, uint64_t ghdr) { \
		/* flip uncle color for test */ \
		uint64_t ghdr_flipped = ghdr ^ ((w->dv & 0x02ULL) ? 0x01ULL : (0x01ULL<<32)); \
		/* check if tree is high enough for rotate (the first condition) and the subtree is unbalanced (the second) */ \
		if(_unlikely((w->b > &w->ibuf[2]) & ((ghdr_flipped & RBT_BMASK) == 0))) { \
			/* \
			 * rotate last three nodes; A might be root                          \
			 * case LL (0):     case LR(1):    case RL(2):    case RR(3):        \
			 *     A              A            A              A                  \
			 *    /       B      /       C      \       C      \         B       \
			 *   B   ->  / \    B   ->  / \      B ->  / \      B   ->  / \      \
			 *  /       C   A    \     B   A    /     A   B      \     A   C     \
			 * C                  C            C                  C              \
			 */ \
			static uint8_t const shuffle_mask[64] __attribute__(( aligned(64) )) = { \
			/* new   (Cl, Cr, Bl, Br, Al, Ar, aidx, 0); payloads are not moved */ \
				_rbt_mask(rbCl, rbCr, rbBl, rbAi, rbBr, rbAr, rbAl),	/* LL; left-rotate */ \
				_rbt_mask(rbAl, rbAi, rbBl, rbCl, rbCr, rbAr, rbBr),	/* LR: right-rotate -> left-rotate */ \
				_rbt_mask(rbAi, rbAr, rbCr, rbBr, rbAl, rbCl, rbBl),	/* RL: left-rotate -> right-rotate */ \
				_rbt_mask(rbCl, rbCr, rbAi, rbBr, rbAl, rbBl, rbAr)		/* RR: right-rotate */ \
			}; \
			static uint8_t const gather_mask[16] = { 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15 }; \
			/* load node pointers (with masking colors out) */ \
			rbt_header_t *ap = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-2]))); \
			rbt_header_t *bp = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]))); \
			rbt_header_t *cp = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[0]))); \
			/* compute shuffle masks */ \
			uint64_t const y = (w->dv<<4) & 0x30, z = ~(3 * y) & 0x20; \
			uint64_t const color_mask = 0x10001ULL<<z; \
			__m128i const gmask = _mm_load_si128((__m128i const *)gather_mask); \
			__m128i const rmask = _mm_load_si128((__m128i const *)&shuffle_mask[y]); \
			/* gather indices */ \
			__m128i bc = _mm_insert_epi64(_mm_loadl_epi64((__m128i const *)cp), _loadu_u64(bp), 1); \
			__m128i an = _mm_insert_epi64(_mm_loadl_epi64((__m128i const *)ap), _rbt_idx(w->b[-2]), 1); \
			/* shuffle */ \
			bc = _mm_shuffle_epi8(bc, gmask); \
			an = _mm_shuffle_epi8(an, gmask); \
			__m128i l = _mm_shuffle_epi8(_mm_unpacklo_epi64(bc, an), rmask); \
			__m128i h = _mm_shuffle_epi8(_mm_unpackhi_epi64(bc, an), rmask); \
			l = _mm_andnot_si128(_mm_cvtsi64_si128(color_mask), l); \
			/* repack */ \
			bc = _mm_unpacklo_epi16(l, h); \
			an = _mm_unpackhi_epi16(l, h); \
			/* store */ \
			_storeu_u64(cp, _mm_cvtsi128_si64(bc)); \
			_storeu_u64(bp, _mm_extract_epi64(bc, 1)); \
			_storeu_u64(ap, _mm_cvtsi128_si64(an)); \
			/* adjust pointer */ \
			w->b -= 2; \
			w->dv >>= 2; \
			/* fixup links for parent node */ \
			w->b[0] = _mm_extract_epi32(an, 2) | 0x01;		/* mark black */ \
		} \
	} \
	static _force_inline \
	void rbt_insert_finalize_##_sfx(rbt_iter_t *w, uint32_t *p) { \
		/* fix grandparent -> parent link */ \
		rbt_header_t *pphdr = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]))); \
		pphdr->children[w->dv & 0x01] = w->b[0];			/* marked black; see above */ \
		/* move upward until hit root */ \
		while(w->b > &w->ibuf[1]) { \
			_update(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1])), _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[0]))); \
			w->b--; \
		} \
		/* recolor root black */ \
		rbt_header_t *pxhdr = _hdr(_rbt_ptr(_bucket_t, p, 0));			/* "parent of root" header pointer */ \
		pxhdr->children[1] |= 0x01;							/* dv[0] for B->C direction */ \
		debug("cidx(%u), pidx(%u), dv(%lx)", _rbt_idx(w->b[0]), _rbt_idx(w->b[-1]), w->dv); \
		debug("update %s child(%u) of last parent(%u), recolor root(%u)", (w->dv & 0x01) ? "right" : "left", _rbt_idx(w->b[0]), _rbt_idx(w->b[-1]), p[1]); \
	} \
	/*  insert functions: new node is located by node_idx, which is always an element of arr */ \
	static _force_inline \
	void rbt_insert_##_sfx(_bucket_t *arr, uint32_t node_idx) { \
		uint32_t *p = (uint32_t *)arr; \
		/* load node */ \
		rbt_iter_t w; \
		rbt_insert_load_node_##_sfx(&w, p, node_idx); \
		/* recolor */ \
		uint64_t ghdr = rbt_insert_recolor_##_sfx(&w, p); \
		/* rebalance */ \
		rbt_insert_rebalance_##_sfx(&w, p, ghdr); \
		/* save index and recolor root black */ \
		rbt_insert_finalize_##_sfx(&w, p); \
		return; \
	}

/* without _update */
#define RBT_INIT(_sfx, _bucket_t, _hdr, _cmp)			RBT_INIT_IVT(_sfx, _bucket_t, _hdr, _cmp, _rbt_nop)

/*
 * node searcher: rbt_find_head initializes iterator and rbt_find_next fetches next one by one
 *
 * rbt_iter_t it;
 * bucket_t b = { .val = v };
 * bucket_t *p = rbt_find_head(&it, arr, &v);
 * while(p != NULL) {
 *     // do something here
 *     rbt_find_next(&it, arr, &v);
 * }
 *
 * range iterator on normal tree:
 *   _cmp_head(p1, p2)			( (p1)->val >= (p2)->val )
 *   _cmp_tail(p1, p2)			( (p1)->val <  (p2)->val )
 *   _cmp_iter(p1, p2)			( 1 )
 *   where head_anchor and tail_anchor containing lb and ub (for [lb, ub) range) respectively.
 *
 * range iterator on interval tree:
 *   max_rval is defined maximum right boundary of all the children of the node:
 *     _update(parent, child)	{ if(child->max_rval > parent->max_rval) { parent->max_rval = child->max_rval; } }
 *
 *   contained range query (p1 contained in p2):
 *     _cmp_head(p1, p2)		( (p1)->val             >= (p2)->val )
 *     _cmp_tail(p1, p2)		( (p1)->val             <  (p2)->val + (p2)->len )
 *     _cmp_iter(p1, p2)		( (p1)->val + (p1)->len <  (p2)->val + (p2)->len )
 *     where head_anchor and tail_anchor are the same; containing lb and ub (for [lb, ub) range).
 *
 *   containing range query (p2 contained in p1):
 *     _cmp_head(p1, p2)		( (p1)->max_rval        >= (p2)->val + (p2)->len )
 *     _cmp_tail(p1, p2)		( (p1)->val             <= (p2)->val )
 *     _cmp_iter(p1, p2)		( (p1)->val + (p1)->len >= (p2)->val + (p2)->len )
 *
 *   intersection query (p1 overlaps with p2):
 *     _cmp_head(p1, p2)		( (p1)->max_rval        > (p2)->val )
 *     _cmp_tail(p1, p2)		( (p1)->val             < (p2)->val + (p2)->len )
 *     _cmp_iter(p1, p2)		( (p1)->val + (p1)->len > (p2)->val )
 */
#define RBT_INIT_ITER(_sfx, _bucket_t, _hdr, _cmp_head, _cmp_tail, _cmp_iter) \
	static _force_inline \
	uint64_t rbt_fetch_head_intl_##_sfx(struct rbt_iter_s *w, uint32_t const *p, uint32_t raw_root_idx, _bucket_t const *head_anchor) { \
		debug("root(%u, %u)", _rbt_idx(raw_root_idx), raw_root_idx & 0x01); \
		/* go left while the condition holds */ \
		for(uint32_t raw_idx = raw_root_idx; ((int64_t)w->dv) > 0 && _rbt_idx(raw_idx) != 0;) { \
			_bucket_t const *node = _rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx)); \
			if(!_cmp_head(node, head_anchor)) { debug("break"); break; } \
			/* go down the tree by one */ \
			debug("go left, idx(%u, %u)", _rbt_idx(raw_idx), raw_idx & 0x01); \
			*w->b++ = raw_idx;					/* save current node index */ \
			w->dv <<= 1;						/* move root flag */ \
			raw_idx = _hdr(node)->children[0];	/* extract left child index */ \
		} \
		debug("done, dv(%lx)", w->dv); \
		return((w->dv & 0x01ULL) == 0);			/* true if we were able to move left no more than once */ \
	} \
	static _force_inline \
	uint64_t rbt_fetch_next_intl_##_sfx(struct rbt_iter_s *w, uint32_t const *p, _bucket_t const *tail_anchor) { \
		/* find right; (w->dv & 0x02, w->b[-1]) points current node */ \
		while(1) { \
			/* try right children */ \
			uint32_t raw_idx = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1])))->children[1]; \
			w->dv |= 0x01ULL;					/* mark next move is right (tentatvely) */ \
			if(!rbt_fetch_head_intl_##_sfx(w, p, raw_idx, tail_anchor)) { \
				/* right child not found */ \
				debug("right child not found"); \
				if(rbt_iter_remove_left_parents(w)) { return(1); }		/* reached tail; not found */ \
			} \
			_bucket_t const *node = _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1])); \
			if(!_cmp_tail(node, tail_anchor)) { debug("node(%u) is out of tail", _rbt_idx(w->b[-1])); return(1); }	/* reached tail; not found */ \
			if(_cmp_iter(node, tail_anchor)) { debug("found, node(%u)", _rbt_idx(w->b[-1])); break; }		/* found */ \
			/* node does not satisfy the condition; poll the next */ \
			debug("poll next node(%u)", _rbt_idx(w->b[-1])); \
		} \
		return(0); \
	} \
	static _force_inline \
	uint64_t rbt_init_iter_##_sfx(struct rbt_iter_s *w, _bucket_t const *arr, _bucket_t const *head_anchor) { \
		uint32_t const *p = (uint32_t const *)arr; \
		/* initialize working buffer */ \
		rbt_iter_init_static(w); \
		/* left children */ \
		uint32_t const raw_root_idx = _hdr(_rbt_ptr(_bucket_t, p, 0))->children[1]; \
		return(rbt_fetch_head_intl_##_sfx(w, p, raw_root_idx, head_anchor)); \
	} \
	static _force_inline \
	_bucket_t const *rbt_fetch_head_##_sfx(struct rbt_iter_s *w, _bucket_t const *arr, _bucket_t const *tail_anchor) { \
		uint32_t const *p = (uint32_t const *)arr; \
		if(w->b <= &w->ibuf[1]) { return(NULL); }		/* is initial state; not found */ \
		/* find first node that satisfy the condition */ \
		_bucket_t const *node = _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1])); \
		if(!_cmp_tail(node, tail_anchor)) { debug("head node(%u) is out of tail", _rbt_idx(w->b[-1])); return(NULL); }	/* first node is out of tail; not found */ \
		if(!_cmp_iter(node, tail_anchor)) { \
			/* first node is inappropriate; find next */ \
			debug("not found; go next"); \
			if(rbt_fetch_next_intl_##_sfx(w, p, tail_anchor)) { return(NULL); } \
		} \
		/* found */ \
		debug("found, node(%u)", _rbt_idx(w->b[-1])); \
		return(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]))); \
	} \
	static _force_inline \
	_bucket_t const *rbt_fetch_next_##_sfx(struct rbt_iter_s *w, _bucket_t const *arr, _bucket_t const *tail_anchor) { \
		uint32_t const *p = (uint32_t const *)arr; \
		if(rbt_fetch_next_intl_##_sfx(w, p, tail_anchor)) { return(NULL); } \
		return(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]))); \
	}

#if 0
/* create node */
unittest()
{
	kvec_t(struct rbt_node_s) v = { 0 };
	kv_pushp(struct rbt_node_s, v);

	rbt_init_static(v.a);
	for(size_t i = 0; i < 16384<<12; i++) {
		debug("i(%lu), v(%p, %lu)", i, v.a, v.n);
		rbt_print(v.a, v.n, 2);
		rbt_foreach(v.a, 2, 64, 128, {
			fprintf(stderr, "%u\n", n->payload);
		});

		struct rbt_iter_s it;
		for(_bucket_t *n = rbt_find_node(v.a, 2, &it, 64);
			n != v.a && n->payload < 128;
			n = rbt_find_next(v.a, 2, &it)) {
			fprintf(stderr, "n(%p), idx(%lu), payload(%u)\n", n, n - v.a, n->payload);
		}

		struct rbt_node_s *p = kv_pushp(struct rbt_node_s, v);
		p->payload = rand() & 255;

		debug("(%p, %u)", p, p->payload);
		rbt_insert(v.a, 2, v.n - 1);
	}
	rbt_print(v.a, v.n, 2);
	free(v.a);
}
#endif

#endif
/**
 * end of rbt.h
 */
