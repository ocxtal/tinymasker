
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


typedef struct rbt_header_s {
	uint32_t children[2];
} rbt_header_t;

#define RBT_SCALING_FACTOR			( 4 )
#define RBT_IMASK					( 0xfffffffc )
#define RBT_BMASK					( 0x100000001ULL )
#define _rbt_ptr(_type, _p, _i)		(_type *)(&(_p)[(_i) * sizeof(_type) / (RBT_SCALING_FACTOR * sizeof(uint32_t))])
#define _rbt_idx(_raw_idx)			( (_raw_idx) & RBT_IMASK )
#define _rbt_nop(_x, _y, _z)		( 0 )

#define _rbt_is_black(_x)			( ((_x) & 0x01) == 0x01 )
#define _rbt_is_red(_x)				( ((_x) & 0x01) == 0x00 )
#define _rbt_black(_x)				( (_x) | 0x01 )
#define _rbt_red(_x)				( (_x) & ~0x01 )

/*
 * node iterator
 * invariant condition:
 *   ibuf[0] contains index of "parent of the root" node, this node always has index zero
 *   ibuf[1] contains index of the root node
 *   b points at ibuf of the current node; can be undefined in the incomplete state
 *   dv[0] is direction to the current (leaf) node from its parent
 */
typedef struct rbt_iter_s {
	uint64_t dv;			/* direction vector */
	uint32_t *b;			/* pointer to ibuf */
	uint32_t ibuf[64];		/* tree traversal array */
} rbt_iter_t;

static _force_inline
void rbt_iter_init_static(rbt_iter_t *w)
{
	/*
	 * initialize direction bitvector
	 * the root is right node of the god node
	 */
	w->dv = 0x01ULL;

	/*
	 * initialize index buffer
	 * w->ibuf[0] is for "parent of the root" node, w->ibuf[1] for the root node
	 */
	w->b = &w->ibuf[1];			/* *w->b uninitialized for now (incomplete state) */
	w->ibuf[0] = 0;
	return;
}

static _force_inline
uint64_t rbt_iter_ascend_left(rbt_iter_t *w)
{
	/* return true if reached "parent of root" */
	if((w->dv & (w->dv + 1)) == 0) { return(0); }

	/* remove the tail contiguous Rs */
	ZCNT_RESULT uint64_t const rcnt = _tzc_u64(~w->dv);
	debug("remove left parent, rcnt(%lu), dv(%lx -> %lx), b(%p)", rcnt, w->dv, w->dv>>rcnt, w->b - rcnt);

	w->dv >>= rcnt;
	w->b   -= rcnt;
	return(1);
}


/* rebalance stuffs */
typedef struct {
	uint32_t aidx;
	uint32_t cis;
} rbt_rebalance_t;

enum rbt_shuffle_e { rbCl = 0, rbCr, rbBl, rbBr, rbAl, rbAr, rbAi, rbR, rbB };
#define _rbt_mask(_cl, _cr, _bl, _br, _al, _ar, _ai) \
	2*(_cl), 2*(_cl)+1, 2*(_cr), 2*(_cr)+1, 2*(_bl), 2*(_bl)+1, 2*(_br), 2*(_br)+1, \
	2*(_al), 2*(_al)+1, 2*(_ar), 2*(_ar)+1, 2*(_ai), 2*(_ai)+1, 0x80, 0x80


/*
 * We assume the nodes of tree is packed in a single vector which is addressable by non-negative index i on `arr'.
 * The nodes are located by uint32_t *ptr and index multiple of 4, which forces the the sizeof(node_t) be multiple of 16 (== 4 * sizeof(uint32_t)).
 * The first node of the array (arr[0]) is reserved for parent node of root, which cannot be used for storing any data.
 * It must keep the minimum value so that the binary tree condition holds anywhere in the tree, otherwise either insert or iterator will be broken.
 * The "parent of root" node is initialized by calling rbt_init_static function.
 *
 * _bucket_t:              container type, must contain rbt_header_t.
 * _hdr(p):                a macro which takes pointer to _bucket_t and returns pointer to rbt_header_t.
 * _cmp(p1, p2):           a macro which takes two pointers to _bucket_t and returns comparison result of p1->val < p2->val.
 * _update(parent, child): finalizer for parent -> child relation. we expect the function returns non-zero when further update needed for parents.
 */
#define RBT_INIT_IVT(_sfx, _bucket_t, _hdr, _cmp, _update) \
	typedef uint64_t (*rbt_callback_##_sfx##_t)(_bucket_t *node, _bucket_t *parent); \
	static _force_inline \
	uint64_t rbt_print_intl_##_sfx(FILE *fp, char *(*formatter)(_bucket_t const *p), uint32_t *depth, uint8_t *check, uint32_t const *p, uint32_t raw_idx, uint32_t tab) { \
		static char const *spaces = "                                                                                                "; \
		if(_rbt_idx(raw_idx) == 0) { \
			if(tab > depth[0]) { depth[0] = tab; } \
			if(tab < depth[1]) { depth[1] = tab; } \
			return(0); \
		} \
		uint64_t err = check[raw_idx>>2] ? 1 : 0; check[raw_idx>>2]++; \
		uint32_t const left = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx)))->children[0]; \
		uint32_t const right = _hdr(_rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx)))->children[1]; \
		/* sanity check */ \
		if(_rbt_is_red(raw_idx) && (_rbt_is_red(left) || _rbt_is_red(right))) { \
			err |= 0x80; check[raw_idx>>2] |= 0x80; \
		} \
		/* print right child */ \
		err |= rbt_print_intl_##_sfx(fp, formatter, depth, check, p, right, tab + 1); \
		/* print self */ { \
			char *s = formatter != NULL ? formatter(_rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx))): NULL; \
			fprintf(fp, "%.*s(%u, %c)%s%s\n", \
				(int)(2 * tab), spaces, \
				_rbt_idx(raw_idx) / RBT_SCALING_FACTOR, (raw_idx & 0x01) ? 'B' : 'R', \
				s == NULL ? "" : ": ", s == NULL ? "" : s \
			); \
			free(s); \
		} \
		/* print left child */ \
		err |= rbt_print_intl_##_sfx(fp, formatter, depth, check, p, left, tab + 1); \
		return(err); \
	} \
	static _force_inline \
	void rbt_print_##_sfx(_bucket_t const *arr, size_t len, char *(*formatter)(_bucket_t const *p)) { \
		uint8_t *check = calloc(1, len + 1); \
		uint32_t depth[2] = { 0, UINT32_MAX }; \
		uint32_t const *p = (uint32_t const *)arr; \
		uint64_t const err = rbt_print_intl_##_sfx(stderr, formatter, depth, check, p, _hdr(_rbt_ptr(_bucket_t, p, 0))->children[1], 0); \
		printf("depth(%u, %u)\n", depth[0], depth[1]); \
		if(err) { \
			fprintf(stderr, "broken\n"); \
			for(size_t i = 0; i < len; i++) { \
				if(check[i] < 2) { continue; } \
				fprintf(stderr, "(%lu, %u, (%u, %c), (%u, %c)), count(%u)\n", i * RBT_SCALING_FACTOR, _loadu_u32(&p[(i * RBT_SCALING_FACTOR) + 2]), \
					p[i * RBT_SCALING_FACTOR] & 0xfffc, (p[i * RBT_SCALING_FACTOR] & 0x01) ? 'B' : 'R', \
					p[(i * RBT_SCALING_FACTOR) + 1] & 0xfffc, (p[(i * RBT_SCALING_FACTOR) + 1] & 0x01) ? 'B' : 'R', \
					check[i] \
				); \
			} \
		} \
		free(check); \
		return; \
	} \
	static _force_inline \
	void rbt_init_static_##_sfx(_bucket_t *arr) { \
		uint32_t *p = (uint32_t *)arr; \
		_storeu_u64(_hdr(_rbt_ptr(_bucket_t, p, 0)), RBT_BMASK);			/* leaves; both black */ \
		return; \
	} \
	/* \
	 * find path to a node that is appropriate for the leaf to be inserted \
	 */ \
	static _force_inline \
	void rbt_find_leaf_intl_##_sfx(rbt_iter_t *w, uint32_t const *p, _bucket_t const *leaf) { \
		uint32_t const raw_root_idx = w->b[-1]; \
		_bucket_t const *root = _rbt_ptr(_bucket_t, p, raw_root_idx); \
		debug("root(%u, %u)", _rbt_idx(raw_root_idx), raw_root_idx & 0x01); \
		/* add 5-bit offset */ \
		w->dv  |= 0x01ULL; \
		w->dv <<= 5; \
		/* dig tree from root; "parent of root" is at [0] and root is right node of the node */ \
		for(uint32_t raw_idx = _hdr(root)->children[1]; ((int64_t)w->dv) > 0 && _rbt_idx(raw_idx) != 0;) { \
			_bucket_t const *node = _rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx));		/* extract pointer for the current node */ \
			uint64_t const c = _loadu_u64(_hdr(node));						/* load (right, left) index pair */ \
			uint64_t const shift = _cmp(node, leaf) ? 0x20 : 0;				/* node->val < leaf->val ? RIGHT : LEFT */ \
			debug("node at raw_idx(%u), left(%lu, %lu), right(%lu, %lu), go %s", \
				_rbt_idx(raw_idx), \
				_rbt_idx(c), c & 0x01, \
				_rbt_idx((c>>32)), (c>>32) & 0x01, shift ? "right" : "left"); \
			/* go down the tree by one */ \
			*w->b++ = raw_idx;			/* save current node index; still incomplete state */ \
			w->dv   = (w->dv<<1) + shift;/* save direction */ \
			raw_idx = c>>shift;			/* extract next node index */ \
		} \
		w->dv >>= 5;					/* remove 5-bit offset (see above) */ \
		return; \
	} \
	/* \
	 * callback for patching metadata; for interval tree \
	 */ \
	static _force_inline \
	uint64_t rbt_insert_update_core_##_sfx(uint32_t *p, _bucket_t *pptr) { \
		uint32_t const lidx = _rbt_idx(_hdr(pptr)->children[0]); \
		uint32_t const ridx = _rbt_idx(_hdr(pptr)->children[1]); \
		_bucket_t *lptr = lidx == 0 ? NULL : _rbt_ptr(_bucket_t, p, lidx); \
		_bucket_t *rptr = ridx == 0 ? NULL : _rbt_ptr(_bucket_t, p, ridx); \
		_unused(lptr);		/* guard for nop */ \
		_unused(rptr); \
		return(_update(pptr, lptr, rptr));			/* might be nop; then pointers are unused */ \
	} \
	/* \
	 * locate the first node that node->val >= anchor->val \
	 * (note: _cmp(p1, p2) must return p1->val < p2->val, not p1->val <= p2->val so that this function work correctly) \
	 */ \
	static _force_inline \
	_bucket_t const *rbt_find_right_##_sfx(_bucket_t const *arr, rbt_iter_t *w, _bucket_t const *anchor) { \
		uint32_t const *p = (uint32_t const *)arr; \
		/* clear iterator */ \
		rbt_iter_init_static(w); \
		/* find leaf location */ \
		rbt_find_leaf_intl_##_sfx(w, p, anchor); \
		/* remove the tail contiguous Rs; returns true if no node is found */ \
		if(rbt_iter_ascend_left(w)) { \
			/* found */ \
			return(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]))); \
		} \
		/* here b[0]->val < anchor->val and b[-1]->val >= anchor->val */ \
		return(NULL); \
	} \
	/* \
	 * initialize header of the node to be inserted, then determine a location to be inserted \
	 */ \
	static _force_inline \
	void rbt_insert_load_node_##_sfx(rbt_iter_t *w, uint32_t *p, uint32_t node_idx) { \
		/* save leaf */ \
		uint32_t const raw_cidx = RBT_SCALING_FACTOR * node_idx;		/* leaf node index */ \
		_bucket_t *nptr = _rbt_ptr(_bucket_t, p, raw_cidx); \
		_storeu_u64(_hdr(nptr), RBT_BMASK);			/* mark two children of the leaf black (both empty) */ \
		rbt_insert_update_core_##_sfx(p, nptr);		/* invoke finalizer if needed; children pointers are always NULL (optimized out) */ \
		/* clear iterator */ \
		rbt_iter_init_static(w); \
		/* find parent; dv[0] for B->C direction */ \
		rbt_find_leaf_intl_##_sfx(w, p, nptr); \
		/* parent of leaf found; save leaf index at the tail */ \
		*w->b = raw_cidx;							/* inserted node is red in the first place; complete state afterward */ \
		debug("b(%p, %zu), ibuf(%p)", w->b, w->b - w->ibuf, w->ibuf); \
		return; \
	} \
	/* \
	 * ascend the tree while the red-black condition is not met \
	 */ \
	static _force_inline \
	uint64_t rbt_insert_recolor_##_sfx(rbt_iter_t *w, uint32_t *p) { \
		uint32_t cidx = _rbt_idx(w->b[0]);			/* current node (leaf) index */ \
		uint32_t pidx = _rbt_idx(w->b[-1]);			/* parent node index */ \
		uint64_t ghdr = 0; \
		/* until (parent, uncle, grandparent) == (red, red, black) breaks */ \
		while(_likely(w->b > &w->ibuf[2])) { \
			/* update parent node */ \
			rbt_header_t *pphdr = _hdr(_rbt_ptr(_bucket_t, p, pidx));				/* load parent points */ \
			pphdr->children[w->dv & 0x01] = cidx;									/* mark red (since color bit is masked out by _rbt_idx); dv[0] for B->C direction */ \
			/* load grandparent */ \
			uint32_t gidx = _rbt_idx(w->b[-2]);										/* grandparent index */ \
			ghdr = _loadu_u64(_hdr(_rbt_ptr(_bucket_t, p, gidx)));					/* load parent and uncle (index, color) tuples */ \
			/* check if we need further fixup; break the loop if either parent or uncle is black; here b points grandparent and dv[0] indicates direction of C */ \
			if(_unlikely(ghdr & RBT_BMASK)) { break; } \
			/* need; mark parent and uncle black and continue */ \
			_storeu_u64(_hdr(_rbt_ptr(_bucket_t, p, gidx)), ghdr | RBT_BMASK);/* overwrite; parent and uncle are black */ \
			/* current and parent are fixed */ \
			rbt_insert_update_core_##_sfx(p, _rbt_ptr(_bucket_t, p, pidx));	/* parent -> current */ \
			rbt_insert_update_core_##_sfx(p, _rbt_ptr(_bucket_t, p, gidx));	/* grandparent -> parent */ \
			/* move up tree by two nodes */ \
			w->b   -= 2;				/* the iterator is kept complete */ \
			w->dv >>= 2; \
			/* load next parent */ \
			cidx = gidx;				/* color is already masked out (see above) */ \
			pidx = _rbt_idx(w->b[-1]);	/* mask out color of the parent node */ \
		} \
		return(ghdr); \
	} \
	/* \
	 * subtree rotation; load three nodes of interest onto xmm registers and flip them with pshufbs. \
	 * flipping patterns are encoded in shuffle_mask and gather_mask tables. \
	 */ \
	static _force_inline \
	rbt_rebalance_t rbt_insert_rebalance_core_##_sfx(rbt_iter_t *w, uint32_t aidx, rbt_header_t *ap, rbt_header_t *bp, rbt_header_t *cp) { \
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
		static uint8_t const gather_mask[16] __attribute__(( aligned(16) )) = { \
			0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15 \
		}; \
		/* compute shuffle masks */ \
		uint64_t const pattern    = (w->dv<<4) & 0x30; \
		uint64_t const is_cis     = ~(3 * pattern) & 0x20;		/* ~(pattern[4] ^ pattern[5]); true when LL or RR */ \
		uint64_t const color_mask = 0x10001ULL<<is_cis;			/* shift mask when LL or RR to mark B children, otherwise mark C children */ \
		__m128i const gmask = _mm_load_si128((__m128i const *)gather_mask); \
		__m128i const rmask = _mm_load_si128((__m128i const *)&shuffle_mask[pattern]); \
		/* gather indices */ \
		__m128i bc = _mm_insert_epi64(_mm_loadl_epi64((__m128i const *)cp), _loadu_u64(bp), 1); \
		__m128i an = _mm_insert_epi64(_mm_loadl_epi64((__m128i const *)ap), aidx, 1); \
		/* shuffle */ \
		bc = _mm_shuffle_epi8(bc, gmask); \
		an = _mm_shuffle_epi8(an, gmask); \
		__m128i l = _mm_shuffle_epi8(_mm_unpacklo_epi64(bc, an), rmask); \
		__m128i h = _mm_shuffle_epi8(_mm_unpackhi_epi64(bc, an), rmask); \
		l = _mm_andnot_si128(_mm_cvtsi64_si128(color_mask), l); \
		/* repack */ \
		bc = _mm_unpacklo_epi16(l, h);		/* chdr in lower 64bit, bhdr in upper 64bit */ \
		an = _mm_unpackhi_epi16(l, h);		/* ahdr in lower 64bit, aidx in upper 64bit */ \
		/* store */ \
		_storeu_u64(cp, _mm_cvtsi128_si64(bc)); \
		_storeu_u64(bp, _mm_extract_epi64(bc, 1)); \
		_storeu_u64(ap, _mm_cvtsi128_si64(an)); \
		/* return grandparent node index; the grandparent is always black */ \
		return((rbt_rebalance_t){ \
			.aidx = _rbt_black(_mm_extract_epi32(an, 2)), \
			.cis  = is_cis \
		}); \
	} \
	/* \
	 * save link: grandparent -> parent \
	 */ \
	static _force_inline \
	void rbt_insert_patch_link_##_sfx(rbt_iter_t *w, uint32_t *p) { \
		/* fixup parent->node link; parent of grandparent may be root */ \
		_bucket_t *x = _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1])); \
		uint32_t *xbin = &_hdr(x)->children[w->dv & 0x01]; \
		*xbin = *w->b;		/* update parent of grandparent; mark black */ \
		return; \
	} \
	/* \
	 * rebalance nodes when recolor is done; the core algorithm is above; \
	 * we determine if rebalance is necessary or not. dispatch appropriate function if needed. \
	 */ \
	static _force_inline \
	uint64_t rbt_insert_rebalance_##_sfx(rbt_iter_t *w, uint32_t *p, uint64_t ghdr) { \
		/* flip uncle color for test */ \
		uint64_t const ghdr_flipped = ghdr ^ ((w->dv & 0x02ULL) ? 0x01ULL : (0x01ULL<<32)); \
		/* compute testing flag */ \
		uint64_t const ghdr_test = ( \
			   (w->b > &w->ibuf[2])					/* if the tree is high enough for rotate */ \
			& ((ghdr_flipped & RBT_BMASK) == 0)		/* if the subtree is unbalanced */ \
		); \
		if(_likely(ghdr_test == 0)) { \
			/* we don't need to rebalance the tree */ \
			debug("skip rebalance"); \
			rbt_insert_patch_link_##_sfx(w, p); \
			return(1); \
		} \
		/* continue to tree rebalancing; first compose pointers */ \
		debug("rebalance, (%u, %u, %u)", _rbt_idx(w->b[-2]), _rbt_idx(w->b[-1]), _rbt_idx(w->b[0])); \
		_bucket_t *a = _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-2]));		/* grandparent */ \
		_bucket_t *b = _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]));		/* parent */ \
		_bucket_t *c = _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[0]));		/* inserted node */ \
		/* rebalance; we get new aidx (grandparent index; black) */ \
		rbt_rebalance_t const r = rbt_insert_rebalance_core_##_sfx(w, _rbt_idx(w->b[-2]), _hdr(a), _hdr(b), _hdr(c)); \
		/* remove rotated subtree from iterator */ \
		w->b   -= 2; \
		w->dv >>= 2; \
		/* save grandparent node index */ \
		*w->b = r.aidx;						/* update tree iterator; kept complete */ \
		rbt_insert_patch_link_##_sfx(w, p);	/* patch link before updating metadata */ \
		/* update metadata for rebalanced nodes */ \
		rbt_insert_update_core_##_sfx(p, a); \
		uint64_t const cont = rbt_insert_update_core_##_sfx(p, b); \
		if(r.cis == 0) { return(cont); }	/* we don't need to update C node when LL or RR */ \
		return(rbt_insert_update_core_##_sfx(p, c)); \
	} \
	/* \
	 * patch metadata after inserting a node. mainly for interval tree. not intended to be exposed as API. \
	 */ \
	static _force_inline \
	void rbt_insert_update_##_sfx(rbt_iter_t *w, uint32_t *p) { \
		/* move upward until hit root */ \
		while(--w->b > &w->ibuf[0]) { \
			_bucket_t *x = _rbt_ptr(_bucket_t, p, _rbt_idx(w->b[0])); \
			if(!rbt_insert_update_core_##_sfx(p, x)) { break; } \
		} \
		w->b++; \
		return; \
	} \
	/* \
	 * patch the first (god) node. \
	 */ \
	static _force_inline \
	void rbt_insert_finalize_##_sfx(rbt_iter_t *w, uint32_t *p) { \
		_unused(w); \
		debug("fixup link, b(%p, %zu), ibuf(%p), dv(%lx), %u -> %u", w->b, w->b - w->ibuf, w->ibuf, w->dv, w->b[-1], w->b[0]); \
		/* recolor root black */ \
		rbt_header_t *rhdr = _hdr(_rbt_ptr(_bucket_t, p, 0));	/* "parent of root" header pointer */ \
		rhdr->children[1] = _rbt_black(rhdr->children[1]);		/* dv[0] for B->C direction */ \
	} \
	/* \
	 * insert functions: new node is located by node_idx, which is always an element of arr \
	 */ \
	static _force_inline \
	void rbt_insert_##_sfx(_bucket_t *arr, uint32_t node_idx) { \
		uint32_t *p = (uint32_t *)arr; \
		/* load node */ \
		rbt_iter_t w; \
		rbt_insert_load_node_##_sfx(&w, p, node_idx);			/* iterator is complete state afterward */ \
		/* recolor */ \
		uint64_t const ghdr = rbt_insert_recolor_##_sfx(&w, p); \
		/* rebalance */ \
		if(rbt_insert_rebalance_##_sfx(&w, p, ghdr)) { \
			rbt_insert_update_##_sfx(&w, p); \
		} \
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
 * three comparison macros, _cmp_head, _cmp_tail, and _cmp_iter, are required to instanciate iterator,
 * where the three are supposed to determine the following conditions:
 *   _cmp_head(node, range):
 *      supposed to return false when the node values (or ranges) of its subtree are
 *      too small to meet the iteration condition.
 *   _cmp_tail(node, range):
 *      supposed to return true when the node value (or range) is
 *      not too large that it meets the condition for now.
 *   _cmp_iter(node, range):
 *      supposed to return true when the node value (or range) meets the condition,
 *      applied for nodes passed the _cmp_tail test to filter the result.
 *
 * combinations of 
 *
 *
 * range iterator on normal tree:
 *   _cmp_head(p1, p2)			( (p1)->val >= (p2)->val )
 *   _cmp_tail(p1, p2)			( (p1)->val <  (p2)->val )
 *   _cmp_iter(p1, p2)			( 1 )
 *   where anchor containing lb and ub (for [lb, ub) range) respectively.
 *
 * key search query:
 *   _cmp_head(p1, p2)			( (p1)->val >= (p2)->val )
 *   _cmp_tail(p1, p2)			( (p1)->val <= (p2)->val )
 *   _cmp_iter(p1, p2)			( 1 )
 *   where anchor having the key searched.
 *
 * range iterator on interval tree:
 *   max_rval is defined maximum right boundary of all the children of the node:
 *     _update(parent, child)	{ if(child->max_rval > parent->max_rval) { parent->max_rval = child->max_rval; } }
 *
 *   contained range query (a set of nodes which are contained in the query range):
 *     _cmp_head(p1, p2)		( (p1)->val             >= (p2)->val )
 *     _cmp_tail(p1, p2)		( (p1)->val             <  (p2)->val + (p2)->len )
 *     _cmp_iter(p1, p2)		( (p1)->val + (p1)->len <  (p2)->val + (p2)->len )
 *     where anchor is containing lb and ub (for [lb, ub) range).
 *
 *   containing range query (a set of nodes each of which contains the query range):
 *     _cmp_head(p1, p2)		( (p1)->max_rval        >= (p2)->val + (p2)->len )
 *     _cmp_tail(p1, p2)		( (p1)->val             <= (p2)->val )
 *     _cmp_iter(p1, p2)		( (p1)->val + (p1)->len >= (p2)->val + (p2)->len )
 *
 *   intersection query (a set of nodes each of which intersects with the query range):
 *     _cmp_head(p1, p2)		( (p1)->max_rval        > (p2)->val )
 *     _cmp_tail(p1, p2)		( (p1)->val             < (p2)->val + (p2)->len )
 *     _cmp_iter(p1, p2)		( (p1)->val + (p1)->len > (p2)->val )
 *
 *
 * the iterator is kept complete state after rbt_fetch_head
 */
#define RBT_INIT_ITER_PATCH(_sfx, _bucket_t, _hdr, _cmp_head, _cmp_tail, _cmp_iter, _update) \
	/* \
	 * locate the leftmost node where _cmp_head (leftside condition for a subtree) becomes true. \
	 * return true if at least one node is available in the current subtree. \
	 */ \
	static _force_inline \
	void rbt_iter_find_leaf_intl_##_sfx(rbt_iter_t *w, uint32_t const *p, _bucket_t const *anchor) { \
		uint32_t const raw_root_idx = w->b[-1]; \
		_bucket_t const *root = _rbt_ptr(_bucket_t, p, _rbt_idx(raw_root_idx)); \
		debug("root(%u) -> right(%u)", raw_root_idx, _hdr(root)->children[1]); \
		/* add 5-bit offset */ \
		w->dv  |= 0x01ULL; \
		w->dv <<= 5; \
		/* dig tree from root; "parent of root" is at [0] and root is right node of the node */ \
		for(uint32_t raw_idx = _hdr(root)->children[1]; ((int64_t)w->dv) > 0 && _rbt_idx(raw_idx) != 0;) { \
			_bucket_t const *node = _rbt_ptr(_bucket_t, p, _rbt_idx(raw_idx));		/* extract pointer for the current node */ \
			uint64_t const c = _loadu_u64(_hdr(node));						/* load (right, left) index pair */ \
			uint64_t const shift = _cmp_head(node, anchor) ? 0 : 0x20;		/* node->val < anchor->val ? LEFT : RIGHT */ \
			debug("node at raw_idx(%u), left(%lu, %lu), right(%lu, %lu), go %s", \
				_rbt_idx(raw_idx), \
				_rbt_idx(c), c & 0x01, \
				_rbt_idx((c>>32)), (c>>32) & 0x01, shift ? "right" : "left"); \
			/* go down the tree by one */ \
			*w->b++ = raw_idx;			/* save current node index; still incomplete state */ \
			w->dv   = (w->dv<<1) + shift;/* save direction */ \
			raw_idx = c>>shift;			/* extract next node index */ \
		} \
		w->dv >>= 5;					/* remove 5-bit offset (see above) */ \
		return; \
	} \
	static _force_inline \
	uint64_t rbt_fetch_intl_##_sfx(rbt_iter_t *w, uint32_t const *p, _bucket_t const *anchor) { \
		while(1) { \
			rbt_iter_find_leaf_intl_##_sfx(w, p, anchor); \
			if(!rbt_iter_ascend_left(w)) { \
				debug("no children found, term."); \
				return(0); \
			} \
			uint32_t const raw_idx = _rbt_idx(w->b[-1]); \
			_bucket_t const *node  = _rbt_ptr(_bucket_t, p, raw_idx); \
			if(!_cmp_tail(node, anchor)) { \
				debug("node(%u) out of tail.", raw_idx); \
				return(0); \
			} \
			if(_cmp_iter(node, anchor)) { \
				debug("found node(%u)", raw_idx); \
				return(1); \
			} \
		} \
	} \
	/* \
	 * locate the first node where _cmp_head becomes true. \
	 */ \
	static _force_inline \
	void rbt_init_iter_##_sfx(rbt_iter_t *w) { \
		rbt_iter_init_static(w);	/* initialize working buffer */ \
		return; \
	} \
	static _force_inline \
	_bucket_t const *rbt_fetch_##_sfx(rbt_iter_t *w, _bucket_t const *arr, _bucket_t const *anchor) { \
		uint32_t const *p = (uint32_t const *)arr; \
		if(rbt_fetch_intl_##_sfx(w, p, anchor)) { \
			/* found */ \
			return(_rbt_ptr(_bucket_t, p, _rbt_idx(w->b[-1]))); \
		} \
		return(NULL); \
	} \
	/* \
	 * callback for patching metadata; for interval tree \
	 */ \
	static _force_inline \
	uint64_t rbt_patch_core_##_sfx(uint32_t *p, _bucket_t *pptr) { \
		/* content the same as rbt_insert_update_core */ \
		uint32_t const lidx = _rbt_idx(_hdr(pptr)->children[0]); \
		uint32_t const ridx = _rbt_idx(_hdr(pptr)->children[1]); \
		_bucket_t *lptr = lidx == 0 ? NULL : _rbt_ptr(_bucket_t, p, lidx); \
		_bucket_t *rptr = ridx == 0 ? NULL : _rbt_ptr(_bucket_t, p, ridx); \
		_unused(lptr);		/* guard for nop */ \
		_unused(rptr); \
		return(_update(pptr, lptr, rptr));			/* might be nop; then pointers are unused */ \
	} \
	static _force_inline \
	void rbt_patch_##_sfx(rbt_iter_t const *w, _bucket_t *arr) { \
		uint32_t *p = (uint32_t *)arr; \
		/* move upward until hit root */ \
		uint32_t const *b = w->b; \
		while(--b > &w->ibuf[0]) { \
			_bucket_t *x = _rbt_ptr(_bucket_t, p, _rbt_idx(*b)); \
			if(!rbt_patch_core_##_sfx(p, x)) { break; } \
		} \
		return; \
	}


#define RBT_INIT_ITER(_sfx, _bucket_t, _hdr, _cmp_head, _cmp_tail, _cmp_iter) 	RBT_INIT_ITER_PATCH(_sfx, _bucket_t, _hdr, _cmp_head, _cmp_tail, _cmp_iter, _rbt_nop)


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

		rbt_iter_t it;
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
