/**
 * @file rhash.h
 * @brief Robinhood hashmap; 64bit key -> arbitrary sized bucket and string hashmap
 */

#ifndef _RHASH_H_INCLUDED
#define _RHASH_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"
#include "arch.h"				/* for _lzc_u64 */


/* rhash.h: Robinhood hash table for integers (much faster than khash.h) */

/* common constants */
#define RH_INIT_SIZE			( 256 )				/* initial table size */
#define RH_THRESH				( 0.85 )			/* max occupancy */
#define RH_DST_MAX				( 48 )				/* max distance from base pos */
#define RH_INIT_VAL				( UINT64_MAX )		/* initial vals */
#define RH_UNTOUCHED			( UINT64_MAX )
#define RH_MOVED				( UINT64_MAX - 1 )
#define RH_FAILED				( UINT32_MAX )
_static_assert(2 * RH_DST_MAX < RH_INIT_SIZE);

/* common types */
// typedef struct rh_bidx_s { size_t idx, n; } rh_bidx_t;
typedef struct rh_bidx_s { size_t idx, cnt, poll_len; } rh_bidx_t;
typedef struct rh_hdr_s { uint32_t size, cnt; } rh_hdr_t;
_static_assert(sizeof(rh_hdr_t) == 8);

#define rh_is_untouched(_key_t, _x)					( (_x) == (_key_t)(RH_UNTOUCHED) )
#define rh_is_moved(_key_t, _x)						( (_x) == (_key_t)(RH_MOVED) )
#define rh_is_empty(_key_t, _x, _msb_mask)			( (_x) >= (_key_t)(RH_MOVED) )
#define rh_is_empty_mask(_key_t, _x, _msb_mask)		( ((_x) & (_msb_mask)) != 0 )
#define rh_cmp(_key_t, _x, _y, _val_mask)			( (_key_t)(_x) - (_key_t)(_y) )
#define rh_cmp_mask(_key_t, _x, _y, _val_mask)		( ((_key_t)(_x) - (_key_t)(_y)) & (_val_mask) )

#define RH_INIT_INTL(_sfx, _bkt_t, _key_t, _key, _val_t, _val, _is_empty, _msb_mask, _cmp, _val_mask) \
	/* size of the array, element count, upper bound */ \
	typedef struct rh_##_sfx##_s { uint32_t mask, cnt; _bkt_t *a; uint32_t ub, max; } rh_##_sfx##_t; \
	/* misc */ \
	static _force_inline _bkt_t const *rh_ptr_##_sfx(rh_##_sfx##_t const *h) { return(h->a); } \
	static _force_inline _bkt_t const *rh_tail_ptr_##_sfx(rh_##_sfx##_t const *h) { return(&h->a[h->mask + 1]); } \
	static _force_inline size_t rh_size_##_sfx(rh_##_sfx##_t const *h) { return(sizeof(*h->a) * (h->mask + 1)); } \
	static _force_inline size_t rh_max_##_sfx(rh_##_sfx##_t const *h) { return(h->mask + 1); } \
	static _force_inline size_t rh_cnt_##_sfx(rh_##_sfx##_t const *h) { return(h->cnt); } \
	/* init and destroy */ \
	static _force_inline \
	void rh_clear_arr_##_sfx(_bkt_t *a, size_t len) { \
		/* init table with invalid key-val pairs */ \
		memset(a, 0xff, len * sizeof(_bkt_t)); \
	} \
	static _force_inline \
	void rh_init_static_##_sfx(rh_##_sfx##_t *h, size_t size) { \
		/* roundup to power of two */ \
		size = 0x8000000000000000>>(_lzc_u64(size - 1) - 1); \
		size = MAX2(size, RH_INIT_SIZE); \
		/* initialize kh object */ \
		*h = (rh_##_sfx##_t){ \
			.mask = size - 1,		/* in-use table size */ \
			.max = size,			/* malloc'd table size */ \
			.cnt = 0, \
			.ub = size * RH_THRESH, \
			.a = malloc(sizeof(_bkt_t) * size) \
		}; \
		rh_clear_arr_##_sfx(h->a, size); \
		return; \
	} \
	static _force_inline \
	rh_##_sfx##_t *rh_init_##_sfx(size_t size) { \
		rh_##_sfx##_t *h = calloc(1, sizeof(rh_##_sfx##_t)); rh_init_static_##_sfx(h, size); \
		return(h); \
	} \
	static _force_inline \
	void rh_destroy_static_##_sfx(rh_##_sfx##_t *h) { \
		if(h != NULL) { free(h->a); } return; \
	} \
	static _force_inline \
	void rh_destroy_##_sfx(rh_##_sfx##_t *h) { \
		rh_destroy_static_##_sfx(h); free(h); return; \
	} \
	static _force_inline \
	void rh_resize_##_sfx(rh_##_sfx##_t *h, size_t size) { \
		h->max = MAX2(size, RH_INIT_SIZE); \
		h->a = realloc(h->a, sizeof(_bkt_t) * h->max); \
	} \
	static _force_inline \
	void rh_dump_##_sfx(rh_##_sfx##_t const *h, void *fp, write_t wfp) { \
		rh_hdr_t hdr = { 0 }; \
		/* dump a mark of zero-sized table */ \
		if(h == NULL || h->a == NULL) { wfp(fp, &hdr, sizeof(rh_hdr_t)); return; } \
		/* dump size */ \
		hdr = (rh_hdr_t){ .size = h->mask + 1, .cnt = h->cnt };			/* table size and occupancy */ \
		wfp(fp, &hdr, sizeof(rh_hdr_t)); \
		wfp(fp, h->a, sizeof(_bkt_t) * hdr.size);		/* dump content */ \
		return; \
	} \
	static _force_inline \
	void rh_load_static_##_sfx(rh_##_sfx##_t *h, void *fp, read_t rfp) { \
		rh_hdr_t hdr = { 0 }; \
		/* read sizes, return NULL if table size is zero */ \
		if((rfp(fp, &hdr, sizeof(rh_hdr_t))) != sizeof(rh_hdr_t) || hdr.size == 0) { \
			*h = (rh_##_sfx##_t){ 0 }; return; \
		} \
		/* create hash object */ \
		*h = (rh_##_sfx##_t){ \
			.mask = hdr.size - 1,						/* in-use table size */ \
			.max = hdr.size,							/* malloc'd table size */ \
			.cnt = hdr.cnt, \
			.ub = hdr.size * RH_THRESH, \
			.a = malloc(sizeof(_bkt_t) * hdr.size) \
		}; \
		/* read table */ \
		if((rfp(fp, h->a, sizeof(_bkt_t) * hdr.size)) != sizeof(_bkt_t) * hdr.size) { \
			free(h->a); h->a = NULL; return;				/* error happened */ \
		} \
		return; \
	} \
	static _force_inline \
	rh_##_sfx##_t *rh_load_##_sfx(void *fp, read_t rfp) { \
		rh_##_sfx##_t *h = calloc(1, sizeof(rh_##_sfx##_t)); \
		rh_load_static_##_sfx(h, fp, rfp); \
		return(h); \
	} \
	static _force_inline \
	void rh_clear_##_sfx(rh_##_sfx##_t *h) { \
		if(h == 0) { return; } \
		/* clear hash table; don't clear max since it holds the malloc'd table size */ \
		h->mask = RH_INIT_SIZE - 1; \
		h->cnt = 0; \
		h->ub = RH_INIT_SIZE * RH_THRESH; \
		rh_clear_arr_##_sfx(h->a, RH_INIT_SIZE); \
		return; \
	} \
	static _force_inline \
	rh_bidx_t rh_reallocate_##_sfx(_bkt_t *a, _key_t k0, _val_t v0, uint64_t mask, uint64_t force_save) { \
		uint64_t const kmask = mask | ((uint64_t)(((_key_t)RH_MOVED)>>1) + 1); \
		size_t const min_bin_dist = -(RH_DST_MAX + 1ULL); \
		size_t max_poll_len = 0; \
		/* nb: base index; ni: offset from base index */ \
		size_t nb = k0 & mask, ni = 0; \
		/* initialize working variables */ \
		size_t b = nb, i = ni - 1; \
		debug("find bin for k(%lx), v(%lx), start from b(%lx)", k0, v0, b); \
		uint64_t k1; \
		do { \
			i++; \
			k1 = (uint64_t)_key(&a[(b + i) & mask]); \
			debug("test k1(%lx) at i(%lx)", k1, (b + i) & mask); \
			if(_unlikely(_is_empty(_key_t, k1, _msb_mask))) { goto _rh_reallodate_found_new; } \
			if(_unlikely(k0 == k1)) { goto _rh_reallocate_found_replace; } \
			/* if the origin of the polled key is larger than the key to be inserted */ \
		} while(b - (k1 & kmask) < min_bin_dist); \
		/* save bin location */ \
		size_t idx = (b + i) & mask; \
		do { \
			uint64_t v1 = (uint64_t)_val(&a[(b + i) & mask]); \
			debug("found swappable bin for k(%lx), v(%lx) at i(%lx), pushed out k(%lx), v(%lx)", k0, v0, (b + i) & mask, k1, v1); \
			_key(&a[(b + i) & mask]) = (_key_t)k0; \
			_val(&a[(b + i) & mask]) = (_val_t)v0; \
			k0 = k1; v0 = v1; \
			/* update max polled length; skip for to-be-moved key */ \
			max_poll_len = MAX2(max_poll_len, i); \
			/* calculate next polling base */ \
			nb = k0 & mask; \
			ni = (b + i + 1 - nb) & mask; \
			/* reset working variables */ \
			b = nb; \
			i = ni - 1; \
			/* load the next key */ \
			do { \
				i++; \
				k1 = (uint64_t)_key(&a[(b + i) & mask]); \
				debug("test k1(%lx) at i(%lx)", k1, (b + i) & mask); \
				if(_unlikely(_is_empty(_key_t, k1, _msb_mask))) { goto _rh_reallodate_found_tail; } \
				if(_unlikely(k0 == k1)) { goto _rh_reallocate_found_replace; } \
			} while(b - (k1 & kmask) < min_bin_dist); \
		} while(1); \
	_rh_reallodate_found_new:; { \
			idx = (b + i) & mask; \
		_rh_reallodate_found_tail:; \
			debug("found last bin for k(%lx), v(%lx) at i(%lx); idx(%lx), max_poll_len(%lu), save(%u)", k0, v0, (b + i) & mask, idx, MAX2(max_poll_len, i), force_save | (k0 != k1)); \
			/* save the last key-value pair if the bin is newly allocated one */ \
			_key(&a[(b + i) & mask]) = (_key_t)k0; \
			_val(&a[(b + i) & mask]) = (_val_t)v0; \
			return((rh_bidx_t){ \
				.idx = idx, \
				.cnt = 1, \
				.poll_len = MAX2(max_poll_len, i)	/* max polled length */ \
			}); \
		} \
	_rh_reallocate_found_replace:; { \
			debug("found replacable bin for k(%lx), v(%lx) at i(%lx); max_poll_len(%lu), save(%u)", k0, v0, (b + i) & mask, MAX2(max_poll_len, i), force_save | (k0 != k1)); \
			/* b and i correspond to key when jumped here */ \
			size_t const nidx = (b + i) & mask; \
			if(force_save) { \
				_key(&a[nidx]) = (_key_t)k0; \
				_val(&a[nidx]) = (_val_t)v0; \
			} \
			return((rh_bidx_t){ \
				.idx = nidx, \
				.cnt = 0, \
				.poll_len = MAX2(max_poll_len, i) \
			}); \
		} \
	} \
	static _force_inline \
	void rh_extend_##_sfx(rh_##_sfx##_t *h) { \
		size_t max_poll_len = 0; \
		do { \
			size_t const prev_size = h->mask + 1, size = 2 * prev_size; \
			uint64_t const mask = size - 1; \
			/* debug("extend, prev_size(%lu), cnt(%u), ub(%u), max(%u)", prev_size, h->cnt, h->ub, h->max); */ \
			h->mask = mask; h->ub = size * RH_THRESH;		/* update size */ \
			if(size > h->max) {								/* double the table if needed */ \
				h->a = realloc(h->a, sizeof(_bkt_t) * size); h->max = size; \
			} \
			/* clear the extended area */ \
			rh_clear_arr_##_sfx(h->a + prev_size, prev_size); \
			/* rehash */ \
			max_poll_len = 0; \
			for(size_t i = 0; i < size; i++) { \
				uint64_t const k0 = (uint64_t)_key(&h->a[i]);		/* load key (with val if compacted) */ \
				if(_is_empty(_key_t, k0, _msb_mask) || (k0 & mask) == i) { continue; }	/* test if rehashing is required */ \
				uint64_t const v0 = (uint64_t)_val(&h->a[i]); \
				/* clear the current bin before re-inserting the key */ \
				_key(&h->a[i]) = RH_INIT_VAL; \
				_val(&h->a[i]) = RH_INIT_VAL; \
				 debug("k(%lx), v(%lx) at i(%lx) needs rehashing", k0, v0, i);  \
				/* re-insert the key; search an appropriate bin */ \
				rh_bidx_t const b = rh_reallocate_##_sfx(h->a, k0, v0, mask, 1); \
				max_poll_len = MAX2(max_poll_len, b.poll_len & (mask>>1)); \
			} \
			/* fprintf(stderr, "extend, h(%p), size(%zx), max_poll_len(%zu)\n", h, size, max_poll_len); */ \
		} while(max_poll_len >= RH_DST_MAX); \
		return; \
	} \
	static _force_inline \
	void rh_put_##_sfx(rh_##_sfx##_t *h, _key_t k, _val_t v) { \
		/* debug("search bin for k(%lx)", k); */ \
		while(1) { \
			if(h->cnt < h->ub) { \
				rh_bidx_t const b = rh_reallocate_##_sfx(h->a, k, v, h->mask, 1);	/* allocate bin for the new key */ \
				if(b.poll_len < RH_DST_MAX) { h->cnt += b.cnt; return; } \
			} \
			/* debug("failed to allocate bin, key(%lx), cnt(%u), max(%u)", k, h->cnt, h->max); trap(); */ \
			rh_extend_##_sfx(h);		/* extend table when cnt exceeded upper bound or empty bin not found within MAX_DST from base */ \
		} \
		return; \
	} \
	static _force_inline \
	_bkt_t *rh_put_ptr_##_sfx(rh_##_sfx##_t *h, _key_t k) { \
		/* debug("search bin for k(%lx)", k); */ \
		while(1) { \
			/* size_t poll_len = 0; */ \
			if(h->cnt < h->ub) { \
				rh_bidx_t const b = rh_reallocate_##_sfx(h->a, k, (_val_t)RH_INIT_VAL, h->mask, 0);		/* allocate bin for the new key */ \
				if(b.poll_len < RH_DST_MAX) { h->cnt += b.cnt; return(&h->a[b.idx]); }			/* &h->a[b.idx].u64[1] */ \
				/* poll_len = b.poll_len; */ \
			} \
			/* debug("failed to allocate bin, key(%lx), cnt(%u), max(%u)", k, h->cnt, h->max); trap(); */ \
			/* fprintf(stderr, "try extend, cnt(%u), ub(%u), poll_len(%zu)\n", h->cnt, h->ub, poll_len); */ \
			rh_extend_##_sfx(h);		/* extend table when cnt exceeded upper bound or empty bin not found within MAX_DST from base */ \
		} \
		return(NULL); \
	} \
	static _force_inline \
	_val_t rh_get_##_sfx(rh_##_sfx##_t const *h, _key_t key) { \
		uint64_t mask = h->mask, k = _key(&h->a[key & mask]); \
		size_t b = key & mask, i = 0; \
		do { \
			if(_cmp(_key_t, k, key, _val_mask) == 0) { return(_val(&h->a[(b + i) & mask])); } \
			i++; k = _key(&h->a[(b + i) & mask]); \
		} while(i <= RH_DST_MAX && !rh_is_untouched(_key_t, k));		/* !rh_is_untouched(k) || rh_is_moved(k) */ \
		return((_val_t)RH_INIT_VAL);									/* not found */ \
	} \
	static _force_inline \
	_bkt_t *rh_get_ptr_##_sfx(rh_##_sfx##_t const *h, _key_t key) { \
		uint64_t mask = h->mask, k = _key(&h->a[key & mask]); \
		size_t b = key & mask, i = 0; \
		do { \
			if(_cmp(_key_t, k, key, _val_mask) == 0) { return(&h->a[(b + i) & mask]); } \
			i++; k = _key(&h->a[(b + i) & mask]); \
		} while(i <= RH_DST_MAX && !rh_is_untouched(_key_t, k));		/* !rh_is_untouched(k) || rh_is_moved(k) */ \
		return(NULL); \
	}

/* wrappers */
#define RH_INIT(_sfx, _bkt_t, _key_t, _key, _val_t, _val) \
	RH_INIT_INTL(_sfx, _bkt_t, _key_t, _key, _val_t, _val, rh_is_empty, 0ULL, rh_cmp, 0ULL)
#define RH_INIT_MASK(_sfx, _bkt_t, _key_t, _key, _val_t, _val, _empty_mask, _val_mask) \
	RH_INIT_INTL(_sfx, _bkt_t, _key_t, _key, _val_t, _val, rh_is_empty, _empty_mask, rh_cmp_mask, _val_mask)


/* freezed (compacted) */
#define RH_INIT_FREEZE(_sfx, _bkt_t, _key_t, _val_t) \
	/* size of the array, element count, upper bound */ \
	typedef struct rh_fz_##_sfx##_s { uint32_t mask, cnt; _bkt_t *a; } rh_fz_##_sfx##_t; \
	/* convert */ \
	rh_fz_##_sfx##_t rh_freeze_##_sfx(rh_##_sfx##_t const *h) { \
		return((rh_fz_##_sfx##_t){ .mask = h->mask, .cnt = h->cnt, .a = h->a, }); \
	} \
	/* misc */ \
	static _force_inline _bkt_t const *rh_fz_ptr_##_sfx(rh_fz_##_sfx##_t const *h) { return(rh_ptr_##_sfx((rh_##_sfx##_t const *)(h))); } \
	static _force_inline _bkt_t const *rh_fz_tail_ptr_##_sfx(rh_fz_##_sfx##_t const *h) { return(rh_tail_ptr_##_sfx((rh_##_sfx##_t const *)(h))); } \
	static _force_inline size_t rh_fz_size_##_sfx(rh_fz_##_sfx##_t const *h) { return(rh_size_##_sfx((rh_##_sfx##_t const *)(h))); } \
	static _force_inline size_t rh_fz_max_##_sfx(rh_fz_##_sfx##_t const *h) { return(rh_max_##_sfx((rh_##_sfx##_t const *)(h))); } \
	static _force_inline size_t rh_fz_cnt_##_sfx(rh_fz_##_sfx##_t const *h) { return(rh_cnt_##_sfx((rh_##_sfx##_t const *)(h))); } \
	/* destroy */ \
	static _force_inline \
	void rh_fz_destroy_static_##_sfx(rh_fz_##_sfx##_t *h) { \
		if(h != NULL) { free(h->a); } return; \
	} \
	static _force_inline \
	void rh_fz_destroy_##_sfx(rh_fz_##_sfx##_t *h) { \
		rh_fz_destroy_static_##_sfx(h); free(h); return; \
	} \
	/* I/O for fz object */ \
	static _force_inline \
	void rh_fz_dump_##_sfx(rh_fz_##_sfx##_t const *h, void *fp, write_t wfp) { \
		rh_dump_##_sfx((rh_##_sfx##_t const *)h, fp, wfp); \
		return; \
	} \
	static _force_inline \
	void rh_fz_load_static_##_sfx(rh_fz_##_sfx##_t *h, void *fp, read_t rfp) { \
		rh_##_sfx##_t rh;			/* on stack */ \
		rh_load_static_##_sfx(&rh, fp, rfp); \
		*h = rh_freeze_##_sfx(&rh); \
		return; \
	} \
	static _force_inline \
	rh_fz_##_sfx##_t *rh_fz_load_##_sfx(void *fp, read_t rfp) { \
		rh_fz_##_sfx##_t *h = calloc(1, sizeof(rh_fz_##_sfx##_t)); \
		rh_fz_load_static_##_sfx(h, fp, rfp); \
		return(h); \
	} \
	/* get */ \
	static _force_inline \
	_val_t rh_fz_get_##_sfx(rh_fz_##_sfx##_t const *h, _key_t key) { \
		return(rh_get_##_sfx((rh_##_sfx##_t const *)(h), key)); \
	} \
	static _force_inline \
	_bkt_t *rh_fz_get_ptr_##_sfx(rh_fz_##_sfx##_t const *h, _key_t key) { \
		return(rh_get_ptr_##_sfx((rh_##_sfx##_t const *)(h), key)); \
	}



#define rh_foreach_mask(_bkt_t, _key_t, _key, _rh, _is_empty, _msb_mask, _body) { \
	_bkt_t *base = (_rh)->a; \
	for(size_t _i = 0; _i < (_rh)->mask + 1; _i++) { \
		_bkt_t *p = &base[_i]; \
		uint64_t key = (uint64_t)_key(p); \
		if(_is_empty(_key_t, key, _msb_mask)) { continue; } \
		{ _body; } \
	} \
}
#define rh_foreach(_bkt_t, _key_t, _key, _rh, _body)	rh_foreach_mask(_bkt_t, _key_t, _key, _rh, rh_is_empty, 0ULL, _body)

/**
 * @struct rh_str_t
 * @brief string -> string hashmap
 */
typedef struct rh_str_idx_bucket_s { uint64_t key; size_t ofs; } rh_str_idx_bucket_t;
#define _rh_key_str_idx(_p)				(_p)->key
#define _rh_val_str_idx(_p)				(_p)->ofs

RH_INIT(str_idx, rh_str_idx_bucket_t,
	uint64_t, _rh_key_str_idx,
	size_t, _rh_val_str_idx
);
typedef struct { uint8_v s; rh_str_idx_t h; } rh_str_t;

#define rh_str_cnt(_h)					( rh_cnt_str_idx(&((rh_str_t *)(_h))->h) )
#define rh_str_ptr(_h)					( rh_ptr_str_idx(&((rh_str_t const *)(_h))->h) )
#define rh_str_init_static(_h, _size)	{ (_h)->s = (uint8_v){ 0 }; rh_init_static_str_idx(&((rh_str_t *)(_h))->h, _size); }
#define rh_str_init(_size)				({ rh_str_t *h = calloc(1, sizeof(rh_str_t)); rh_init_static_str_idx(&(_h)->h, _size); h; })
#define rh_str_destroy_static(_h)		{ rh_destroy_static_str_idx(&((rh_str_t *)(_h))->h); free(((rh_str_t *)(_h))->s.a); }
#define rh_str_destroy(_h)				{ rh_destroy_static_str_idx(_h); free(_h); }
#define rh_str_clear(_h)				{ rh_clear_str_idx(&((rh_str_t *)(_h))->h); }

/* string hash */
static _force_inline
uint32_t rh_str_hash(char const *p, size_t l)
{
	uint64_t a = 0xcafe;
	while(*p && l) { a ^= (a<<5) ^ (a>>18) ^ *p++; l--; }
	return(a);
}

/**
 * @fn rh_str_put, rh_str_get
 */
static _force_inline
void rh_str_put(rh_str_t *h, char const *k, size_t klen, char const *v, size_t vlen)
{
	if(klen == 0) { klen = strlen(k); }
	if(vlen == 0) { vlen = strlen(v); }

	uint64_t key = ((uint64_t)rh_str_hash(k, klen)) | (klen<<32);
	while(1) {
		rh_str_idx_bucket_t *p = rh_put_ptr_str_idx(&h->h, key);
		if(p->ofs == (size_t)RH_INIT_VAL || strncmp(k, (char const *)kv_ptr(h->s) + p->ofs, klen) == 0) {
			/* bin found */
			p->ofs = kv_size(h->s);		/* overwrite offset; this never cause memory leak */
			kv_pushm(uint8_t, h->s, k, klen);
			kv_push(uint8_t, h->s, '\0');
			kv_pushm(uint8_t, h->s, v, vlen);
			kv_push(uint8_t, h->s, '\0');
			break;
		}
		key += 7;						/* collision detected; rehash */
	}
	return;
}
static _force_inline
char const *rh_str_get(rh_str_t const *h, char const *k, size_t klen)
{
	if(klen == 0) { klen = strlen(k); }

	uint64_t key = ((uint64_t)rh_str_hash(k, klen)) | (klen<<32);
	while(1) {
		rh_str_idx_bucket_t const *p = rh_get_ptr_str_idx(&h->h, key);
		if(p == NULL || p->ofs == (size_t)RH_INIT_VAL) { return(NULL); }	/* not found */
		if(strncmp(k, (char const *)kv_ptr(h->s) + p->ofs, klen) == 0) {
			return((char const *)kv_ptr(h->s) + (p->ofs + klen + 1));			/* found */
		}
		key += 7;		/* collision */
	}
	return(NULL);		/* never reach here */
}


#endif
/**
 * end of rhash.h
 */
