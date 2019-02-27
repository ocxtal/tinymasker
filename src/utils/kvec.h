/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:

#include "kvec.h"
int main() {
	kvec_t(int) array;
	kv_init(array);
	kv_push(int, array, 10); // append
	kv_a(int, array, 20) = 5; // dynamic
	kv_A(array, 20) = 4; // static
	kv_destroy(array);
	return 0;
}
*/

/*
  2008-09-22 (0.1.0):

	* The initial version.

*/

#ifndef _KVEC_H_INCLUDED
#define _KVEC_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"


#define kv_roundup32(x)				_roundup(x, 32)
#define kvec_t(type) 				struct { type *a; size_t n, m; }
#define kv_init(v)					((v).n = (v).m = 0, (v).a = 0)
#define kv_build(v, _p, _n, _m)		((v).n = (_n), (v).m = (_m), (v).a = (_p))
#define kv_inits(type)				((kvec_t(type)){ .n = 0, .m = 0, .a = NULL })
#define kv_destroy(v)				free((v).a)
#define kv_clear(v)					({ (v).n = 0; })
#define kv_A(v, i)					((v).a[(i)])
#define kv_pop(v)					((v).a[--(v).n])
#define kv_cnt(v)					((v).n)
#define kv_size(v)					( sizeof(*(v).a) * (v).n )
#define kv_max(v)					((v).m)
#define kv_ptr(v)					((v).a)
#define kv_tail(v)					( &(v).a[(v).n - 1] )

#define kv_resize(type, v, s) ({ \
	(v).m = (s); (v).a = realloc((v).a, sizeof(type) * (v).m); \
})
#define kv_expand(type, v, s1, s2) ({ \
	(v).m = kv_roundup32(MAX2((s1), (s2))); \
	kv_resize(type, v, (v).m); \
})
#define kv_reserve(type, v, s) ({ \
	if(_unlikely((v).m < (s))) { kv_expand(type, v, 2 * (v).m, (s)); } (v).a; \
})
#define kv_copy(type, v1, v0) ({ \
	if(_unlikely((v1).m < (v0).n)) { kv_expand(type, v1, (v0).n, 0); } \
	(v1).n = (v0).n; \
	memcpy((v1).a, (v0).a, sizeof(type) * (v0).n); \
}
#define kv_push(type, v, x) ({ \
	if(_unlikely((v).n >= (v).m)) { kv_expand(type, v, 2 * (v).m, 2); } \
	(v).a[(v).n++] = (x); \
	(v).n - 1; \
})
#define kv_pushp(type, v) ({ \
	if(_unlikely((v).n >= (v).m)) { kv_expand(type, v, 2 * (v).m, 2); } \
	&(v).a[(v).n++]; \
})
#define kv_pushm(type, v, arr, cnt) ({ \
	if(_unlikely(((v).m - (v).n) < (size_t)(cnt))) { \
		kv_expand(type, v, 2 * (v).m, (v).n + (cnt)); \
	} \
	memcpy(&(v).a[(v).n], (arr), (cnt) * sizeof(type)); \
	(v).n += (size_t)(cnt); \
	(v).n - (size_t)(cnt); \
})
#define kv_pushmp(type, v, arr, cnt) ({ \
	if(_unlikely(((v).m - (v).n) < (size_t)(cnt))) { \
		kv_expand(type, v, 2 * (v).m, (v).n + (cnt)); \
	} \
	(v).n += (size_t)(cnt); \
	&(v).a[(v).n - (size_t)(cnt)]; \
})
#define kv_a(type, v, i) ((v).m <= (size_t)(i)?						\
						  ((v).m = (v).n = (i) + 1, kv_roundup32((v).m), \
						   (v).a = (type*)realloc((v).a, sizeof(type) * (v).m), 0) \
						  : (v).n <= (size_t)(i)? (v).n = (i)			\
						  : 0), (v).a[(i)]

#define kv_reverse(type, v, start) do { \
	if ((v).m > 0 && (v).n > (start)) { \
		size_t __i, __end = (v).n - (start); \
		type *__a = (v).a + (start); \
		for (__i = 0; __i < __end>>1; ++__i) { \
			type __t = __a[__end - 1 - __i]; \
			__a[__end - 1 - __i] = __a[__i]; __a[__i] = __t; \
		} \
	} \
} while (0)
#define kv_foreach(type, v, _body) { \
	type *p = (type *)(v).a - 1; \
	type *_t = (type *)(v).a + (v).n; \
	while(++p < _t) { _body; } \
}

#define kv_dump(type, v, filename, ...) { \
	FILE *_fp = fopen(filename, "w"); \
	kv_foreach(type, v, ({ fprintf(_fp, __VA_ARGS__); })); \
	fclose(_fp); \
}

/* margined vector */
#define kvecm_t(type)				struct { type *a; size_t n, m; uint32_t h, t; }
#define kvm_init(v, _h, _t)			({ (v).n = (v).m = 0; (v).a = NULL; (v).h = (_h); (v).t = (_h) + (_t); 0; })
#define kvm_destroy(v)				{ if((v).a) { free(_sub_offset((v).a, (v).h)); } }
#define kvm_clear(v)				({ (v).n = 0; })
#define kvm_pop(v)					kv_pop(v)
#define kvm_cnt(v)					kv_cnt(v)
#define kvm_size(v)					kv_size(v)
#define kvm_max(v)					kv_max(v)
#define kvm_ptr(v)					( (v).a )
#define kvm_base_ptr(v)				( _sub_offset((v).a, (v).h) )
#define kvm_base_size(v)			( sizeof(*(v).a) * (v).n + (v).t )


#define kvm_resize(type, v, s) ({ \
	(v).m = (s); \
	(v).a = (type *)_add_offset(realloc((v).a == NULL ? NULL : _sub_offset((v).a, (v).h), sizeof(type) * (v).m + (v).t), (v).h); \
})
#define kvm_expand(type, v, s1, s2) ({ \
	(v).m = kv_roundup32(MAX2((s1), (s2))); \
	kvm_resize(type, v, (v).m); \
})
#define kvm_reserve(type, v, s) ({ \
	if(_unlikely((v).m < (s))) { kvm_expand(type, v, 2 * (v).m, (s)); } 0; \
})
#define kvm_push(type, v, x) ({ \
	if(_unlikely((v).n >= (v).m)) { kvm_expand(type, v, 2 * (v).m, 2); } \
	(v).a[(v).n++] = (x); \
	(v).n - 1; \
})
#define kvm_pushp(type, v) ({ \
	if(_unlikely((v).n >= (v).m)) { kvm_expand(type, v, 2 * (v).m, 2); } \
	&(v).a[(v).n++]; \
})
#define kvm_pushm(type, v, arr, size) ({ \
	if(_unlikely(((v).m - (v).n) < (size_t)(size))) { \
		kvm_expand(type, v, 2 * (v).m, (v).n + (size)); \
	} \
	memcpy(&(v).a[(v).n], (arr), (size) * sizeof(type)); \
	(v).n += (size_t)(size); \
	(v).n - (size_t)(size); \
})
#define kvm_fill_margin(type, v, c) { \
	memset(_sub_offset((v).a, (v).h), (c), (v).h); \
	memset((void *)&(v).a[(v).n], (c), (v).t - (v).h); \
}
#define kvm_foreach(type, v, _b)	kv_foreach(type, v, _b)

/** heap queue : elements in v must be orderd in heap */
#define kvechq_t(type)		struct { type *a; size_t n, m, offset; }
#define kv_hq_init(v)		({ (v).n = 0; (v).m = 0; (v).offset = 0; (v).a = NULL; (v); })
#define kv_hq_destroy(v)	kv_destroy(v)
#define kv_hq_clear(v)		( (v).offset = 0, (v).n = 0 )
#define kv_hq_cnt(v)		( kv_cnt(v) - (v).offset )
#define kv_hq_size(v)		( kv_size(v) )

#define kv_hq_head(v)		( (v).a[(v).offset] )
#define kv_hq_push(type, __comp, v, x) { \
	kv_push(type, v, x); \
	/* debug("pushed at i(%lu)", (v).n - 1); */ \
	size_t _ofs = ((size_t)(v).offset) - 1, i = (v).n - (v).offset; \
	while(i > 1 && __comp((v).a[(i>>1) + _ofs], (v).a[i + _ofs]) > 0) { \
		/* debug("swap, i(%lu, %lu)", (i>>1) + _ofs, i + _ofs); */ \
		type _tmp = (v).a[(i>>1) + _ofs]; \
		(v).a[(i>>1) + _ofs] = (v).a[i + _ofs]; \
		(v).a[i + _ofs] = _tmp; \
		i >>= 1; \
	} \
}
#define kv_hq_pop(type, __comp, v) ({ \
	size_t _ofs = ((size_t)(v).offset) - 1, i = 1, j = 2; \
	type _popped = (v).a[i + _ofs]; \
	(v).a[i + _ofs] = (v).a[(v).n + _ofs]; (v).n--; \
	while(j <= (v).n) { \
		/*k = (j + 1 < (v).n && kv_hq_n(v, j + 1) < kv_hq_n(v, j)) ? (j + 1) : j; */ \
		size_t k = (j + 1 <= (v).n && __comp((v).a[j + 1 + _ofs], (v).a[j + _ofs]) < 0) ? j + 1 : j; \
		/*k = (kv_hq_n(v, k) < kv_hq_n(v, i)) ? k : 0; */ \
		if(__comp((v).a[k + _ofs], (v).a[i + _ofs]) >= 0) { break; } \
		type _tmp = (v).a[k + _ofs]; \
		(v).a[k + _ofs] = (v).a[i + _ofs]; \
		(v).a[i + _ofs] = _tmp; \
		i = k; j = k<<1; \
	} \
	_popped; \
})
#define kv_hq_offset(v)		( (v).offset )
#define kv_hq_flush(v)		( (v).offset = (v).n )

#endif
/**
 * end of kvec.h
 */
