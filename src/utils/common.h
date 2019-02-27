
/**
 * @file common.h
 */

#ifndef _COMMON_H_INCLUDED
#define _COMMON_H_INCLUDED

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif
#if defined(__darwin__) && !defined(_DARWIN_C_FULL)
#  define _DARWIN_C_SOURCE		_DARWIN_C_FULL
#endif


/* use sched.h and madvise */
#if defined(__linux__)
#  define _GNU_SOURCE
#  define _DEFAULT_SOURCE
#  include <features.h>
#endif


#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>			/* size_t, ptrdiff_t, ... */


/* miscellaneous */
/* max, min, roundup */
#define MAX2(x,y)			 	( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 				( (x) < (y) ? (x) : (y) )
#define MAX3(x,y,z)				MAX2(x, MAX2(y, z))
#define MIN3(x,y,z)				MIN2(x, MIN2(y, z))
#define MAX4(w,x,y,z)			MAX2(MAX2(w, x), MAX2(y, z))
#define MIN4(w,x,y,z)			MIN2(MIN2(w, x), MIN2(y, z))
#define _roundup(x, base)		( ((x) + (base) - 1) & ~((base) - 1) )
#define _cutdown(x, base)		( (x) & ~((base) - 1) )

#define _pack_ptr_flag(p, b)	( (((intptr_t)(p))<<1) | (((intptr_t)(b)) & 0x01) )
#define _unpack_ptr(u)			( (void *)(((intptr_t)(u))>>1) )
#define _unpack_flag(u)			( (uint64_t)(((intptr_t)(u)) & 0x01) )

#define _add_offset(p, ofs)		( (void *)(((uint8_t *)(p)) + (ptrdiff_t)((size_t)(ofs))) )
#define _sub_offset(p, ofs)		( (void *)(((uint8_t *)(p)) - (ptrdiff_t)((size_t)(ofs))) )
#define _ptr_diff(p, q)			( (size_t)((uint8_t *)(p) - (uint8_t *)(q)) )

/* _likely, _unlikely */
#define _likely(x)				__builtin_expect(!!(x), 1)
#define _unlikely(x)			__builtin_expect(!!(x), 0)

/* _unused */
#define _unused(x)				{ (void)(x); }

/* _force_inline */
#define _force_inline	inline

/**
 * static assertion macros
 */
#define _sa_cat_intl(x, y)		x##y
#define _sa_cat(x, y)			_sa_cat_intl(x, y)
#define _static_assert(expr)	typedef char _sa_cat(_st, __LINE__)[(expr) ? 1 : -1]

/* misc vector types */
typedef union { uint64_t u64[1]; uint32_t u32[2]; } v2u32_t;
typedef struct { uint32_t u32[3]; } v3u32_t;
typedef union { uint64_t u64[2]; uint32_t u32[4]; } v4u32_t;
typedef struct { uint64_t u64[2]; } v2u64_t;
_static_assert(sizeof(v2u32_t) == 8);
_static_assert(sizeof(v3u32_t) == 12);
_static_assert(sizeof(v4u32_t) == 16);
_static_assert(sizeof(v2u64_t) == 16);

typedef struct { v2u32_t *a; size_t n, m; } v2u32_v;
typedef struct { v3u32_t *a; size_t n, m; } v3u32_v;
typedef struct { v4u32_t *a; size_t n, m; } v4u32_v;
typedef struct { uint64_t *a; size_t n, m; } uint64_v;
typedef struct { uint32_t *a; size_t n, m; } uint32_v;
typedef struct { uint16_t *a; size_t n, m; } uint16_v;
typedef struct { uint8_t *a; size_t n, m; } uint8_v;
typedef struct { uint8_t *a; size_t n, m; uint32_t h, t; } uint8m_v;		/* margined */
typedef struct { void **a; size_t n, m; } ptr_v;

/**
 * @type read_t, write_t
 * @brief abstract I/O function types for rhash, rbt, and pstream
 */
typedef size_t (*read_t)(void *fp, void *buf, size_t size);
typedef size_t (*write_t)(void *fp, void const *buf, size_t size);

#endif
/**
 * end of common.h
 */
