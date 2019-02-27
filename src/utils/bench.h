
/**
 * @file bench.h
 * @brief time measurement utilities
 */

#ifndef _BENCH_H_INCLUDED
#define _BENCH_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"

/* benchmarking macros */
#ifdef BENCH
#include <sys/time.h>


#ifdef __x86_64__

#include <x86intrin.h>

#define rdtscp() ({ \
	uint32_t x; \
	__rdtscp(&x); \
})
#define microbench(x) ({ \
	uint64_t ts = rdtscp(); \
	x; \
	uint64_t te = rdtscp(); \
	te - ts; \
})
#endif


typedef struct bench_s { struct timeval s; uint64_t a; } bench_t;
#define bench_init(b)		{ memset(&(b).s, 0, sizeof(struct timeval)); (b).a = 0; }
#define bench_start(b)		{ gettimeofday(&(b).s, NULL); }
#define bench_end(b) { \
	struct timeval _e; \
	gettimeofday(&_e, NULL); \
	(b).a += ( (_e.tv_sec  - (b).s.tv_sec ) * 1000000000 \
	         + (_e.tv_usec - (b).s.tv_usec) * 1000); \
}
#define bench_get(b)		( (b).a )

#else /* #ifdef BENCH */

/** disable bench */
struct _bench { uint64_t _unused; };
typedef struct _bench bench_t;
#define bench_init(b) 		;
#define bench_start(b) 		;
#define bench_end(b)		;
#define bench_get(b)		( 0ULL )

#endif /* #ifdef BENCH */

#endif
/**
 * end of bench.h
 */
