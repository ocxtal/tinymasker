
/**
 * @file xorshift.h
 * @brief xorshift pseudo-random generator
 */

#ifndef _XORSHIFT_H_INCLUDED
#define _XORSHIFT_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"


typedef struct xorshift_s { uint32_t v[4]; } xorshift_t;

static inline
void xorshift_init_static(xorshift_t *x, uint32_t seed) {
	x->v[0] = seed<<13;
	x->v[1] = (seed>>9) ^ (seed<<19);
	x->v[2] = ((seed>>9) ^ (seed<<19))>>7;
	x->v[3] = seed;
	return;
}
#define xorshift_init(_seed)		({ xorshift_t *x = calloc(1, sizeof(xorshift_t)); xorshift_init_static(x, _seed); x; })
#define xorshift_destroy_static(_x)	;
#define xorshift_destroy(_x)		{ free(_x); }

static inline
uint32_t xorshift_rand(xorshift_t *x) {
	uint32_t s = x->v[0], t = x->v[1], u = x->v[2], v = x->v[3];
	uint32_t p = s ^ (s<<11), q = (v ^ (v>>19)) ^ (p ^ (p>>8));
	x->v[0] = t; x->v[1] = u; x->v[2] = v; x->v[3] = q;
	return(t);
}

#endif
/**
 * end of xorshift.h
 */
