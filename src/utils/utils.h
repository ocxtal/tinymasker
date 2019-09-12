
/**
 * @file utils.h
 * @brief collection of utils
 */

#ifndef _UTILS_H_INCLUDED
#define _UTILS_H_INCLUDED


/* include global header *before* we include individual dependencies */
#include "common.h"

/* arch/ needed for SIMD and arch-dependent macros */
#include "arch.h"


/* debugging and logging */
#include "bench.h"
#include "debug.h"
#include "log.h"

/* memory management */
#include "wmalloc.h"
#include "lmm.h"

/* testing */
#include "unittest.h"

/* data structures */
#include "ksort.h"
#include "kvec.h"
#include "rbt.h"
#include "rhash.h"

/* supplementary utils */
#include "mmstring.h"
#include "xorshift.h"
#include "opt.h"

/* I/O */
#include "pstream.h"
#include "rbread.h"


#endif	/* _UTILS_H_INCLUDED */

/**
 * end of utils.h
 */
