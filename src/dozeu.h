
// $(CC) -O3 -march=native -DMAIN -o dozeu dozeu.c
// #define DEBUG
// #define DZ_PRINT_VECTOR

#if !defined(_DZ_H_INCLUDED)
#define _DZ_H_INCLUDED


/* print_vector for debugging */
#ifdef DZ_PRINT_VECTOR
#define print_vector(v) { \
	debug("%s (%d, %d, %d, %d, %d, %d, %d, %d)", #v, \
	dz_rm_ofs(_mm_extract_epi16(v, 7)), \
	dz_rm_ofs(_mm_extract_epi16(v, 6)), \
	dz_rm_ofs(_mm_extract_epi16(v, 5)), \
	dz_rm_ofs(_mm_extract_epi16(v, 4)), \
	dz_rm_ofs(_mm_extract_epi16(v, 3)), \
	dz_rm_ofs(_mm_extract_epi16(v, 2)), \
	dz_rm_ofs(_mm_extract_epi16(v, 1)), \
	dz_rm_ofs(_mm_extract_epi16(v, 0))); \
}
#define print_vector_raw(v) { \
	debug("%s (%d, %d, %d, %d, %d, %d, %d, %d)", #v, \
	_mm_extract_epi16(v, 7), \
	_mm_extract_epi16(v, 6), \
	_mm_extract_epi16(v, 5), \
	_mm_extract_epi16(v, 4), \
	_mm_extract_epi16(v, 3), \
	_mm_extract_epi16(v, 2), \
	_mm_extract_epi16(v, 1), \
	_mm_extract_epi16(v, 0)); \
}

#define print_vector_8(v) { \
	debug("%s (%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", #v, \
	_mm_extract_epi8(v, 8 + 7), \
	_mm_extract_epi8(v, 8 + 6), \
	_mm_extract_epi8(v, 8 + 5), \
	_mm_extract_epi8(v, 8 + 4), \
	_mm_extract_epi8(v, 8 + 3), \
	_mm_extract_epi8(v, 8 + 2), \
	_mm_extract_epi8(v, 8 + 1), \
	_mm_extract_epi8(v, 8 + 0), \
	_mm_extract_epi8(v, 7), \
	_mm_extract_epi8(v, 6), \
	_mm_extract_epi8(v, 5), \
	_mm_extract_epi8(v, 4), \
	_mm_extract_epi8(v, 3), \
	_mm_extract_epi8(v, 2), \
	_mm_extract_epi8(v, 1), \
	_mm_extract_epi8(v, 0)); \
}


#else
#define print_vector(v) ;
#define print_vector_raw(v) ;
#define print_vector_8(v) ;

#endif


/**
 * @file dozeu.h
 * @brief SIMD X-drop DP for read-to-graph alignment
 *
 * @author Hajime Suzuki
 * @date 18/03/20 at Kyoto
 * @license MIT
 *
 * @datail
 * minimum requirements:
 *   SSE4.1 (Core 2 / Bobcat or later)
 *   { gcc, clang, icc } x { Linux, Mac OS }
 * input sequence encoding to be either of the three (not required to be terminated with '\0'):
 *   1. {'A','C','G','T','U','a','c','g','t','u'}      (ASCII nucleotide)
 *   2. { 0, 1, 2, 3 }                                 (2-bit nucleotide; optionally diagonal)
 *   3. { 0x01, 0x02, 0x04, 0x08 } or their union      (4-bit nucleotide)
 *   4. { 0, 1, 2, 3, ..., 19 }                        (integer protein)
 *
 * References:
 *   1. BLAST X-drop DP
 *     https://github.com/elucify/blast-docs/wiki/The-Developer's-Guide-to-BLAST (the most datailed document that describes the original X-drop DP algorithm)
 *   2. Myers bit vector (describes vertical (horizontal) tiling)
 *     Gene Myers, A fast bit-vector algorithm for approximate string matching based on dynamic programming, JACM (1999)
 */
#ifdef __cplusplus
extern "C" {
#endif

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif

#if defined(__darwin__) && !defined(_BSD_SOURCE)
#  define _BSD_SOURCE
#endif

// #include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <x86intrin.h>



/* default encoding: ascii */
#if !defined(DZ_WRAPPED_API) || (DZ_WRAPPED_API != 0)
#  if !defined(DZ_PROTEIN) && !defined(DZ_NUCL_ASCII) && !defined(DZ_NUCL_2BIT) && !defined(DZ_NUCL_4BIT) && !defined(DZ_NUCL_DIAG) && !defined(DZ_REALIGNER)
#    define DZ_NUCL_ASCII
#  endif
#endif



/* DZ_NUCL_DIAG suppose 2bit-encoded input */
#if defined(DZ_NUCL_DIAG)
#  define DZ_NONGRAPH
#endif
#if defined(DZ_NONGRAPH)
#  if defined(DZ_PROTEIN) || defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_4BIT) || defined(DZ_REALIGNER)
#    error "DZ_NONGRAPH is not compatibile with protein, nucleotide ASCII, nucleotide 4bit, or realigner."
#  endif
#endif

#if defined(DZ_NUCL_DIAG) && defined(DZ_FULL_LENGTH_BONUS)
#  error "DZ_NUCL_DIAG is not compatibile with full-length bonus."
#endif

#if defined(DZ_NUCL_DIAG) && defined(DZ_OFFSET)
#  error "DZ_NUCL_DIAG is not compatible with offset DP."
#endif



/* scoring and traceback configurations */
#ifndef DZ_CIGAR_OP
#  define DZ_CIGAR_OP				0x03040201
#endif
#define DZ_CIGAR_INTL				( (((uint32_t)DZ_CIGAR_OP)>>16) | (((uint32_t)DZ_CIGAR_OP)<<16) )

#ifndef dz_cmp_max
#  define dz_cmp_max(x, y)			( (x) > (y) )
#endif

#ifndef DZ_OVERFLOW_MARGIN
#  ifdef DEBUG
#    define DZ_OVERFLOW_MARGIN		( 256 )
#  else
#    define DZ_OVERFLOW_MARGIN		( 32768 - 256 )
#  endif
#endif



/* activate wrapped API (enabled by default) */
#ifndef DZ_WRAPPED_API
#  define DZ_WRAPPED_API			( 1 )
#endif



/* vectorization directive */
#ifndef __x86_64__
#  error "x86_64 is required"
#endif
#ifndef __SSE4_1__
#  warning "SSE4.1 is automatically enabled in dozeu.h, please check compatibility to the system."
#  define __dz_vectorize			__attribute__(( target( "sse4.1" ) )) inline
#else
#  define __dz_vectorize			inline		/* follow the compiler options */
#endif



/* inlining directive (FIXME: add appropriate __force_inline flag) */
#define __dz_force_inline			inline



/* warning suppressor */
#ifdef __cplusplus
#  define dz_unused(x)
#else
#  define dz_unused(x)				(void)(x)
#endif



/* unittest and debug configuration */
#if (defined(DZ_UNITTEST) && !defined(__cplusplus))
/* nothing */
#else
#  ifndef UNITTEST_H_INCLUDED
#    define unittest(...)			static void dz_pp_cat(dz_unused_, __LINE__)(void)
#    define ut_assert(...)			;
#  endif
#endif

#if !defined(_DEBUG_H_INCLUDED) && !defined(DEBUG)
#  define debug(...)				;
#  define debugblock(...)			;
#  define trap()					;
#endif



/*
 * scoring configuration
 */
typedef struct dz_score_conf_s {
	/*
	 * dimension is determined by macro at compile time;
	 * DZ_REF_MAT_SIZE for row length and DZ_QUERY_MAT_SIZE for #columns
	 */
	int8_t const *score_matrix;		/* negative for mismatch */

	/* gap penalties */
	uint16_t ins_open;				/* gap penalties in positive */
	uint16_t ins_extend;
	uint16_t del_open;				/* gap penalties in positive */
	uint16_t del_extend;

	/* end controls */
	uint16_t scan_ref;				/* do not penalize head gap on reference if scan_ref != 0 */
	uint16_t full_length_bonus;		/* end-to-end mapping bonus; only activated when compiled with -DDZ_FULL_LENGTH_BONUS */

	/* X-drop threshold */
	uint16_t max_gap_len;			/* as X-drop threshold */
	uint16_t unused;
} dz_score_conf_t;

/* external (custom) allocator for profile constructor */
typedef void *(*dz_malloc_t)(void *ctx, size_t size);
typedef void (*dz_free_t)(void *ctx, void *ptr);

typedef struct dz_allocator_s {
	void *ctx;
	dz_malloc_t fp;
} dz_allocator_t;

typedef struct dz_destructor_s {
	void *ctx;
	dz_free_t fp;
} dz_destructor_t;


/*
 * type declarations
 */
typedef struct dz_arena_s dz_arena_t;
typedef struct dz_profile_s dz_profile_t;
typedef struct dz_forefront_s dz_forefront_t;
typedef struct dz_state_s dz_state_t;
typedef struct dz_alignment_s dz_alignment_t;

/* pointer casting */
#define dz_state(_p)				( (dz_state_t *)(_p) )
#define dz_cstate(_p)				( (dz_state_t const *)(_p) )
#define dz_pstate(_p)				( (dz_state_t **)(_p) )
#define dz_cpstate(_p)				( (dz_state_t const **)(_p) )



#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)

/* dz_t: wrapper context */
typedef struct dz_s {
	dz_arena_t *mem;
	dz_profile_t *profile;
} dz_t;

#endif



/*
 * interface abstraction:
 * abstract query packer
 *
 * score_matrix_row == &profile->matrix[ridx * DZ_QUERY_MAT_SIZE]
 * dir == 0 for forward, 1 for reverse
 */
typedef size_t (*dz_query_calc_dim_t)(size_t qlen);						/* ref-side size */
typedef __m128i (*dz_query_conv_t)(int8_t const *score_matrix_row, uint32_t dir, __m128i v);

typedef struct {
	uint32_t dir;
	uint8_t invalid;						/* invalid base */

	/* we do not provide opaque pointer because we suppose all the information we need is contained in profile object */
	dz_query_calc_dim_t calc_dim;
	dz_query_conv_t conv;
} dz_pack_query_t;

typedef struct dz_query_s dz_query_t;		/* forward declaration */

static size_t dz_query_row_len(dz_query_t const *query);				/* query->blen * DZ_L */
static uint8_t const *dz_query_seq(dz_query_t const *query);				/* query->q */
static uint8_t const *dz_query_packed_array(dz_query_t const *query);	/* query->arr */

/*
 * core pack function forward declaration
 */
// static dz_query_t *dz_pack_query_core(dz_arena_t *mem, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack);



/*
 * interface abstraction:
 * abstract sequence fetcher
 */
typedef struct {
	uint32_t rch;
	uint32_t is_term;
} dz_fill_fetch_t;

/* exact banding for realigner */
typedef struct {
	uint32_t sidx, eidx;
} dz_fill_bound_t;

typedef struct {
	uint16_t rch;
	uint16_t hidden;
	uint32_t idx;
} dz_trace_link_t;

typedef struct {
	uint32_t match;
	int32_t score;
} dz_trace_match_t;

typedef dz_fill_fetch_t (*dz_fetch_next_t)(void *opaque, int8_t const *score_matrix, dz_query_t const *query);
typedef __m128i (*dz_get_profile_t)(void const *opaque, int8_t const *score_matrix, dz_query_t const *query, size_t qidx);
typedef dz_fill_bound_t (*dz_get_bound_t)(void const *opaque);
typedef dz_trace_match_t (*dz_get_match_t)(void const *opaque, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link);

typedef struct {
	void *opaque;
	dz_fetch_next_t const fetch_next;
	dz_get_profile_t const get_profile;
	dz_get_bound_t const get_bound;		/* can be NULL */
} dz_fetcher_t;

typedef struct {
	void const *opaque;
	dz_get_match_t const get_match;
} dz_getter_t;




/* implementations */

#if defined(DZ_NUCL_ASCII)
/*
 * matrix size is 4 x 4
 * query sequence is converted to 2bit representation in pack_query
 */
#define DZ_REF_MAT_SIZE			( 4 )
#define DZ_QUERY_MAT_SIZE		( 4 )
#define DZ_TRANSPOSE_MATRIX		( 0 )		/* do not transpose */

enum dz_alphabet {
	rA = 0x00, rC = 0x01, rG = 0x02, rT = 0x03, rU = 0x03,
	qA = 0x00, qC = 0x04, qG = 0x08, qT = 0x0c, qU = 0x0c,
	#ifdef DZ_N_AS_UNMATCHING_BASE
		rN = 0x04, qN = 0x02, qX = 0x02
	#else
		rN = 0x90, qN = 0x90, qX = 0x02		/* pshufb instruction clears the column when the 7-th bit of the index is set */
	#endif
};

static
size_t dz_query_calc_dim_ascii(size_t qlen)
{
	dz_unused(qlen);

	return(1);
}

static __dz_vectorize
__m128i dz_query_conv_ascii(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	/*
	 * ASCII to 2-bit conversion table (shifted by two bits to be used as the upper (row) shuffle index)
	 * @ABC_DEFG_HIJK_LMNO
	 * PQRS_TUVW_XYZ
	 *
	 * ['A' & 0x0f] = A
	 * ['C' & 0x0f] = C
	 * ['G' & 0x0f] = G
	 * ['T' & 0x0f] = T
	 * ['U' & 0x0f] = T
	 * ['N' & 0x0f] = N
	 */
	static uint8_t const table[16 * 2] __attribute__(( aligned(16) )) = {
		0, qA, 0, qC, qT, qU, 0, qG, 0, 0, 0, 0, 0, 0, qN, 0,
		0, qT, 0, qG, qA, qA, 0, qC, 0, 0, 0, 0, 0, 0, qN, 0
	};
	__m128i const cv = _mm_load_si128((__m128i const *)&table[16 * dir]);
	// print_vector_8(v);
	// print_vector_8(cv);
	// print_vector_8(_mm_shuffle_epi8(cv, v));
	return(_mm_shuffle_epi8(cv, v));
}


typedef struct {
	__m128i rv;
	uint8_t conv[16];
	uint8_t const *p, *t;
	int64_t inc;
} dz_fetcher_ascii_t;

static __dz_vectorize
void dz_fetcher_init_ascii(dz_fetcher_ascii_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	debug("ref(%p), rlen(%zu), p(%p), t(%p), dir(%u), inc(%d)", ref, rlen, self->p, self->t, dir, (int32_t)self->inc);

	static uint8_t const conv[16 * 2] __attribute__(( aligned(16) )) = {
		/* ['A' & 0x0f] = A, ['C' & 0x0f] = C, ['G' & 0x0f] = G, ['T' & 0x0f] = T, ['U' & 0x0f] = T, ['N' & 0x0f] = N */
		0, rA, 0, rC, rT, rU, 0, rG, 0, 0, 0, 0, 0, 0, rN, 0,
		0, rT, 0, rG, rA, rA, 0, rC, 0, 0, 0, 0, 0, 0, rN, 0
	};

	/* create conversion table (accessed by indirect %rbp) */
	__m128i const cv = _mm_load_si128((__m128i const *)&conv[16 * dir]);
	_mm_store_si128((__m128i *)self->conv, cv);
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fetch_next_ascii(dz_fetcher_ascii_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);
	dz_unused(query);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.rch     = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = self->conv[c & 0x0f];
	self->rv = _mm_set1_epi8(e);
	debug("c(%c), e(%x)", c, e);
	print_vector_8(self->rv);

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.rch     = e
	});
}

static __dz_vectorize
__m128i dz_get_profile_ascii(dz_fetcher_ascii_t const *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	uint8_t const *packed = dz_query_packed_array(query);

	__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
	__m128i const mv = _mm_load_si128((__m128i const *)score_matrix);
	__m128i const sc = _mm_shuffle_epi8(mv, _mm_or_si128(self->rv, qv));

	print_vector_8(qv);
	// print_vector_8(mv);
	// print_vector_8(sc);
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_get_match_ascii(void const *unused, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link)
{
	dz_unused(unused);

	uint32_t const rch = link->rch;
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = score_matrix[(packed[qidx] | rch) & 0x1f],
		.match = packed[qidx] == (rch<<2)
	});
}




#elif defined(DZ_NUCL_2BIT)
/*
 * 4 x 4 matrix for 2bit encoding
 */
#define DZ_REF_MAT_SIZE			( 4 )
#define DZ_QUERY_MAT_SIZE		( 4 )
#define DZ_TRANSPOSE_MATRIX		( 0 )		/* do not transpose */

enum dz_alphabet {
	rA = 0x00, rC = 0x01, rG = 0x02, rT = 0x03, rU = 0x03,
	qA = 0x00, qC = 0x04, qG = 0x08, qT = 0x0c, qU = 0x0c,
	#ifdef DZ_N_AS_UNMATCHING_BASE
		rN = 0x04, rX = 0x04, qN = 0x02, qX = 0x02
	#else
		rN = 0x90, rX = 0x04, qN = 0x90, qX = 0x02		/* pshufb instruction clears the column when the 7-th bit of the index is set */
	#endif
};

static
size_t dz_query_calc_dim_2bit(size_t qlen)
{
	dz_unused(qlen);

	return(1);
}

static __dz_vectorize
__m128i dz_query_conv_2bit(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	static uint8_t const table[16 * 2] __attribute__(( aligned(16) )) = {
		qA, qC, qG, qT, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN,		/* shift up by 2bits */
		qT, qG, qC, qA, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
	};
	__m128i const cv = _mm_load_si128((__m128i const *)&table[16 * dir]);

	return(_mm_shuffle_epi8(cv, v));
}


typedef struct {
	__m128i rv;
	uint8_t const *p, *t;
	int64_t inc;
} dz_fetcher_2bit_t;

static __dz_vectorize
void dz_fetcher_init_2bit(dz_fetcher_2bit_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fetch_next_2bit(dz_fetcher_2bit_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);
	dz_unused(query);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.rch     = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = ((self->inc>>1) ^ c) & 0x03;
	self->rv = _mm_set1_epi8(e);

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.rch     = e
	});
}

static __dz_vectorize
__m128i dz_get_profile_2bit(dz_fetcher_2bit_t const *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	uint8_t const *packed = dz_query_packed_array(query);

	__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
	__m128i const mv = _mm_load_si128((__m128i const *)score_matrix);
	__m128i const sc = _mm_shuffle_epi8(mv, _mm_or_si128(self->rv, qv));
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_get_match_2bit(void const *unused, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link)
{
	dz_unused(unused);

	uint32_t const rch = link->rch;
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = score_matrix[(packed[qidx] | rch) & 0x1f],
		.match = packed[qidx] == (rch<<2)
	});
}




#elif defined(DZ_NUCL_DIAG)
/*
 * 4 x 4 with buffer, 2bit encoding
 */
#define DZ_REF_MAT_SIZE			( 4 )
#define DZ_QUERY_MAT_SIZE		( 4 )
#define DZ_TRANSPOSE_MATRIX		( 1 )

enum dz_alphabet {
	rA = 0x00, rC = 0x01, rG = 0x02, rT = 0x03, rU = 0x03,
	qA = 0x00, qC = 0x04, qG = 0x08, qT = 0x0c, qU = 0x0c,
	#ifdef DZ_N_AS_UNMATCHING_BASE
		rN = 0x04, rX = 0x04, qN = 0x02, qX = 0x02
	#else
		rN = 0x90, rX = 0x04, qN = 0x90, qX = 0x02		/* pshufb instruction clears the column when the 7-th bit of the index is set */
	#endif
};

static
size_t dz_query_calc_dim_diag(size_t qlen)
{
	dz_unused(qlen);

	return(1);
}

static __dz_vectorize
__m128i dz_query_conv_diag(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	/* same as 2bit */
	static uint8_t const table[16 * 2] __attribute__(( aligned(16) )) = {
		qA, qC, qG, qT, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN,		/* shift up by 2bits */
		qT, qG, qC, qA, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
	};
	__m128i const cv = _mm_load_si128((__m128i const *)&table[16 * dir]);

	return(_mm_shuffle_epi8(cv, v));
}


typedef struct {
	__m128i rev, comp;		/* flip */

	uint8_t const *p, *t;
	int64_t inc;
} dz_fetcher_diag_t;

static __dz_vectorize
void dz_fetcher_init_diag(dz_fetcher_diag_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	static uint8_t const rev[32] __attribute__(( aligned(16) )) = {
		7, 6, 5, 4, 3, 2, 1, 0, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
		0, 1, 2, 3, 4, 5, 6, 7, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80
	};
	__m128i const rv = _mm_loadu_si128((__m128i const *)&rev[0x10 * dir]);

	static uint8_t const comp[32] __attribute__(( aligned(16) )) = {
		0, 0, 0, 0, 0, 0, 0, 0, rN, rN, rN, rN, rN, rN, rN, rN,
		3, 3, 3, 3, 3, 3, 3, 3, rN, rN, rN, rN, rN, rN, rN, rN
	};
	__m128i const cv = _mm_loadu_si128((__m128i const *)&comp[0x10 * dir]);

	self->rev  = rv;
	self->comp = cv;

	debug("ref(%p), rlen(%zu), dir(%u)", ref, rlen, dir);

	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fetch_next_diag(dz_fetcher_diag_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);
	dz_unused(query);

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.rch     = 0
	});
}

static __dz_vectorize
__m128i dz_get_profile_diag(dz_fetcher_diag_t const *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	/* fetch query */
	uint8_t const *packed = dz_query_packed_array(query);
	__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);

	/* fetch reference */
	size_t const fidx  = 1ULL - 8ULL - qidx;
	size_t const ridx  = qidx;
	uint8_t const *ref = &self->p[self->inc < 0 ? ridx : fidx];
	__m128i const x = _mm_loadl_epi64((__m128i const *)ref);
	debug("qidx(%zu), ref(%p, %p), inc(%zx), idx(%zx, %zx)", qidx, self->p, ref, self->inc, fidx, ridx);
	print_vector_8(x);

	/* take revcomp if needed */
	__m128i const y  = _mm_shuffle_epi8(x, self->rev);
	__m128i const rv = _mm_xor_si128(y, self->comp);

	print_vector_8(qv);
	print_vector_8(rv);

	/* shuffle */
	__m128i const mv = _mm_load_si128((__m128i const *)score_matrix);
	__m128i const sc = _mm_shuffle_epi8(mv, _mm_or_si128(rv, qv));

	print_vector_8(sc);
	return(_mm_cvtepi8_epi16(sc));
}

typedef struct {
	uint8_t const *p;
	int64_t inc;
} dz_getter_diag_t;

static __dz_vectorize
void dz_getter_init_diag(dz_getter_diag_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p   = &ref[dir ? rlen - 1ULL : 0];
	self->inc = dir ? -1LL : 1LL;
	return;
}

static __dz_vectorize
dz_trace_match_t dz_get_match_diag(dz_getter_diag_t const *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link)
{
	size_t const ridx = 1ULL + (size_t)link->idx - qidx;
	uint8_t const *ref = &self->p[self->inc * ridx];
	uint32_t const rch = ((self->inc>>1) ^ *ref) & 0x03;

	uint8_t const *packed = dz_query_packed_array(query);

	debug("idx(%u, %zu, %zu), p(%p), rch(%x), qch(%x), score(%d), match(%u)", link->idx, ridx, qidx, ref, rch, packed[qidx], score_matrix[(packed[qidx] | rch) & 0x1f], packed[qidx] == (rch<<2));
	return((dz_trace_match_t){
		.score = score_matrix[(packed[qidx] | rch) & 0x1f],
		.match = packed[qidx] == (rch<<2)
	});
}




#elif defined(DZ_NUCL_4BIT)
/*
 * matrix size is 16 x 16
 */
#define DZ_REF_MAT_SIZE			( 16 )
#define DZ_QUERY_MAT_SIZE		( 16 )
#define DZ_TRANSPOSE_MATRIX		( 1 )		/* do not transpose */

/* 4bit nucleotide (with ambiguous bases) */
enum dz_alphabet_reference {
	rA = 0x01, rC = 0x02, rG = 0x04, rT = 0x08, rU = 0x08,

	rR = rA | rG, rY = rC | rT, rS = rG | rC,
	rW = rA | rT, rK = rG | rT, rM = rA | rC,

	rB = rC | rG | rT, rD = rA | rG | rT,
	rH = rA | rC | rT, rV = rA | rC | rG,

	rN = 0x00
};

enum dz_alphabet_query {
	qA = 0x01, qC = 0x02, qG = 0x04, qT = 0x08, qU = 0x08,

	qR = qA | qG, qY = qC | qT, qS = qG | qC,
	qW = qA | qT, qK = qG | qT, qM = qA | qC,

	qB = qC | qG | qT, qD = qA | qG | qT,
	qH = qA | qC | qT, qV = qA | qC | qG,

	qN = 0x00, qX = 0x00
};

static
size_t dz_query_calc_dim_4bit(size_t qlen)
{
	dz_unused(qlen);

	return(DZ_REF_MAT_SIZE);
}

static __dz_vectorize
__m128i dz_query_conv_4bit(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	if(dir == 0) { return(v); }		/* do nothing for forward, just copy */

	static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
		0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
		0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
	};
	__m128i const cv = _mm_load_si128((__m128i const *)comp);
	__m128i const w = _mm_shuffle_epi8(cv, v);
	return(w);			/* complemented */
}


typedef struct {
	__m128i mv;
	uint8_t const *p, *t;
	int64_t inc;
	uint64_t conv;
} dz_fetcher_4bit_t;

static __dz_vectorize
void dz_fetcher_init_4bit(dz_fetcher_4bit_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	self->conv = dir ? 0xf7b3d591e6a2c480 : 0xfedcba9876543210;
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fetch_next_4bit(dz_fetcher_4bit_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(query);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.rch     = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = (self->conv>>(4 * c)) & 0x0f;
	self->mv = _mm_load_si128((__m128i const *)&score_matrix[e * DZ_QUERY_MAT_SIZE]);

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.rch     = e
	});
}

static __dz_vectorize
__m128i dz_get_profile_4bit(dz_fetcher_4bit_t const *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	dz_unused(score_matrix);

	uint8_t const *packed = dz_query_packed_array(query);

	__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
	__m128i const sc = _mm_shuffle_epi8(self->mv, qv);
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_get_match_4bit(void const *unused, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link)
{
	dz_unused(unused);

	uint32_t const rch = link->rch;
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = score_matrix[rch * DZ_QUERY_MAT_SIZE + packed[qidx]],
		.match = packed[qidx] == rch
	});
}




#elif defined(DZ_PROTEIN)
/*
 * matrix size is 32 x 32
 */
#define DZ_REF_MAT_SIZE			( 32 )
#define DZ_QUERY_MAT_SIZE		( 32 )
#define DZ_TRANSPOSE_MATRIX		( 1 )		/* do not transpose */

static
size_t dz_query_calc_dim_protein(size_t qlen)
{
	dz_unused(qlen);

	return(DZ_REF_MAT_SIZE);
}

static __dz_vectorize
__m128i dz_query_conv_protein(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(dir);

	__m128i const lcv = _mm_loadu_si128((__m128i const *)&score_matrix[0]);		/* lower conversion table */
	__m128i const hcv = _mm_loadu_si128((__m128i const *)&score_matrix[16]);	/* upper conversion table */

	/* apply lower and upper conversion mask */
	__m128i const lv = _mm_shuffle_epi8(lcv, v);
	__m128i const hv = _mm_shuffle_epi8(hcv, v);

	/* blend lower and upper with bit [4] */
	__m128i const tv = _mm_blendv_epi8(lv, hv, _mm_slli_epi32(v, 3));
	return(tv);
}


typedef struct {
	uint8_t const *p, *t;
	int64_t inc;
	uint8_t const *parr;
} dz_fetcher_protein_t;

static __dz_vectorize
void dz_fetcher_init_protein(dz_fetcher_protein_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	self->parr = NULL;
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fetch_next_protein(dz_fetcher_protein_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.rch     = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const e = *self->p & 0x1f;
	uint8_t const *packed = dz_query_packed_array(query);
	self->parr = &packed[e * dz_query_row_len(query)];

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.rch     = e
	});
}

static __dz_vectorize
__m128i dz_get_profile_protein(dz_fetcher_protein_t const *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	dz_unused(score_matrix);
	dz_unused(query);

	__m128i const sc = _mm_loadl_epi64((__m128i const *)&self->parr[qidx]);
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_get_match_protein(void const *unused, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link)
{
	dz_unused(unused);
	dz_unused(score_matrix);

	uint32_t const rch = link->rch;
	uint8_t const *seq = dz_query_seq(query);
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score =  packed[rch * dz_query_row_len(query) + qidx],
		.match = (seq[qidx - 1] & 0x1f) == rch
	});
}




#elif defined(DZ_REALIGNER)

#define DZ_SCALING_INDEX		( 6 )
#define DZ_TRANSPOSE_MATRIX		( 0 )		/* do not transpose */

/*
 * we expect the input sequences are 2bit encoded
 */
static
size_t dz_query_calc_dim_realigner(size_t qlen)
{
	dz_unused(qlen);

	return(1);
}

static __dz_vectorize
__m128i dz_query_conv_realigner(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	static uint8_t const mask[16] __attribute__(( aligned(16) )) = {
		0x00, 0x02, 0x04, 0x06, 0x08, 0x08, 0x08, 0x08,
		0x06, 0x04, 0x02, 0x00, 0x08, 0x08, 0x08, 0x08
	};

	/* take complement and double */
	__m128i const mv = _mm_load_si128(mask);
	__m128i const cv = _mm_set1_epi8(dir ? 0x08 : 0x00);
	__m128i const x  = _mm_shuffle_epi8(mv, _mm_or_si128(cv, v));
	return(x);
}


typedef struct dz_realigner_node_s {
	uint32_t a, c, g, t, gap;
	uint32_t consensus_mask;

	/* range on query */
	dz_fill_bound_t bound;

	/* links */
	struct dz_realigner_node_s *next, *prev;
} dz_realigner_node_t;

typedef struct {
	__m128i sv;
	dz_realigner_node_t const *p;
	uint32_t dir;
} dz_fetcher_realigner_t;

static __dz_vectorize
void dz_fetcher_init_realigner(dz_fetcher_realigner_t *self, uint32_t dir, dz_realigner_node_t const *ref)
{
	self->p   = ref;
	self->dir = dir;
	return;
}


static __dz_vectorize
__m128i dz_fill_calc_delta_c(dz_realigner_node_t const *node)
{
	static uint8_t const mask[16] __attribute__(( aligned(16) )) = {
		0x01, 0x01, 0x02, 0x02, 0x04, 0x04, 0x08, 0x08,
		0x10, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
	};
	__m128i const m = _mm_load_si128((__m128i const *)mask);

	__m128i const z = _mm_setzero_si128();
	__m128i const v = _mm_set1_epi8(node->consensus_mask);
	__m128i const w = _mm_cmpeq_epi8(_mm_and_si128(v, m), z);
	__m128i const x = _mm_sub_epi16(z, w);
	__m128i const y = _mm_slli_epi16(x, DZ_SCALING_INDEX);
	return(y);
}

static __dz_vectorize
__m128i dz_fill_calc_n(__m128i s, __m128i t)
{
	__m128i const x = _mm_add_epi32(s, _mm_shuffle_epi32(s, 0x4e));
	__m128i const y = _mm_add_epi32(x, _mm_shuffle_epi32(x, 0xb1));
	__m128i const z = _mm_add_epi32(y, _mm_shuffle_epi32(t, 0x00));
	return(z);
}

static __dz_vectorize
__m128i dz_fill_calc_delta_a(dz_realigner_node_t const *node)
{
	/* it's possible to use ymm here for better performance */
	__m128i const x1 = _mm_loadu_si128((__m128i const *)&node->a);		/* (t, g, c, a) */
	__m128i const x2 = _mm_loadu_si128((__m128i const *)&node->gap);	/* (-, -, -, gap) */
	__m128i const n  = dz_fill_calc_n(x1, x2);							/* a + c + g + t + gap */

	__m128i const y1 = _mm_slli_epi32(_mm_sub_epi32(n, x1), DZ_SCALING_INDEX);	/* (..., a + g + t + gap, c + g + t + gap)<<k */
	__m128i const y2 = _mm_slli_epi32(_mm_sub_epi32(n, x2), DZ_SCALING_INDEX);

	__m128 const m  = _mm_cvtepi32_ps(n);
	__m128 const z1 = _mm_div_ps(_mm_cvtepi32_ps(y1), m);
	__m128 const z2 = _mm_div_ps(_mm_cvtepi32_ps(y2), m);

	__m128i const w1 = _mm_cvtps_epi32(z1);
	__m128i const w2 = _mm_cvtps_epi32(z2);
	return(_mm_packs_epi32(w1, w2));		/* convert to 16bit */
}

static __dz_vectorize
dz_fill_fetch_t dz_fetch_next_realigner(dz_fetcher_realigner_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);

	if(self->p == NULL) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.rch     = 0
		});
	}

	/* fetch profile */
	dz_realigner_node_t const *node = self->p;

	/* compute deltas; both scaled */
	__m128i const dc = dz_fill_calc_delta_c(node);
	__m128i const da = dz_fill_calc_delta_a(node);
	self->sv = _mm_add_epi16(dc, da);

	/* forward pointer */
	dz_realigner_node_t const *link = &node->next;
	self->p = link[self->dir];
	return((dz_fill_fetch_t){
		.is_term = 0,
		.rch     = node->consensus_mask
	});
}

static __dz_vectorize
__m128i dz_get_profile_realigner(dz_fetcher_realigner_t const *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	dz_unused(score_matrix);

	uint8_t const *packed = dz_query_packed_array(query);
	__m128i const v = _mm_loadu_si128((__m128i const *)&packed[qidx]);

	/* convert to (2x + 1, 2x) pairs */
	__m128i const w = _mm_sub_epi8(v, _mm_cmpeq_epi8(self->sv, self->sv));	/* +1 */
	__m128i const x = _mm_unpacklo_epi8(v, w);

	/* shuffle 16bit array */
	__m128i const y = _mm_shuffle_epi8(self->sv, x);
	return(y);
}

static __dz_vectorize
dz_fill_bound_t dz_get_bound_realigner(dz_fetcher_realigner_t *self)
{
	dz_realigner_node_t const *node = self->p;
	return(node->bound);
}

static __dz_vectorize
dz_trace_match_t dz_get_match_realigner(void const *unused, int8_t const *score_matrix, dz_query_t const *query, size_t qidx, dz_trace_link_t const *link)
{
	dz_unused(unused);
	dz_unused(score_matrix);

	uint32_t const rch = link->rch;
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = packed[2 * qidx],
		.match = packed[2 * qidx] == 2 * rch
	});
}

#endif




/*
 * core implementations:
 * utility macros
 */
#define dz_pp_cat_intl(x, y)		x##y
#define dz_pp_cat(x, y)				dz_pp_cat_intl(x, y)
#define dz_static_assert(expr)		typedef char dz_pp_cat(_st_, __LINE__)[(expr) ? 1 : -1]
#define dz_trap()					{ *((volatile uint8_t *)NULL); }



#if defined(DZ_NUCL_ASCII)
#  define DZ_UNITTEST_INDEX			0
#elif defined(DZ_NUCL_2BIT)
#  define DZ_UNITTEST_INDEX			1
#elif defined(DZ_NUCL_DIAG)
#  define DZ_UNITTEST_INDEX			1		/* tests shared with 2bit */
#elif defined(DZ_NUCL_4BIT)
#  define DZ_UNITTEST_INDEX			2
#elif defined(DZ_PROTEIN)
#  define DZ_UNITTEST_INDEX			3
#elif defined(DZ_REALIGNER)
#  define DZ_UNITTEST_INDEX			4		/* FIXME */
#endif
#define dz_ut_sel(a, b, c, d)		( (DZ_UNITTEST_INDEX == 0) ? (a) : (DZ_UNITTEST_INDEX == 1) ? (b) : (DZ_UNITTEST_INDEX == 2) ? (c) : (d) )



/* control */
#define dz_likely(x)				__builtin_expect(!!(x), 1)
#define dz_unlikely(x)				__builtin_expect(!!(x), 0)
#define dz_roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )
#define dz_rounddown(x, base)		( (x) & ~((base) - 1) )
#define dz_max2(x, y)				( (x) < (y) ? (y) : (x) )
#define dz_min2(x, y)				( (x) < (y) ? (x) : (y) )



/* memory and pointer handling */
#define dz_add_ptr(_x, _ofs)		( (void *)((uint8_t *)(_x) + (size_t)(_ofs)) )
#define dz_sub_ptr(_x, _ofs)		( (void *)((uint8_t *)(_x) - (size_t)(_ofs)) )

#define DZ_MEM_MARGIN_SIZE			( (size_t)256 )
#define DZ_MEM_ALIGN_SIZE			( (size_t)16 )
#define DZ_MEM_INIT_SIZE			( (size_t)16 * 1024 * 1024 )

#define dz_inside(x, y, z)			( ((uint64_t)(y) - (uint64_t)(x)) < ((uint64_t)(z) - (uint64_t)(x)) )
#define dz_loadu_u64(p)				({ uint8_t const *_p = (uint8_t const *)(p); *((uint64_t const *)_p); })
#define dz_storeu_u64(p, e)			{ uint8_t *_p = (uint8_t *)(p); *((uint64_t *)(_p)) = (e); }



/* score handling */
#define dz_add_ofs(_x)				( (uint16_t)((_x) + 0x8000) )
#define dz_rm_ofs(_x)				( (int32_t)((_x) - 0x8000) )

#define DZ_CELL_MIN					( INT16_MIN )
#define DZ_CELL_MAX					( INT16_MAX )
#define DZ_L						( sizeof(__m128i) / sizeof(uint16_t) )

#define DZ_SCORE_OFS				( 64 )
#define DZ_HEAD_RCH					( 0x40 )
#define DZ_ROOT_RCH					( 0x80 )



/* vectorization supplementary */
#ifdef __SSE4_1__
#  define dz_is_all_zero(x)			( _mm_test_all_zeros((x), (x)) == 1 )
#else
// #  define dz_is_all_zero(x)			( _mm_movemask_epi8((x)) == 0 )
#endif




/*
 * memory allocation (arena):
 */
typedef struct dz_arena_block_s {
	struct dz_arena_block_s *next;
	size_t size;
} dz_arena_block_t;

typedef struct dz_stack_s {
	dz_arena_block_t *curr;
	uint8_t *top, *end;
	uint64_t _pad;
} dz_stack_t;

/* typedef */
struct dz_arena_s {
	dz_arena_block_t blk;			/* root node of linked list of blocks */
	dz_stack_t stack;				/* points a block in the linked list */
	uint64_t _pad[2];
} /* dz_arena_t */;
dz_static_assert(sizeof(dz_arena_t) == 64);
#define dz_arena_stack_rem(_mem)		( (size_t)((_mem)->stack.end - (_mem)->stack.top) )


/* aligned and margined malloc-free pair */
static __dz_force_inline
void *dz_malloc(size_t size)
{
	void *ptr = NULL;

	/* roundup to align boundary, add margin at the head and forefront */
	size = dz_roundup(size, DZ_MEM_ALIGN_SIZE);
	if(posix_memalign(&ptr, DZ_MEM_ALIGN_SIZE, size + 2 * DZ_MEM_MARGIN_SIZE) != 0) {
		debug("posix_memalign failed");
		dz_trap(); return(NULL);
	}
	debug("posix_memalign(%p), size(%lu)", ptr, size);
	return((void *)((uint8_t *)ptr + DZ_MEM_MARGIN_SIZE));
}

static __dz_force_inline
void *dz_malloc_wrap(void *unused, size_t size)
{
	dz_unused(unused);
	return(dz_malloc(size));
}

static __dz_force_inline
void dz_free(void *ptr)
{
	debug("free(%p)", (uint8_t const *)ptr - DZ_MEM_MARGIN_SIZE);
	free((void *)((uint8_t *)ptr - DZ_MEM_MARGIN_SIZE));
	return;
}

static __dz_force_inline
void dz_free_wrap(void *unused, void *ptr)
{
	dz_unused(unused);
	return(dz_free(ptr));
}

unittest() {
	size_t size[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
	for(size_t i = 0; i < sizeof(size) / sizeof(size_t); i++) {
		uint8_t *p = (uint8_t *)dz_malloc(size[i]);
		ut_assert(p != NULL);
		volatile uint8_t a = p[-32], b = p[0], c = p[size[i]], d = p[size[i] + 32];
		dz_unused(a); dz_unused(b); dz_unused(c); dz_unused(d);
		memset(p, 0, size[i]);
		dz_free(p);
	}
}

/* chained stack */
static __dz_force_inline
void dz_arena_flush(dz_arena_t *mem)
{
	mem->stack.curr = &mem->blk;
	mem->stack.top = (uint8_t *)dz_roundup((uintptr_t)(mem + 1), DZ_MEM_ALIGN_SIZE);	/* space for stack always reserved */
	mem->stack.end = (uint8_t *)mem + mem->blk.size - DZ_MEM_MARGIN_SIZE;
	return;
}

static __dz_force_inline
void dz_arena_init(dz_arena_t *mem, size_t size)
{
	/* init mem object then stack pointers */
	mem->blk = (dz_arena_block_t){
		.next = NULL,
		.size = size
	};

	dz_arena_flush(mem);
	return;
}

static __dz_force_inline
dz_arena_t *dz_arena_create(size_t size)
{
	size = dz_min2(DZ_MEM_INIT_SIZE, sizeof(dz_arena_t) + dz_roundup(size, DZ_MEM_ALIGN_SIZE));
	dz_arena_t *mem = (dz_arena_t *)dz_malloc(size);
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* init mem object then stack pointers */
	dz_arena_init(mem, size);
	return(mem);
}

static __dz_force_inline
void dz_arena_destroy(dz_arena_t *mem)
{
	dz_arena_block_t *blk = mem->blk.next;
	debug("cleanup memory chain, blk(%p)", blk);
	while(blk != NULL) {
		dz_arena_block_t *next = blk->next;
		debug("free blk(%p), next(%p)", blk, next);
		dz_free(blk); blk = next;
	}
	dz_free(mem);
	return;
}
static __dz_force_inline
uint64_t dz_arena_push_stack(dz_arena_t *mem, size_t size)
{
	debug("add_stack, ptr(%p)", mem->stack.curr->next);
	if(dz_unlikely(mem->stack.curr->next == NULL)) {
		/* current stack is the forefront of the memory block chain, add new block */
		size = dz_max2(
			size + dz_roundup(sizeof(dz_arena_block_t), DZ_MEM_ALIGN_SIZE),
			2 * mem->stack.curr->size
		);
		dz_arena_block_t *blk = (dz_arena_block_t *)dz_malloc(size);
		debug("malloc called, blk(%p)", blk);
		if(blk == NULL) { return(1); }

		/* link new node to the forefront of the current chain */
		mem->stack.curr->next = blk;

		/* initialize the memory block */
		blk->next = NULL;
		blk->size = size;
	}

	/* follow a forward link and init stack pointers */
	mem->stack.curr = mem->stack.curr->next;
	mem->stack.top = (uint8_t *)dz_roundup((uintptr_t)(mem->stack.curr + 1), DZ_MEM_ALIGN_SIZE);
	mem->stack.end = (uint8_t *)mem->stack.curr + mem->stack.curr->size;		/* margin already added at the tail by dz_malloc */
	return(0);
}
static __dz_force_inline
void *dz_arena_malloc(dz_arena_t *mem, size_t size)
{
	size = dz_roundup(size, sizeof(__m128i));
	debug("size(%zu), stack(%p, %p, %zu)", size, mem->stack.top, mem->stack.end, dz_arena_stack_rem(mem));

	if(dz_unlikely(size > dz_arena_stack_rem(mem) || dz_arena_stack_rem(mem) < 4096)) {
		dz_arena_push_stack(mem, size);
	}
	void *ptr = (void *)mem->stack.top;
	mem->stack.top += size;

	debug("size(%zu), stack(%p, %p, %zu)", size, mem->stack.top, mem->stack.end, dz_arena_stack_rem(mem));
	return(ptr);
}
static __dz_force_inline
uint64_t dz_arena_pop_stack(dz_arena_t *mem)
{
	dz_arena_block_t *p = &mem->blk, *q = NULL;
	while(p->next != NULL) {
		q = p;
		p = p->next;
	}
	if(q == NULL) { return(0); }

	/* free the last block and fixup link */
	dz_free(p);
	q->next = NULL;
	return(1);
}




/*
 * save and restore malloc state; NOTE: unstable
 */
typedef struct dz_freeze_s {
	dz_stack_t stack;
} dz_freeze_t;
dz_static_assert(sizeof(dz_freeze_t) == 32);

static __dz_force_inline
dz_freeze_t *dz_arena_freeze(dz_arena_t *mem)
{
	/*
	 * be careful! freezed stack is invalid after restored once.
	 */
	dz_stack_t const stack = mem->stack;

	dz_freeze_t *fz = dz_arena_malloc(mem, sizeof(dz_freeze_t));
	fz->stack = stack;
	return(fz);
}
static __dz_force_inline
void dz_arena_restore(dz_arena_t *mem, dz_freeze_t const *fz)
{
	mem->stack = fz->stack;
	return;
}


unittest() {
	size_t size[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
	for(size_t i = 0; i < sizeof(size) / sizeof(size_t); i++) {
		dz_arena_t *mem = dz_arena_create(size[i]);
		uint8_t *p = (uint8_t *)dz_arena_malloc(mem, size[i]);
		ut_assert(p != NULL);
		memset(p, 0, size[i]);
		dz_arena_destroy(mem);
	}
}

static __dz_force_inline
uint8_t *dz_reserve_stack(dz_arena_t *mem, size_t size)
{
	/* allocate from heap */
	if(dz_unlikely(dz_arena_stack_rem(mem) < size)) {
		dz_arena_push_stack(mem, 0);
	}
	return(mem->stack.top);
}

static __dz_force_inline
void dz_save_stack(dz_arena_t *mem, void *ptr)
{
	mem->stack.top = (uint8_t *)ptr;
	return;
}




/*
 * score profile object
 */
/* typedef */
struct dz_profile_s {
	uint16_t iiv[8], iev[8];
	uint16_t div[8], dev[8];	/* deletion: open and extend */

	/* root column */
	struct dz_forefront_s const *root;

	/* constants and miscellaneous */
	uint32_t size;				/* object size for copying */
	uint16_t max_gap_len;
	uint16_t bonus;

	/* constants for fill-in */
	uint32_t ofs;				/* dz_add_ofs(0); never zero */
	uint32_t init;				/* dz_add_ofs(init) */
	uint32_t xth;				/* xth; never zero */
	uint32_t one;				/* 1; never zero */

	int8_t matrix[];
} /* dz_profile_t */;
dz_static_assert(sizeof(dz_profile_t) % sizeof(__m128i) == 0);


static __dz_vectorize
__m128i dz_profile_is_scan(dz_profile_t const *profile)
{
	/* the same constant used in dz_state_merge_core */
	__m128i const negv = _mm_cvtsi64_si128(INT32_MAX);	/* INT32_MAX, 0, 0, 0 */

	__m128i const x = _mm_load_si128((__m128i const *)&profile->ofs);
	__m128i const y = _mm_cmpeq_epi32(x, negv);			/* (0, is_scan, 0, 0) */
	__m128i const z = _mm_and_si128(x, negv);			/* (dz_add_ofs(0), 0, 0, 0) */
	__m128i const w = _mm_add_epi32(y, z);

	/* assertions */
	debugblock({
		assert(_mm_extract_epi32(w, 0) == dz_add_ofs(0));
		assert(_mm_extract_epi32(w, 1) == (profile->init != (uint16_t)dz_add_ofs(DZ_CELL_MIN) ? 0 : -1));
		assert(_mm_extract_epi32(w, 2) == 0);
		assert(_mm_extract_epi32(w, 3) == 0);
	});

	/* [1] is -1 if not scan, otherwise 0; [0] is always dz_add_ofs(0) */
	return(w);
}




/*
 * query packing:
 *
 * blen = roundup(qlen, DZ_L) / DZ_L
 * array must have 16-byte-length margin at the tail
 */
/* typedef */
struct dz_query_s {
	uint32_t bmax;		/* always INT32_MAX */
	uint32_t blen;		/* #vectors */
	uint8_t const *seq;	/* query string (encoded and aligned) */
	int16_t bonus[2 * DZ_L];
	uint8_t arr[];
} /* dz_query_t */;
dz_static_assert(sizeof(dz_query_t) % sizeof(__m128i) == 0);

/* full-length bonus */
static __dz_vectorize
__m128i dz_query_get_bonus(dz_query_t const *query, uint64_t is_tail)
{
	#ifdef DZ_FULL_LENGTH_BONUS
		__m128i const *pbonus = &((__m128i const *)query->arr)[-2];
		__m128i const bv = _mm_load_si128(&pbonus[is_tail]);

	#else
		dz_unused(query);
		dz_unused(is_tail);

		__m128i const bv = _mm_setzero_si128();

	#endif
	return(bv);
}

static __dz_vectorize
size_t dz_query_row_len(dz_query_t const *query)
{
	return(query->blen * DZ_L);
}

static __dz_vectorize
uint8_t const *dz_query_seq(dz_query_t const *query)
{
	return(query->seq);
}

static __dz_vectorize
uint8_t const *dz_query_packed_array(dz_query_t const *query)
{
	return(query->arr);
}


#if !defined(DZ_NUCL_DIAG)
static __dz_vectorize
size_t dz_pack_query_clen(dz_profile_t const *profile, size_t qlen)
{
	dz_unused(profile);

	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	return(clen);
}

#else
static __dz_vectorize
size_t dz_pack_query_clen(dz_profile_t const *profile, size_t qlen)
{
	size_t const clen = dz_roundup(qlen + 1 + profile->max_gap_len, DZ_L);
	return(clen);
}

#endif


/* profile vector is calculated on-the-fly */
static __dz_vectorize
dz_query_t *dz_pack_query_alloc_mem(dz_arena_t *mem, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	size_t const clen = dz_pack_query_clen(profile, qlen);
	size_t const dim  = pack->calc_dim(clen);

	/* allocate mem */
	dz_query_t *q = (dz_query_t *)dz_arena_malloc(mem,
		  sizeof(dz_query_t)	/* header */
		+ clen * dim			/* packed sequence (body) */
		+ sizeof(__m128i)		/* margin at the tail */
	);

	/* save length and seq pointers */
	*q = (dz_query_t){
		.bmax  = INT32_MAX,
		.blen  = qlen == 0 ? 0 : (clen / DZ_L),
		.seq   = query,
		.bonus = { 0 }
	};

	/*
	 * tentative support of full-length bonus;
	 * precaclulated bonus score is saved at q->bonus[0] and q->bonus[1];
	 * [0] for non-tail vector and [1] for tail vector
	 */
	q->bonus[DZ_L + (qlen % DZ_L)] = profile->bonus;
	return(q);
}

static __dz_vectorize
size_t dz_pack_query_forward_core(dz_query_t *q, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	/*
	 * determine lengths
	 *                         clen               blen
	 * qlen == k * DZ_L        (k + 1) * DZ_L      k      * DZ_L
	 * qlen == k * DZ_L - 1     k      * DZ_L     (k - 1) * DZ_L
	 * qlen == k * DZ_L - l     k      * DZ_L     (k - 1) * DZ_L
	 */
	size_t const clen = dz_pack_query_clen(profile, qlen);
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));
	size_t const dim = pack->calc_dim(clen);
	debug("clen(%zu), blen(%zu), dim(%zu)", clen, blen, dim);

	/* calculate profile for each possible reference-side base (residue) */
	for(size_t j = 0; j < dim; j++) {
		debug("j(%zu)", j);

		/* extract array pointers for this row */
		int8_t const *score_matrix = &profile->matrix[j * DZ_QUERY_MAT_SIZE];
		uint8_t *a = &q->arr[j * clen];

		/* until the end of the query sequence */
		__m128i pv = _mm_setzero_si128();	/* as -DZ_SCORE_OFS */
		for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
			__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
			// print_vector_8(qv);
			__m128i const tv = pack->conv(score_matrix, 0, qv);

			/* shift by one to make room for the top row */
			_mm_storeu_si128((__m128i *)&a[i], _mm_alignr_epi8(tv, pv, 15));

			/* save for the next iteration */
			pv = tv;
		}

		/* continue the same conversion on the remainings */ {
			char qbuf[16] = { 0 };
			for(size_t i = 0; i < qlen - blen; i++) {	/* qlen - blen < 16 */
				qbuf[i] = query[i + blen];
			}

			__m128i const qv = _mm_loadu_si128((__m128i const *)&qbuf[0]);
			__m128i const tv = pack->conv(score_matrix, 0, qv);
			_mm_storeu_si128((__m128i *)&a[blen], _mm_alignr_epi8(tv, pv, 15));
		}

		/* clear tail remainders */
		for(size_t i = qlen + 1; i < clen; i++) {
			a[i] = pack->invalid;		/* as -DZ_SCORE_OFS */
		}
	}
	// debug("qlen(%lu), q(%s)", qlen, query);
	return(clen);
}

static __dz_vectorize
size_t dz_pack_query_reverse_core(dz_query_t *q, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	/* constants */
	size_t const clen = dz_pack_query_clen(profile, qlen);
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));
	size_t const dim = pack->calc_dim(clen);
	debug("clen(%zu), blen(%zu), dim(%zu)", clen, blen, dim);

	static uint8_t const rev[16] __attribute__(( aligned(16) )) = {
		15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
	};
	__m128i const rv = _mm_load_si128((__m128i const *)rev);


	/* calculate profile for each possible reference-side base (residue) */
	for(size_t j = 0; j < dim; j++) {
		debug("j(%zu)", j);

		/* extract array pointers for this row */
		int8_t const *score_matrix = &profile->matrix[j * DZ_QUERY_MAT_SIZE];
		uint8_t *a = &q->arr[j * clen];

		/* toward the head of the query sequence */
		__m128i pv = _mm_setzero_si128();
		for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
			__m128i const qv = _mm_loadu_si128((__m128i const *)&query[qlen - 16 - i]);
			__m128i const tv = pack->conv(score_matrix, 1, _mm_shuffle_epi8(qv, rv));		/* bswapdq before apply conversion */

			/* shift by one to make room for the top row */
			_mm_storeu_si128((__m128i *)&a[i], _mm_alignr_epi8(tv, pv, 15));

			/* save for the next iteration */
			pv = tv;
		}

		/* continue the same conversion on the remainings */ {
			char qbuf[16] = { 0 };
			for(size_t i = 0; i < qlen - blen; i++) {
				qbuf[i] = query[qlen - blen - 1 - i];
				debug("i(%zu), c(%c, %u), c(%c, %u)", i, query[qlen - blen - 1 - i], query[qlen - blen - 1 - i], qbuf[i], qbuf[i]);
			}

			__m128i const qv = _mm_loadu_si128((__m128i const *)&qbuf[0]);
			__m128i const tv = pack->conv(score_matrix, 1, qv);		/* already reversed */
			_mm_storeu_si128((__m128i *)&a[blen], _mm_alignr_epi8(tv, pv, 15));
		}

		/* clear tail remainders */
		for(size_t i = qlen + 1; i < clen; i++) {
			a[i] = pack->invalid;		/* as -DZ_SCORE_OFS */
		}
	}
	// debug("qlen(%lu), q(%s)", qlen, query);
	return(clen);
}

/* direction-folded API */
static __dz_vectorize
dz_query_t *dz_pack_query_core(dz_arena_t *mem, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	debugblock({
		assert(qlen / DZ_L < INT32_MAX, "qlen(%zu)", qlen);
	});
	if(qlen >= (size_t)DZ_L * (size_t)INT32_MAX) { return(NULL); }


	dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, query, qlen, pack);
	if(pack->dir) {
		dz_pack_query_reverse_core(q, profile, query, qlen, pack);
	} else {
		dz_pack_query_forward_core(q, profile, query, qlen, pack);
	}
	return(q);
}




/*
 * DP matrix internal objects
 */

/* Smith-Waterman-Gotoh layer indices */
enum dz_layer {
	DZ_E_MATRIX = 0,
	DZ_F_MATRIX = 1,
	DZ_S_MATRIX = 2
};

/* score matrix block; followed by dz_cap_s */
typedef struct dz_swgv_s {
	__m128i e, f, s;	/* horizontal gap, vertical gap, and score */
} dz_swgv_t;
typedef struct {
	uint16_t e[DZ_L], f[DZ_L], s[DZ_L];
} dz_swgv_raw_t;
dz_static_assert(sizeof(dz_swgv_t) % sizeof(__m128i) == 0);
#define dz_swgv(_p)					( (dz_swgv_t *)(_p) )
#define dz_cswgv(_p)				( (dz_swgv_t const *)(_p) )


static __dz_vectorize
dz_swgv_t dz_load_swgv(dz_swgv_t const *p)
{
	dz_swgv_t v;
	v.e = _mm_loadu_si128(&p->e);
	v.f = _mm_loadu_si128(&p->f);
	v.s = _mm_loadu_si128(&p->s);
	return(v);
}

static __dz_vectorize
dz_swgv_t *dz_store_swgv(dz_swgv_t *p, dz_swgv_t v)
{
	_mm_storeu_si128(&p->e, v.e);
	_mm_storeu_si128(&p->f, v.f);
	_mm_storeu_si128(&p->s, v.s);
	return(p);
}

static __dz_vectorize
dz_swgv_t dz_subs_swgv(dz_swgv_t v, __m128i ofs)
{
	v.e = _mm_subs_epu16(v.e, ofs);
	v.f = _mm_subs_epu16(v.f, ofs);
	v.s = _mm_subs_epu16(v.s, ofs);
	return(v);
}

static __dz_vectorize
dz_swgv_t dz_max_swgv(dz_swgv_t v, dz_swgv_t w)
{
	return((dz_swgv_t){
		.e = _mm_max_epu16(v.e, w.e),
		.f = _mm_max_epu16(v.f, w.f),
		.s = _mm_max_epu16(v.s, w.s)
	});
}




/*
 * DP matrix outline info: consist of column range, column count, and score tracker
 * exposed to API layer so that users can control DP extension seeing scores and ranges.
 * merged object created at the beginning of dz_extend and converted to internal representations.
 */

/* score tracker */
typedef struct {
	uint32_t inc;		/* with offset */
	int32_t abs;		/* without offset */
} dz_score_t;

static __dz_vectorize
int32_t dz_score_init(dz_score_t *score)
{
	/*
	 * score->abs = 0;
	 * score->inc = dz_add_ofs(0);
	 */

	dz_storeu_u64(score, dz_add_ofs(0));
	return(0);
}

static __dz_vectorize
int32_t dz_score_abs(dz_score_t const *score)
{
	return(score->abs);
}

static __dz_vectorize
int32_t dz_score_base(dz_score_t const *score)
{
	return(score->abs - dz_rm_ofs(score->inc));
}

static __dz_vectorize
int32_t dz_score_flush(dz_score_t *score)
{
	uint32_t const inc = score->inc;
	score->abs += dz_rm_ofs(inc);		/* with offset; FIXME: fold tracker->adj ?? */
	return(score->abs);
}

static __dz_vectorize
uint32_t dz_score_rebase(dz_score_t *score, uint32_t delta)
{
	#if 0
	uint32_t const inc = score->inc - delta;
	score->abs += dz_rm_ofs(inc);
	score->inc = dz_add_ofs(0);
	debugblock({
		// assert(dz_rm_ofs(inc) >= 0, "inc(%d)", dz_rm_ofs(inc));
		assert(dz_rm_ofs(inc) <= UINT16_MAX, "inc(%d)", dz_rm_ofs(inc));
	});
	return((uint32_t)dz_rm_ofs(inc));	/* always positive */
	#endif

	score->abs += delta;
	score->inc -= delta;
	return((uint32_t)dz_rm_ofs(score->inc));
}




/*
 * coordinate and match information for traceback
 */
typedef struct {
	uint16_t rch;		/* fetched char (after conversion) */
	uint16_t adj;		/* difference of offset to the previous column */
	uint32_t idx;		/* forwarded length */
} dz_link_t;
dz_static_assert(sizeof(dz_link_t) == 8);
dz_static_assert(offsetof(dz_link_t, rch) == offsetof(dz_trace_link_t, rch));
dz_static_assert(offsetof(dz_link_t, idx) == offsetof(dz_trace_link_t, idx));
#define dz_ctlink(x)		( (dz_trace_link_t const *)(x) )

static __dz_vectorize
void dz_link_init_root(dz_link_t *link)
{
	/*
	 * link->rch = DZ_HEAD_RCH | DZ_ROOT_RCH;
	 * link->adj = 0;
	 * link->idx = 0;	// no incoming vector for the root
	 */

	dz_storeu_u64(link, DZ_HEAD_RCH | DZ_ROOT_RCH);
	return;
}

static __dz_vectorize
void dz_link_init_intl(dz_link_t *link)
{
	/*
	 * // still regarded as head
	 * link->rch = DZ_HEAD_RCH | 0x01;		// to make the value different from head
	 * link->adj = 0;
	 * link->idx = 1;		// #incoming vectors == 1
	 */

	dz_storeu_u64(link, DZ_HEAD_RCH | 0x01 | (0x01ULL<<32));
	return;
}

static __dz_vectorize
void dz_link_init_head(dz_link_t *link, size_t fcnt)
{
	/*
	 * link->rch = DZ_HEAD_RCH;
	 * link->adj = 0;
	 * link->idx = fcnt;	// save incoming vector count
	 */

	dz_storeu_u64(link, DZ_HEAD_RCH | ((uint64_t)fcnt<<32));
	return;
}

static __dz_vectorize
size_t dz_link_fcnt(dz_link_t const *link)
{
	return(link->idx);
}

static __dz_vectorize
uint64_t dz_link_is_root(dz_link_t const *link)
{
	/* test for root tail */
	return((link->rch & DZ_ROOT_RCH) != 0);
}

static __dz_vectorize
uint64_t dz_link_is_intl(dz_link_t const *link)
{
	return(link->rch != DZ_HEAD_RCH);
}

static __dz_vectorize
uint64_t dz_link_is_head(dz_link_t const *link)
{
	/* root or head */
	return((link->rch & DZ_HEAD_RCH) != 0);
}




/*
 * column range;
 * placed just after every score vector to indicate the length
 */
typedef struct dz_range_s {
	uint32_t sblk, eblk;	/* uint32_t pair so that we can take max of ranges with pmaxud */
} dz_range_t;
#define dz_range(_p)				( (dz_range_t *)(_p) )
#define dz_crange(_p)				( (dz_range_t const *)(_p) )

static __dz_vectorize
dz_range_t dz_range_estimate_max(size_t max_gap_len)
{
	/* roundup */
	max_gap_len = dz_roundup(max_gap_len, DZ_L);
	return((dz_range_t){
		.sblk = 0,
		.eblk = max_gap_len / DZ_L + 1
	});
}

/*
 * range manipulation functions
 */
static __dz_vectorize
uint64_t dz_range_is_head(dz_range_t const *r)
{
	return(r->sblk == 0);
}

static __dz_vectorize
uint64_t dz_range_is_term(dz_range_t const *r)
{
	debug("test range(%u, %u)", r->sblk, r->eblk);
	return(r->sblk >= r->eblk);
}

static __dz_vectorize
uint64_t dz_range_is_bottom(dz_range_t const *r, uint32_t blen)
{
	debug("test bottom(%u, %u)", r->eblk, blen);
	return(r->eblk >= blen);
}

static __dz_vectorize
uint64_t dz_range_shrink(dz_range_t *r)
{
	r->sblk++;
	return(dz_range_is_term(r));
}

static __dz_vectorize
uint64_t dz_range_extend(dz_range_t *r, uint32_t blen)
{
	r->eblk++;
	return(dz_range_is_bottom(r, blen));
}




/*
 * max tracker
 */
typedef struct {
	dz_score_t score;					/* score.abs is maximum of incoming scores (or zero in scan), score.inc be dz_add_ofs(0) */
	struct dz_cap_s const *cap;			/* NULL on init */
} dz_max_t;

#if 0
static __dz_vectorize
uint64_t dz_max_is_overflow(dz_max_t const *max)
{
	debug("inc(%u), thresh(%u)", max->score.inc, UINT16_MAX - DZ_OVERFLOW_MARGIN);
	return(max->score.inc > (UINT16_MAX - DZ_OVERFLOW_MARGIN));
}
#endif

#if defined(DZ_OFFSET)
static __dz_vectorize
void dz_max_forward(dz_max_t *max, size_t fwd)
{
	dz_unused(max);
	dz_unused(fwd);

	max->score.abs -= 8 * fwd;
	max->score.inc += 8 * fwd;
	return;
}

#else		/* DZ_OFFSET */
static __dz_vectorize
void dz_max_forward(dz_max_t *max, size_t fwd)
{
	dz_unused(max);
	dz_unused(fwd);

	return;
}

#endif		/* DZ_OFFSET */

static __dz_vectorize
uint32_t dz_max_rebase(dz_max_t *max, uint32_t inc)
{
	return(dz_score_rebase(&max->score, inc));
}

static __dz_vectorize
uint32_t dz_max_update(dz_max_t *max, struct dz_cap_s const *cap, uint32_t xth, uint32_t diff_inc)
{
	/* update max; the actual score (absolute score accumulated from the origin) is expressed as max + inc; 2 x cmov */
	debug("max(%u), (%u, %u)", diff_inc, max->score.inc, max->score.abs);

	/* with offset */
	if(!dz_cmp_max(diff_inc, xth)) {
		return(0);
	}

	uint32_t const next_inc = diff_inc + max->score.inc - xth;
	debug("max updated, inc(%u, %d, %u)", next_inc, dz_rm_ofs(next_inc), diff_inc);

	max->score.inc = next_inc;
	max->cap = cap;
	return(dz_rm_ofs(next_inc));
}




/*
 * DP matrix scaffolds
 */

/* placed at the head */
typedef struct dz_head_s {
	dz_score_t base;	/* base score; always positive */
	dz_link_t link;
} dz_head_t;
#define dz_head(_p)					( (dz_head_t *)(_p) )
#define dz_chead(_p)				( (dz_head_t const *)(_p) )

/* followed by dz_tail_t; sblk and eblk are shared to forefront_s */
typedef struct dz_cap_s {
	dz_range_t range;	/* column range */
	dz_link_t link;		/* rch, adj, idx */
} dz_cap_t;
dz_static_assert(sizeof(dz_cap_t) % sizeof(__m128i) == 0);
#define dz_cap(_p)					( (dz_cap_t *)(_p) )
#define dz_ccap(_p)					( (dz_cap_t const *)(_p) )

/* dz_head_t and dz_cap_t are compatibile for head two elements */
dz_static_assert(offsetof(dz_head_t, link) == offsetof(dz_cap_t, link));
dz_static_assert(offsetof(dz_head_t, base) == offsetof(dz_cap_t, range));



#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
#  if defined(DZ_BENCH_CNT)
#    define DZ_TAIL_MARGIN				( 32 + 32 + 32 )
#  else
#    define DZ_TAIL_MARGIN				( 32 + 32 )		/* dz_meta_t */
#  endif
#else
#  define DZ_TAIL_MARGIN				( 32 )
#endif



/* scaffold handling */
static __dz_force_inline
size_t dz_calc_column_size(dz_range_t const *range)
{
	size_t const span = range->eblk - range->sblk;
	size_t const column_size = sizeof(dz_swgv_t) * (span + 8);
	size_t const size = (
		  sizeof(dz_cap_t)	/* header size */
		+ column_size		/* estimated column size */
		+ DZ_TAIL_MARGIN	/* tail cap and metadata */
	);
	return(size);
}

static __dz_vectorize
void dz_save_incoming(dz_head_t *head, dz_state_t const **ff, size_t fcnt)
{
	dz_state_t const **q = dz_cpstate(head) - fcnt;
	for(size_t i = 0; i < fcnt; i++) {
		q[i] = ff[i];
		debug("i(%zu), fcnt(%zu), ff(%p)", i, fcnt, ff[i]);
	}
	return;
}

static __dz_vectorize
dz_head_t *dz_slice_head(dz_arena_t *mem, dz_range_t const *range, dz_state_t const **ff, size_t fcnt)
{
	dz_unused(ff);

	/* forefront pointer array size */
	size_t const fsize = dz_roundup(
		sizeof(dz_state_t *) * fcnt,
		sizeof(__m128i)
	);

	/* column size (including head and internal caps) */
	size_t const csize = dz_calc_column_size(range);

	/* allocate mem */
	dz_head_t *head = dz_head(dz_add_ptr(dz_reserve_stack(mem, fsize + csize), fsize));
	debug("csize(%zu), fsize(%zu), stack(%p, %p, %zu), head(%p)", csize, fsize, mem->stack.top, mem->stack.end, dz_arena_stack_rem(mem), head);
	return(head);
}

static __dz_vectorize
dz_cap_t *dz_init_head(dz_head_t *head, dz_max_t const *max, dz_state_t const **ff, size_t fcnt)
{
	dz_unused(ff);

	/* save current max as base */
	uint64_t const x = dz_loadu_u64(&max->score);
	dz_storeu_u64(&head->base, x);		/* inc contains offset */

	/* mark head */
	dz_link_init_head(&head->link, fcnt);
	return(dz_cap(head));
}

static __dz_vectorize
dz_cap_t *dz_init_root_head(dz_head_t *head)
{
	/* clear score */
	dz_score_init(&head->base);

	/* mark root */
	dz_link_init_root(&head->link);
	return(dz_cap(head));
}

static __dz_vectorize
uint64_t dz_is_head(dz_cap_t const *cap)
{
	return(dz_link_is_head(&cap->link));
}

static __dz_vectorize
uint64_t dz_is_root(dz_cap_t const *cap)
{
	return(dz_link_is_root(&cap->link));
}




/* ping-pong pointer conversion between dz_cap_t and dz_swgv_t */
static __dz_vectorize
dz_swgv_t *dz_anchor_col(dz_swgv_t *col, size_t sblk)
{
	return(col - sblk);
}

static __dz_vectorize
dz_swgv_t *dz_slice_col(dz_cap_t *prev_cap, size_t sblk)
{
	dz_swgv_t *col = dz_swgv(prev_cap + 1);
	return(dz_anchor_col(col, sblk));
}

static __dz_vectorize
dz_swgv_t const *dz_restore_col(dz_cap_t const *cap)
{
	return(dz_cswgv(cap) - cap->range.eblk);
}

static __dz_vectorize
dz_cap_t *dz_slice_cap(dz_swgv_t *col, size_t eblk)
{
	/* put cap at the tail of current column */
	return(dz_cap(&col[eblk]));
}

static __dz_vectorize
dz_cap_t const *dz_unwind_cap(dz_cap_t const *cap)
{
	dz_swgv_t const *col = dz_restore_col(cap);
	debug("cap(%p), range(%u, %u), col(%p), prev_cap(%p)", cap, cap->range.sblk, cap->range.eblk, col, dz_ccap(col + cap->range.sblk) - 1);
	return(dz_ccap(col + cap->range.sblk) - 1);
}

static __dz_vectorize
uint16_t *dz_init_cap(dz_cap_t *cap, dz_range_t const *range, dz_link_t const *link)
{
	/* link can be NULL (for root) */
	// dz_link_t const lk = { .rch = 0, .idx = 0 };
	// if(link == NULL) { link = &lk; }

	/*
	 * cap->range = *range;
	 * cap->link  = *link;
	 */
	uint64_t const x = dz_loadu_u64(range);
	uint64_t const y = dz_loadu_u64(link);

	dz_storeu_u64(&cap->range, x);
	dz_storeu_u64(&cap->link,  y);

	/* returns pointer to adj */
	return(&cap->link.adj);
}




/* internal link (ilink) */
static __dz_vectorize
dz_swgv_t *dz_slice_ilink(uint8_t *ptr, dz_cap_t const *prev_cap)
{
	/* save link */
	dz_state_t const **ff = ((dz_state_t const **)ptr) + 2;
	ff[-1] = dz_cstate(&prev_cap->range);

	/* slice head */
	dz_head_t *head = dz_head(ff);

	/* score and link */
	dz_score_init(&head->base);
	dz_link_init_intl(&head->link);		/* #incoming vectors == 1 */

	/* slice column */
	return(dz_swgv(head + 1));
}

static __dz_vectorize
uint64_t dz_is_ilink(dz_cap_t const *cap)
{
	return(dz_link_is_intl(&cap->link));
}

static __dz_vectorize
dz_cap_t const *dz_rewind_ilink(dz_cap_t const *cap)
{
	dz_state_t const **ff = dz_cpstate(cap);
	return(dz_ccap(ff[-1]));
}




/*
 * column tracker
 * dz_ref_cnt_t is only used in tail caps, which is converted from dz_link_t
 */
typedef struct {
	uint32_t section;		/* section count */
	uint32_t column;		/* column count */
} dz_ref_cnt_t;

/*
 * dz_state_t is internally used as DP matrix object.
 */
struct dz_state_s {
	dz_range_t range;
	dz_ref_cnt_t cnt;
	dz_max_t max;

	#if defined(DZ_BENCH_CNT)
		size_t vcnt, unused[3];
	#endif
} /* dz_state_t */;

/* dz_state_t and dz_cap_t are compatibile */
dz_static_assert(offsetof(dz_state_t, range)      == offsetof(dz_cap_t, range));
dz_static_assert(offsetof(dz_state_t, cnt.column) == offsetof(dz_cap_t, link.idx));


/*
 * to take min between eblk and blen
 */
dz_static_assert(offsetof(dz_query_t, blen) == offsetof(dz_state_t, range.eblk));


static __dz_vectorize
void dz_init_tail(dz_state_t *tail, dz_range_t const *range)
{
	/*
	 * initialize state:
	 * tail->cnt.section = 0;
	 * tail->cnt.column  = 0;
	 * tail->range.sblk  = 0;
	 * tail->range.eblk  = eblk;
	 */
	uint64_t const x = dz_loadu_u64(range);
	dz_storeu_u64(&tail->cnt,   0);
	dz_storeu_u64(&tail->range, x);

	/* max tracker */
	dz_score_init(&tail->max.score);
	tail->max.cap = NULL;
	return;
}

static __dz_vectorize
dz_state_t *dz_slice_tail(dz_arena_t *mem, dz_swgv_t *col, size_t eblk)
{
	dz_state_t *tail = dz_state(dz_slice_cap(col, eblk));

	/* write back stack pointer */
	size_t const size = dz_roundup(
		DZ_TAIL_MARGIN,		/* make room for tail and metadata */
		sizeof(__m128i)
	);
	dz_save_stack(mem, dz_add_ptr(tail, size));
	return(tail);
}

static __dz_vectorize
dz_swgv_t const *dz_restore_tail_column(dz_state_t const *tail)
{
	debug("tail(%p), range(%u, %u), col(%p)", tail, tail->range.sblk, tail->range.eblk, dz_cswgv(tail) - tail->range.eblk);
	return(dz_cswgv(tail) - tail->range.eblk);
}




/*
 * dz_forefront_t is alias of dz_state_t, only for backward compatibility.
 */
struct dz_forefront_s {
	dz_range_t range;
	dz_ref_cnt_t cnt;

	/* score exposed outside */
	uint32_t inc;
	int32_t max;

	dz_cap_t const *cap;

	#if defined(DZ_BENCH_CNT)
		size_t vcnt, unused[3];
	#endif
} /* dz_forefront_t */;

/*
 * for backward compatibility
 */
dz_static_assert(sizeof(dz_forefront_t) == sizeof(dz_state_t));
dz_static_assert(offsetof(dz_forefront_t, cnt)   == offsetof(dz_state_t, cnt));
dz_static_assert(offsetof(dz_forefront_t, range) == offsetof(dz_state_t, range));
dz_static_assert(offsetof(dz_forefront_t, max)   == offsetof(dz_state_t, max.score.abs));


/* pointer casting */
#define dz_cff(_p)					( (dz_forefront_t const *)(_p) )
#define dz_cpff(_p)					( (dz_forefront_t const **)(_p) )

/* APIs: control */
#define dz_profile_root(_profile)	( dz_cpff(&(_profile)->root) )
#define dz_root(_self)				( dz_profile_root((_self)->profile) )

/* APIs: X-drop */
#define dz_is_terminated(_ff)		( dz_cff(_ff)->range.sblk >= dz_cff(_ff)->range.eblk )




/*
 * dz_state_merge_core_t is internal representation of dz_state_t while DP filling.
 */
typedef struct {
	__m128i ccap;		/* dz_range_t and dz_ref_cnt_t pair on xmm */
	__m128i max;		/* dz_max_t on xmm */
} dz_state_merge_core_t;

typedef struct {
	__m128i ccap;		/* dz_range_t and dz_ref_cnt_t pair on xmm */
	__m128i max;		/* dz_max_t on xmm */
	__m128i blen;
} dz_state_merge_t;


/*
 * assumptions for dz_state_merge_core, dz_state_patch_score, and dz_state_patch_blen
 */
dz_static_assert(offsetof(dz_state_t, range.sblk) == 0);
dz_static_assert(offsetof(dz_state_t, range.eblk) == 4);
dz_static_assert(offsetof(dz_state_t, max.score.inc) == 16);
dz_static_assert(offsetof(dz_state_t, max.score.abs) == 20);


static __dz_vectorize
__m128i dz_state_merge_fold(__m128i v, __m128i mask, __m128i x)
{
	__m128i const y = _mm_xor_si128(x, mask);		/* flip sign if needed, removed otherwise */
	__m128i const z = _mm_max_epi32(y, v);
	return(z);
}

#if !defined(DZ_NONGRAPH)
/*
 * merge incoming vectors; using SIMD max to avoid test-and-branch;
 */
static __dz_vectorize
dz_state_merge_core_t dz_state_merge_core(dz_state_t const **ff, size_t fcnt)
{
	/* for negating sblk to take min */
	__m128i const negv = _mm_cvtsi64_si128(INT32_MAX);	/* INT32_MAX, 0, 0, 0 */
	__m128i const zero = _mm_setzero_si128();

	/* init accumulators */
	__m128i ccap = zero;			/* flattened dz_cap_t */
	__m128i max  = negv;			/* INT32_MAX, 0, 0, 0 */

	/* fold max */
	for(size_t i = 0; i < fcnt; i++) {
		if(dz_is_terminated(ff[i])) { continue; }

		/*
		 * (sblk, eblk, ccnt, scnt) and
		 * negate sblk so that we can use pmind to take the minimum of them
		 * unsigned max is better? (FIXME)
		 */
		ccap = dz_state_merge_fold(ccap, negv,
			_mm_loadu_si128((__m128i const *)&ff[i]->range.sblk)
		);

		/*
		 * (max.score, max.inc, -, -); signed max for scores
		 */
		max  = dz_state_merge_fold(max, zero,
			_mm_loadl_epi64((__m128i const *)&ff[i]->max.score)		/* clear out upper 64bits */
		);
	}

	/* flip sblk again before save */
	return((dz_state_merge_core_t){
		.ccap = _mm_xor_si128(ccap, negv),
		.max  = max					/* [0] is kept INT32_MAX */
	});
}

#else
/* just copy for non-graph */
static __dz_vectorize
dz_state_merge_core_t dz_state_merge_core(dz_state_t const **ff, size_t fcnt)
{
	dz_unused(fcnt);

	debugblock({
		assert(fcnt == 1);
		assert(!dz_is_terminated(ff[0]));
		// assert(ff[0] == dz_root(self));
	});

	__m128i const negv = _mm_cvtsi64_si128(INT32_MAX);	/* INT32_MAX, 0, 0, 0 */
	__m128i const ccap = _mm_loadu_si128((__m128i const *)&ff[0]->range.sblk);
	__m128i const max  = _mm_loadl_epi64((__m128i const *)&ff[0]->max.score);

	return((dz_state_merge_core_t){
		.ccap = ccap,
		.max  = _mm_max_epi32(max, negv)	/* [0] is always INT32_MAX */
	});
}

#endif

static __dz_vectorize
__m128i dz_state_patch_blen(__m128i ccap, __m128i blv)
{
	/*
	 * clip range by query length;
	 * blen = (0, 0, blen, INT32_MAX) (see dz_query_t)
	 */
	__m128i const x = _mm_shuffle_epi32(blv, 0x04);		/* (INT32_MAX, INT32_MAX, blen, INT32_MAX) */
	__m128i const patched = _mm_min_epi32(ccap, x);		/* clipped */
	return(patched);
}

static __dz_vectorize
__m128i dz_state_patch_score(__m128i max, __m128i is_scan)
{
	/*
	 * is_scan[1] is -1 when not scan, otherwise zero
	 * is_scan[0] is always dz_add_ofs(0)
	 */

	/* reset inc and cap pointer */
	__m128i const patched = _mm_and_si128(max, is_scan);
	return(patched);
}

static __dz_vectorize
dz_state_merge_t dz_state_merge(dz_query_t const *query, dz_state_t const **ff, size_t fcnt, __m128i is_scan)
{
	/* merge */
	dz_state_merge_core_t const merged = dz_state_merge_core(ff, fcnt);

	/* load blen for clipping */
	__m128i const blen = _mm_loadl_epi64((__m128i const *)query);/* clear out upper 64bits */

	/* clip range and reset score */
	return((dz_state_merge_t){
		.ccap = dz_state_patch_blen(merged.ccap, blen),
		.max  = dz_state_patch_score(merged.max, is_scan),
		.blen = blen
	});
}




/*
 * working space
 */
typedef struct {
	/* 16bytes */
	dz_range_t range;
	dz_ref_cnt_t cnt;

	/* 16bytes */
	dz_max_t max;

	/* 16bytes */
	uint32_t bmax, blen;
	dz_link_t link;

	/* 16bytes; constants */
	uint32_t ofs;		/* dz_add_ofs(0) */
	uint32_t init;		/* dz_add_ofs(init) */
	uint32_t xth;		/* xdrop threshold */
	uint32_t one;		/* we expect we have 0x01ULL here */
} dz_ctx_t;


/*
 * assumptions for dz_ctx_init and dz_ctx_finalize
 */

/* range, cnt, and max have the same layout as dz_state_t */
dz_static_assert(offsetof(dz_ctx_t, range) - offsetof(dz_ctx_t, range) == offsetof(dz_state_t, range));
dz_static_assert(offsetof(dz_ctx_t, cnt)   - offsetof(dz_ctx_t, range) == offsetof(dz_state_t, cnt));
dz_static_assert(offsetof(dz_ctx_t, max)   - offsetof(dz_ctx_t, range) == offsetof(dz_state_t, max));

dz_static_assert(offsetof(dz_ctx_t, blen)  - offsetof(dz_ctx_t, bmax) == 4);
dz_static_assert(offsetof(dz_ctx_t, link)  - offsetof(dz_ctx_t, bmax) == 8);
dz_static_assert(offsetof(dz_ctx_t, init)  - offsetof(dz_ctx_t, ofs) == 4);
dz_static_assert(offsetof(dz_ctx_t, xth)   - offsetof(dz_ctx_t, ofs) == 8);
dz_static_assert(offsetof(dz_ctx_t, one)   - offsetof(dz_ctx_t, ofs) == 12);


static __dz_vectorize
void dz_ctx_init(dz_ctx_t *ctx, dz_state_merge_t const merged, dz_profile_t const *profile)
{
	/* transmute dz_state_merge_t -> dz_ctx_t */
	_mm_store_si128((__m128i *)&ctx->range, merged.ccap);		/* (range, cnt) */
	_mm_store_si128((__m128i *)&ctx->max,   merged.max);		/* (max.inc, max.abs, 0, 0) */
	_mm_store_si128((__m128i *)&ctx->bmax,  merged.blen);		/* (INT32_MAX, blen, 0, 0); clear out upper 64bits for dz_link_t */

	/* load constants */
	__m128i const ofs = _mm_load_si128((__m128i const *)&profile->ofs);
	_mm_store_si128((__m128i *)&ctx->ofs, ofs);
	return;
}

static __dz_vectorize
__m128i dz_ctx_finalize_ccap(__m128i ccap, __m128i link, __m128i ofs)
{
	/*
	 * accumulate link.cnt over cnt.column; increment cnt.section
	 */
	__m128i const z = _mm_setzero_si128();
	__m128i const x = _mm_unpackhi_epi32(ofs, link);	/* (-, -, 1, link.cnt) */
	__m128i const y = _mm_unpackhi_epi64(z, x);			/* (0, 0, 1, link.cnt) */
	return(_mm_add_epi32(ccap, y));		/* (range.sblk, range.eblk, cnt.section, cnt.column) */
}

static __dz_vectorize
__m128i dz_ctx_finalize_max(__m128i max, __m128i ofs)
{
	/*
	 * accumulate score.inc over score.abs; score.inc at max[0] and ofs(0) at cv[0]
	 */
	__m128i const z = _mm_setzero_si128();
	__m128i const x = _mm_sub_epi32(max, ofs);		/* subtract offset from inc */
	__m128i const y = _mm_unpacklo_epi32(z, x);		/* (0, inc, 0, -) */
	__m128i const w = _mm_unpacklo_epi64(y, z);		/* (0, inc, 0, 0) */
	return(_mm_add_epi32(max, w));
}

static __dz_vectorize
void dz_ctx_finalize(dz_ctx_t const *ctx, dz_state_t *tail)
{
	__m128i const ccap = _mm_load_si128((__m128i const *)&ctx->range);		/* (range, cnt) */
	__m128i const max  = _mm_load_si128((__m128i const *)&ctx->max);		/* (max.inc, max.abs, -, -) */
	__m128i const link = _mm_load_si128((__m128i const *)&ctx->bmax);		/* (INT32_MAX, blen, -, link.cnt) */
	__m128i const ofs  = _mm_load_si128((__m128i const *)&ctx->ofs);		/* (dz_add_ofs(0), dz_add_ofs(init), xth, 1) */


	dz_state_merge_core_t const v = {
		.ccap = dz_ctx_finalize_ccap(ccap, link, ofs),
		.max  = dz_ctx_finalize_max(max, ofs)
	};
	_mm_storeu_si128((__m128i *)&tail->range, v.ccap);
	_mm_storeu_si128((__m128i *)&tail->max,   v.max);

	// tail->max.score.abs -= 8 * tail->range.sblk;

	debug("rebase, inc(%d, %d), max(%d, %d)",
		dz_rm_ofs(_mm_extract_epi32(max, 0)), _mm_extract_epi32(ccap, 0),
		_mm_extract_epi32(max, 1), _mm_extract_epi32(v.max, 1)
	);
	return;
}




/*
 * fill-in vector handling
 */
typedef struct {
	uint16_t _pad1[7];
	uint16_t min;			/* init S value for internal vector */
	uint16_t init;			/* init S value for the first vector */
	uint16_t _pad2[7];

	__m128i div, dev;

	#if defined(DZ_NUCL_DIAG)
		__m128i iiv, iev;
	#elif defined(DZ_OFFSET)
		__m128i iiv, ievl, iev8;
	#else		/* vertical without offset DP */
		__m128i iiv, iev1, iev2, iev4, iev8;
	#endif
} dz_cvec_t;
dz_static_assert((sizeof(dz_cvec_t) % 16) == 0);

#if defined(DZ_OFFSET)
static __dz_vectorize
__m128i dz_cvec_fold_ievl(__m128i iev)
{
	__m128i v = iev;
	v = _mm_adds_epu16(v, _mm_slli_si128(v, 2));
	v = _mm_adds_epu16(v, _mm_slli_si128(v, 4));
	v = _mm_adds_epu16(v, _mm_slli_si128(v, 8));
	v = _mm_subs_epu16(v, iev);

	print_vector_raw(v);
	return(v);
}
#endif

static __dz_vectorize
dz_cvec_t dz_cvec_init_intl(dz_profile_t const *profile, uint16_t init)
{
	/* insertion penalties */
	__m128i const iiv = _mm_load_si128((__m128i const *)profile->iiv);
	__m128i const iev = _mm_load_si128((__m128i const *)profile->iev);
	__m128i const div = _mm_load_si128((__m128i const *)profile->div);
	__m128i const dev = _mm_load_si128((__m128i const *)profile->dev);

	return((dz_cvec_t){
		#if defined(DZ_NUCL_DIAG)
			.iiv = iiv,
			.iev = iev,

		#elif defined(DZ_OFFSET)
			.iiv  = iiv,
			.ievl = dz_cvec_fold_ievl(iev),
			.iev8 = _mm_slli_epi16(iev, 3),

		#else	/* vertical without offset DP */
			.iiv  = iiv,
			.iev1 = iev,
			.iev2 = _mm_add_epi16(iev, iev),
			.iev4 = _mm_slli_epi16(iev, 2),
			.iev8 = _mm_slli_epi16(iev, 3),
		#endif

		/* deletion penalties */
		.div  = div,
		.dev  = dev,

		/* X-drop threshold */
		.min  = dz_add_ofs(DZ_CELL_MIN),
		.init = init
	});
}

static __dz_vectorize
dz_cvec_t dz_cvec_init(dz_profile_t const *profile)
{
	return(dz_cvec_init_intl(profile, profile->init));
}




/*
 * max tracker
 */
/* X-drop and maximum tracking */
#if defined(DZ_OFFSET)
#define dz_maxv_context_intl(v, inc, xth) \
	register __m128i maxv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
	register __m128i xthv =  _mm_adds_epu16(_mm_set1_epi16((inc) - (xth)), (v)->ievl);		/* FIXME? */ \
	register __m128i iev8 = (v)->iev8; \
;

#define dz_maxv_forward() { \
	/* do nothing; for compatibility between offset DP */ \
	xthv = _mm_adds_epu16(xthv, iev8); \
}

#else
#define dz_maxv_context_intl(v, inc, xth) \
	register __m128i maxv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
	register __m128i xthv = _mm_set1_epi16((inc) - (xth)); \
;
#define dz_maxv_forward() ;

#endif

#define dz_maxv_context(v, t, xth) \
	/* actual DP value = base (implicit) + j * ge + cell value */ \
	dz_maxv_context_intl(v, (t)->score.inc, xth);		/* inc with offset */ \
;

#define dz_maxv_is_drop(curr_sv) ({ \
	/* maxv contains clip(curr_max - (prev_max - xdrop_threshold)) where clip(x) = x if x > 0 else 0 */ \
	__m128i const xtest_drop = _mm_subs_epu16(curr_sv, xthv); \
	/* all columns become zero when all the cells dropped X from the current max */ \
	print_vector_raw(xtest_drop); \
	debug("test(%u)", dz_is_all_zero(xtest_drop)); \
	(uint64_t)dz_is_all_zero(xtest_drop); \
})

#ifdef DZ_FULL_LENGTH_BONUS
#define dz_maxv_update(curr_sv, curr_bv) ({ \
	/* \
	 * update maximum-tracking vector; uses xdrop-adjusted vector again \
	 * for collecting columnwise max values on x->maxv. \
	 * appropriate full-length bonus fed to bv in positive. \
	 */ \
	__m128i const xtest = _mm_subs_epu16(curr_sv, xthv);		/* reused in is_drop */ \
	__m128i const xmaxv = _mm_add_epi16(xtest, curr_bv); \
	/* update max */ \
	maxv = _mm_max_epu16(maxv, xmaxv); \
	/* test X-drop */ \
	dz_maxv_is_drop(curr_sv); \
})

#else
#define dz_maxv_update(curr_sv, curr_bv) ({ \
	dz_unused(curr_bv); \
	/* \
	 * update maximum-tracking vector; uses xdrop-adjusted vector again \
	 * for collecting columnwise max values on x->maxv. \
	 * appropriate full-length bonus fed to bv in positive. \
	 */ \
	__m128i const xtest = _mm_subs_epu16(curr_sv, xthv);		/* reused in is_drop */ \
	/* update max */ \
	maxv = _mm_max_epu16(maxv, xtest); \
	/* test X-drop */ \
	dz_maxv_is_drop(curr_sv); \
})

#endif


typedef struct {
	uint32_t peak;
	uint32_t diff;
} dz_maxv_inc_t;

#define dz_maxv_max_inc() ({ \
	/* xthv after accumulation mostly correspond to the maximum value of DP cells */ \
	uint32_t const sample = (uint16_t)_mm_extract_epi16(xthv, 7);		/* previous max */ \
	sample; \
})
#define dz_maxv_fold_diff() ({ \
	/* \
	 * folding max over the max-tracking vector; returns increase on score and inc \
	 */ \
	__m128i t = maxv; \
	t = _mm_max_epu16(t, _mm_srli_si128(t, 8)); \
	t = _mm_max_epu16(t, _mm_srli_si128(t, 4)); \
	t = _mm_max_epu16(t, _mm_srli_si128(t, 2)); \
	uint16_t const folded_inc = _mm_extract_epi16(t, 0); \
	debug("max(%u)", folded_inc); \
	(uint32_t)folded_inc; \
})
#define dz_maxv_finalize() ({ \
	debug("finalize maxv"); \
	((dz_maxv_inc_t){ \
		/* both with offset */ \
		.peak = dz_maxv_max_inc(), \
		.diff = dz_maxv_fold_diff() \
	}); \
})




/* fetcher */

typedef struct {
	int8_t const *matrix;
	dz_query_t const *query;

	/*
	 * we don't explicitly copy fetcher object onto stack,
	 * because doing so somewhat scrapes metadata off from the original object,
	 * making unable the fetcher functions to be inlined.
	 */
	// dz_fetcher_t const fetcher;

	#ifdef DZ_NUCL_DIAG
		uint32_t max_eblk;
	#endif
} dz_fetch_work_t;

static __dz_vectorize
dz_fetch_work_t dz_fetch_init(dz_profile_t const *profile, dz_query_t const *query)
{
	return((dz_fetch_work_t){
		.matrix  = profile->matrix,
		.query   = query
	});
}

static __dz_vectorize
uint64_t dz_fetch_next(dz_fetch_work_t *f, dz_fetcher_t const *fetcher, dz_link_t *link)
{
	/* we expect this function call is inlined */
	dz_fill_fetch_t const next = fetcher->fetch_next(fetcher->opaque, f->matrix, f->query);
	debug("fetched, rch(%x), is_term(%u)", next.rch, next.is_term);
	if(next.is_term) { return(1); }

	/* FIXME: update adj?? */
	link->rch = next.rch;
	link->idx++;

	#ifdef DZ_NUCL_DIAG
		f->max_eblk = dz_roundup(link->idx, DZ_L) / DZ_L;
		debug("idx(%u), max_eblk(%u)", link->idx, f->max_eblk);
	#endif
	return(0);
}

static __dz_vectorize
uint64_t dz_fetch_clip_blen(dz_fetch_work_t *f, uint32_t blen)
{
	/* compared *after* p incremented */
	#if !defined(DZ_NUCL_DIAG)
		dz_unused(f);
		return(blen);
	#else
		return(dz_min2(f->max_eblk, blen));
	#endif
}

static __dz_vectorize
__m128i dz_fetch_get_profile(dz_fetch_work_t const *f, dz_fetcher_t const *fetcher, uint32_t p)
{
	return(fetcher->get_profile(fetcher->opaque,
		f->matrix,
		f->query,
		p * DZ_L
	));
}

static __dz_vectorize
__m128i dz_fetch_get_bonus(dz_fetch_work_t const *f, uint64_t is_tail)
{
	return(dz_query_get_bonus(f->query, is_tail));
}




/*
 * DP fill-in macros
 * all vectors and control variables are declared `flat' on stack
 * so that compilers won't generate redundant spill-reloads between
 * function calls.
 */
#if defined(DZ_BENCH_CNT)
#  define dz_bench_context()	size_t vcnt = 0;
#  define dz_bench_inc()		{ vcnt++; }
#  define dz_bench_acc(w)		{ (w)->cnt += vcnt; }
#else
#  define dz_bench_context()	;
#  define dz_bench_inc()		;
#  define dz_bench_acc(w)		;
#endif


#if defined(DZ_NUCL_DIAG)

#define dz_update_context_cvec(cvec) \
	register __m128i dev = (cvec)->dev; \
	register __m128i iev = (cvec)->iev; \
	register __m128i div = (cvec)->div; \
	register __m128i iiv = (cvec)->iiv; \
;
#define dz_update_context_uvec() \
	register __m128i e, f, s; \
	register __m128i ps;		/* prev s */ \
	register __m128i pf;		/* prev f */ \
	register __m128i next_sc; \
	dz_bench_context(); \
;
#define dz_update_context(w, ctx) \
	dz_update_context_uvec(); \
	dz_update_context_cvec(&(w)->cvec); \
	dz_range_t range = ctx->range; \
;
#define dz_update_finalize(w, ctx) { \
	ctx->range = range; \
}

#define dz_update_load_svec(second) { \
	s = _mm_load_si128((__m128i const *)&(second)->s); \
	print_vector(s); \
}
#define dz_update_load_gvec(first) { \
	e = _mm_load_si128((__m128i const *)&(first)->e); \
	f = _mm_load_si128((__m128i const *)&(first)->f); \
	print_vector(e); \
	print_vector(f); \
}
#define dz_update_reset_svec() { \
	s = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
}
#define dz_update_reset_head() {}
#define dz_update_reset_tail() { \
	/* min == zero */ \
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
	/* clear leftward dependency */ \
	e = minv; \
	f = minv; \
	s = minv; \
}
#define dz_update_reset_prev(v, is_head) { \
	/* ofs(0) if scan, min otherwise; FIXME: is cmovq -> pinsertw faster? */ \
	__m128i const psv = _mm_loadu_si128((__m128i const *)&(&(v)->min)[(size_t)(is_head) - 7]); \
	__m128i const pfv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
	ps = psv; \
	pf = pfv; \
}

#else		/* DZ_NUCL_DIAG */

#if defined(DZ_OFFSET)
#define dz_update_context_cvec(cvec) \
	register __m128i dev = (cvec)->dev; \
	register __m128i div = (cvec)->div; \
	register __m128i iiv = (cvec)->iiv; \
;

#else
#define dz_update_context_cvec(cvec) \
	register __m128i dev = (cvec)->dev; \
	register __m128i iev1 = (cvec)->iev1; \
	register __m128i iev2 = (cvec)->iev2; \
	register __m128i iev4 = (cvec)->iev4; \
	register __m128i iev8 = (cvec)->iev8; \
	register __m128i div = (cvec)->div; \
	register __m128i iiv = (cvec)->iiv; \
;
#endif

#define dz_update_context_uvec() \
	register __m128i e, f, s; \
	register __m128i ps;		/* prev s */ \
	register __m128i next_sc; \
	dz_bench_context(); \
;
#define dz_update_context(w, ctx) \
	dz_update_context_uvec(); \
	dz_update_context_cvec(&(w)->cvec); \
	dz_range_t range = ctx->range; \
;
#define dz_update_finalize(w, ctx) { \
	ctx->range = range; \
}

#define dz_update_load_svec(prev_col) { \
	s = _mm_load_si128((__m128i const *)&(prev_col)->s); \
	print_vector(s); \
}
#define dz_update_load_gvec(prev_col) { \
	debug("load, prev_col(%p)", (prev_col)); \
	/* we don't touch f vector */ \
	e = _mm_load_si128((__m128i const *)&(prev_col)->e); \
	print_vector(e); \
	print_vector(f); \
}
#define dz_update_reset_head() { \
	/* min == zero; clear upward dependency */ \
	f = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
}
#define dz_update_reset_tail() { \
	/* min == zero; clear leftward dependency */ \
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
	e = minv; \
	s = minv; \
}
#define dz_update_reset_prev(v, is_head) { \
	/* ofs(0) if scan, min otherwise; FIXME: is cmovq -> pinsertw faster? */ \
	ps = _mm_loadu_si128((__m128i const *)&(&(v)->min)[(size_t)(is_head) - 7]); \
}

#endif		/* DZ_NUCL_DIAG */

#define dz_update_get_svec()		( s )	/* s is global in this scope */
#define dz_update_store(col) ({ \
	_mm_store_si128((__m128i *)&(col)->e, e); \
	_mm_store_si128((__m128i *)&(col)->f, f); \
	_mm_store_si128((__m128i *)&(col)->s, s); \
	print_vector(s); print_vector(e); print_vector(f); \
	s; \
})
#define dz_update_load_empty(v, is_head) { \
	dz_update_reset_head(); \
	dz_update_reset_tail(); \
	dz_update_reset_prev(v, is_head);		/* prev S differs depending on coordinate */ \
}

#define dz_update_load_sc(fetch, fetcher, idx) { \
	next_sc = dz_fetch_get_profile(fetch, fetcher, idx); \
}

#if defined(DZ_NUCL_DIAG)
/* anti-diagonal slicing (for ASCII or 2bit) */
#define dz_update_body() ({ \
	dz_bench_inc(); \
	/* \
	 * E'[i - 1, j] - Ge \
	 * F'[i, j - 1] - Ge \
	 */ \
	__m128i const te = _mm_subs_epu16(e, dev); \
	__m128i const tf = _mm_subs_epu16(_mm_alignr_epi8(f, pf, 14), iev); \
	__m128i const ty = _mm_max_epu16(te, tf); \
	pf = f; \
	/* \
	 * S'[i, j] = max{ S'[i - 1, j - 1] + sc(i, j) \
	 * E'[i - 1, j] - Ge, F'[i, j - 1] - Ge } \
	 */ \
	__m128i const ov = _mm_set1_epi16(DZ_SCORE_OFS);		/* zeros */ \
	__m128i const tx = _mm_adds_epu16(_mm_alignr_epi8(s, ps, 14), next_sc); \
	__m128i const us = _mm_max_epu16(ty, _mm_subs_epu16(tx, ov));		/* to avoid overflow; FIXME: can be removed? */ \
	ps = s; \
	/* \
	 * E'[i, j] = max{ S'[i, j] - Gi, E'[i - 1, j] - Ge } \
	 * F'[i, j] = max{ S'[i, j] - Gi, F'[i, j - 1] - Ge } \
	 */ \
	__m128i const ue = _mm_max_epu16(te, _mm_subs_epu16(us, div)); \
	__m128i const uf = _mm_max_epu16(tf, _mm_subs_epu16(us, iiv)); \
	/* done */ \
	e = ue; \
	f = uf; \
	s = us; \
	/* return */ \
	s; \
})
#define dz_update_tail(v) ({ \
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN)); \
	f = minv; \
	s = minv; \
	s; \
})

#elif defined(DZ_OFFSET)
/* vertical slicling (default) */
#define dz_update_body() ({ \
	dz_bench_inc(); \
	/* \
	 * E'[i - 1, j] - Ge \
	 */ \
	__m128i const te = _mm_subs_epu16(e, dev);		/* FIXME: can be removed */ \
	/* \
	 * U[i, j] = max{ S'[i - 1, j - 1] + sc(i, j), E'[i - 1, j] - Ge } \
	 */ \
	__m128i const ov = _mm_set1_epi16(DZ_SCORE_OFS);		/* zeros */ \
	__m128i const tx = _mm_adds_epu16(next_sc, _mm_alignr_epi8(s, ps, 14)); \
	__m128i const ts = _mm_max_epu16(te, _mm_subs_epu16(tx, ov)); \
	ps = s; \
	/* \
	 * fold F' \
	 * F'[i, j] = max{ U[i, j] - Gi, F'[i, j - 1] - Ge } \
	 */ \
	__m128i tf = _mm_max_epu16( \
		_mm_subs_epu16(ts, iiv), \
		_mm_srli_si128(f, 14) \
	); \
	tf = _mm_max_epu16(tf, _mm_slli_si128(tf, 2)); \
	tf = _mm_max_epu16(tf, _mm_slli_si128(tf, 4)); \
	tf = _mm_max_epu16(tf, _mm_slli_si128(tf, 8)); \
	__m128i const uf = tf; \
	/* E'[i, j] = max{ S'[i, j] - Gi, E'[i - 1, j] - Ge } */ \
	__m128i const ue = _mm_max_epu16(te, _mm_subs_epu16(ts, div)); \
	/* S'[i, j] = max{ U[i, j], F'[i, j] } */ \
	__m128i const us = _mm_max_epu16(tf, ts); \
	/* done */ \
	e = ue; \
	f = uf; \
	s = us; \
	/* return */ \
	s; \
})
#define dz_update_tail(v) ({ \
	/* nothing to do in offset DP */ \
	s; \
})

#else	/* !defined(DZ_NUCL_DIAG) && !defined(DZ_OFFSET) */
/* vertical slicling (default) */
#define dz_update_body() ({ \
	dz_bench_inc(); \
	/* E'[i - 1, j] - Ge */ \
	__m128i const te = _mm_subs_epu16(e, dev); \
	/* U[i, j] = max{ S'[i - 1, j - 1] + sc(i, j), E'[i - 1, j] - Ge } */ \
	__m128i const ov = _mm_set1_epi16(DZ_SCORE_OFS);		/* zeros */ \
	__m128i const tx = _mm_adds_epu16(next_sc, _mm_alignr_epi8(s, ps, 14)); \
	__m128i const ts = _mm_max_epu16(te, _mm_subs_epu16(tx, ov));		/* to avoid overflow; FIXME: can be removed? */ \
	ps = s; \
	/* \
	 * fold F' \
	 * F'[i, j] = max{ U[i, j] - Gi, F'[i, j - 1] - Ge } \
	 */ \
	__m128i tf = _mm_max_epu16( \
		_mm_subs_epu16(ts, iiv), \
		_mm_subs_epu16(_mm_srli_si128(f, 14), iev1) \
	); \
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 2), iev1)); \
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 4), iev2)); \
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 8), iev4)); \
	__m128i const uf = tf; \
	/* E'[i, j] = max{ S'[i, j] - Gi, E'[i - 1, j] - Ge } */ \
	__m128i const ue = _mm_max_epu16(te, _mm_subs_epu16(ts, div)); \
	/* S'[i, j] = max{ U[i, j], F'[i, j] } */ \
	__m128i const us = _mm_max_epu16(tf, ts); \
	/* done */ \
	e = ue; \
	f = uf; \
	s = us; \
	/* return */ \
	s; \
})
#define dz_update_tail(v) ({ \
	f = _mm_subs_epu16(f, iev8); \
	s = f; \
	s; \
})

#endif




/*
 * working buffers
 */
typedef struct {
	/* constants */
	dz_fetch_work_t fetch;
	dz_cvec_t cvec;

	#if defined(DZ_BENCH_CNT)
		size_t cnt;
	#endif
} dz_work_t;




/*
 * col (array of swgv) handling
 */
#if 0
static __dz_vectorize
uint16_t *dz_col_finalize(dz_swgv_t *col, dz_range_t const *range, dz_link_t const *link)
{
	dz_cap_t *cap = dz_slice_cap(col, range->eblk);

	/*
	 * cap->range = *range;
	 * cap->link  = *link;
	 */
	uint64_t const x = dz_loadu_u64(range);
	uint64_t const y = dz_loadu_u64(link);

	dz_storeu_u64(&cap->range, x);
	dz_storeu_u64(&cap->link,  y);

	/* returns pointer to adj */
	return(&cap->link.adj);
}

static __dz_vectorize
dz_swgv_t *dz_reserve_col(dz_arena_t *mem, dz_swgv_t *prev_col, dz_range_t const *range)
{
	/* update stack pointer */
	dz_cap_t *cap = dz_slice_cap(prev_col, range->eblk);
	dz_save_stack(mem, cap + 1);

	/* add stack if memory starved */
	size_t const size = 2 * sizeof(dz_state_t *) + dz_calc_column_size(range);

	/* allocate memory from stack */
	uint8_t *ptr = dz_reserve_stack(mem, size);
	debug("cap(%p, %p)", cap, dz_cap(ptr) - 1);
	if(dz_likely(dz_cap(ptr) == cap + 1)) {
		return(dz_slice_col(cap, range->sblk));
		// return(dz_swgv(cap + 1));
	}

	/* new stack allocated (columns disjoint); create internal link */
	return(dz_slice_ilink(ptr, cap));
}

static __dz_vectorize
dz_swgv_t *dz_col_reserve_second(dz_arena_t *mem, dz_swgv_t *prev_col, dz_range_t const *range)
{
	dz_link_t const lk = { .idx = 1 };
	dz_col_finalize(prev_col, range, &lk);
	return(dz_reserve_col(mem, prev_col, range));
}
#endif

typedef struct {
	uint16_t *adj;
	dz_swgv_t *col;
} dz_reserve_col_t;

static __dz_vectorize
dz_reserve_col_t dz_reserve_col(dz_arena_t *mem, dz_swgv_t *prev_col, dz_range_t const *range, dz_link_t const *link)
{
	/* update stack pointer */
	dz_cap_t *cap = dz_slice_cap(prev_col, range->eblk);
	uint16_t *adj = dz_init_cap(cap, range, link);
	dz_save_stack(mem, cap + 1);

	/* add stack if memory starved */
	size_t const size = 2 * sizeof(dz_state_t *) + dz_calc_column_size(range);

	/* allocate memory from stack */
	uint8_t *ptr = dz_reserve_stack(mem, size);
	debug("cap(%p, %p)", cap, dz_cap(ptr) - 1);
	if(dz_likely(dz_cap(ptr) == cap + 1)) {
		return((dz_reserve_col_t){
			.adj = adj,
			.col = dz_anchor_col(dz_swgv(cap + 1), range->sblk)
		});
		// return(dz_swgv(cap + 1));
	}

	/* new stack allocated (columns disjoint); create internal link */
	return((dz_reserve_col_t){
		.adj = adj,
		.col = dz_anchor_col(dz_slice_ilink(ptr, cap), range->sblk)
	});
}

#if 0
static __dz_vectorize
dz_swgv_t *dz_col_anchor(dz_swgv_t *col, size_t sblk)
{
	return(col - sblk);
}
#endif




/*
 * previous column container
 */
typedef struct {
	#ifdef DZ_NUCL_DIAG
		dz_range_t range;
		dz_swgv_t *second;
	#endif

	/* previous one is always activated */
	dz_swgv_t *first;
} dz_prev_col_t;


static __dz_vectorize
dz_prev_col_t dz_prev_col_push(dz_prev_col_t col, dz_swgv_t *curr, dz_range_t range)
{
	dz_unused(col);
	dz_unused(range);

	return((dz_prev_col_t){
		#ifdef DZ_NUCL_DIAG
			.range  = range,
			.second = col.first,
		#endif

		.first = curr
	});
}


#if !defined(DZ_NONGRAPH)
/*
 * merge incoming vectors
 */
static __dz_vectorize
void dz_col_fold_init(dz_swgv_t *col, dz_range_t const *range)
{
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN));
	dz_swgv_t const v = { minv, minv, minv };

	for(size_t p = range->sblk; p < range->eblk; p++) {
		dz_store_swgv(&col[p], v);
	}
	return;
}

typedef struct {
	__m128i adjv;
	uint64_t underflow;
} dz_col_adjv_t;

static __dz_vectorize
dz_col_adjv_t dz_col_compose_adjv(dz_max_t const *max, dz_state_t const *incoming)
{
	/* base score at the head of this segment */
	int32_t const base = dz_score_base(&incoming->max.score);

	/* calc decrease amount from the current max */
	int32_t const gap  = dz_score_abs(&max->score) - base;
	debug("gap(%d)", gap);

	/* we are sure the value is always >= 0 for ordinary extension; always <= 0 for scan */
	__m128i const gapv = _mm_set1_epi32(gap);

	/* shrink bitwidth; clip if overflowed */
	__m128i const adjv = _mm_packus_epi32(gapv, gapv);
	return((dz_col_adjv_t){
		.adjv = adjv,
		.underflow = gap > UINT16_MAX
	});
}

static __dz_vectorize
void dz_col_fold_swgv_core(dz_swgv_t *dst, dz_swgv_t const *src, dz_range_t const *range, __m128i adjv)
{
	for(size_t p = range->sblk; p < range->eblk; p++) {
		dz_swgv_t const v = dz_load_swgv(&src[p]);
		dz_swgv_t const x = dz_load_swgv(&dst[p]);

		/* adjust offset; unsigned saturated subtraction; clipped when underflowed */
		dz_swgv_t const u = dz_max_swgv(x, dz_subs_swgv(v, adjv));
		print_vector(u.s);
		dz_store_swgv(&dst[p], u);
	}
	return;
}

static __dz_vectorize
void dz_col_fold_swgv(dz_swgv_t *col, dz_max_t const *max, dz_state_t const **ff, size_t fcnt)
{
	dz_unused(col);

	/* for each incoming vectors */
	for(size_t i = 0; i < fcnt; i++) {
		dz_state_t const *s = ff[i];
		dz_col_adjv_t const v = dz_col_compose_adjv(max, s);
		debug("adj(%d)", _mm_extract_epi16(v.adjv, 0));
		if(v.underflow) { continue; }

		/* bare and fold the last column */
		dz_swgv_t const *incoming = dz_restore_tail_column(s);
		dz_col_fold_swgv_core(col, incoming, &s->range, v.adjv);
	}
	return;
}

typedef struct {
	dz_swgv_t *col;
	size_t i;
} dz_col_shrink_t;

static __dz_vectorize
size_t dz_col_shrink_head(dz_swgv_t const *col, dz_range_t const *range)
{
	/* skip head */
	for(size_t i = range->sblk; i < range->eblk; i++) {
		dz_swgv_t const v = dz_load_swgv(&col[i]);
		if(dz_is_all_zero(v.s)) { continue; }

		return(i);
	}

	/* all vectors are zero */
	dz_trap();
	return(range->eblk);
}

static __dz_vectorize
size_t dz_col_shrink_tail(dz_swgv_t *dst, dz_swgv_t const *col, dz_range_t const *range)
{
	for(size_t i = range->sblk; i < range->eblk; i++) {
		dz_swgv_t const v = dz_load_swgv(&col[i]);
		if(dz_is_all_zero(v.s)) { return(i); }

		dz_store_swgv(&dst[i], v);
	}
	return(range->eblk);
}

static __dz_vectorize
dz_swgv_t *dz_col_shrink(dz_swgv_t *col, dz_range_t *range)
{
	/* clip head */
	size_t const sblk = dz_col_shrink_head(col, range);
	dz_swgv_t *dst = &col[range->sblk - sblk];
	range->sblk = sblk;

	/* clip tail */
	range->eblk = dz_col_shrink_tail(dst, col, range);

	// fprintf(stderr, "range(%u, %u)\n", range->sblk, range->eblk);
	return(dst);
}


static __dz_vectorize
dz_prev_col_t dz_col_merge(dz_arena_t *mem, dz_range_t *range, dz_max_t const *max, dz_state_t const **ff, size_t fcnt)
{
	/* range obtained in dz_load_state */
	dz_head_t *head = dz_slice_head(mem, range, ff, fcnt);

	/* save incoming vector array */
	dz_save_incoming(head, ff, fcnt);

	/* save head-cap info; convert to column pointer */
	dz_cap_t *cap  = dz_init_head(head, max, ff, fcnt);
	dz_swgv_t *col = dz_slice_col(cap, range->sblk);

	/* then merge */
	dz_col_fold_init(col, range);			/* clear the column */
	dz_col_fold_swgv(col, max, ff, fcnt);

	/* shrink and pack */
	dz_swgv_t *shrinked = dz_col_shrink(col, range);
	return((dz_prev_col_t){
		.first = shrinked
	});
}


#else	/* DZ_NONGRAPH */
/*
 * just copy incoming vectors
 */
static __dz_vectorize
void dz_col_copy_swgv(dz_swgv_t *dst, dz_swgv_t const *src, dz_range_t const *range)
{
	for(size_t p = range->sblk; p < range->eblk; p++) {
		dz_swgv_t const v = dz_load_swgv(&src[p]);
		dz_store_swgv(&dst[p], v);
	}
	return;
}

static __dz_vectorize
dz_prev_col_t dz_col_merge(dz_arena_t *mem, dz_range_t const *range, dz_max_t const *max, dz_state_t const **ff, size_t fcnt)
{
	debugblock({
		assert(fcnt == 1);
	});

	/* source vectors */
	dz_swgv_t const *csrc = dz_restore_tail_column(ff[0]);

	/* slice head; save incoming vector array */
	dz_head_t *head = dz_slice_head(mem, range, ff, fcnt);
	dz_save_incoming(head, ff, fcnt);

	/* allocate first vector */
	dz_swgv_t *ccol = dz_slice_col(
		dz_init_head(head, max, ff, fcnt),
		range->sblk
	);

	#if !defined(DZ_NUCL_DIAG)
		dz_col_copy_swgv(ccol, csrc, range);
		return((dz_prev_col_t){
			.first = ccol
		});

	#else
		dz_link_t const lk = { .idx = 0 };

		/* push vector; allocate another for second */
		dz_swgv_t const *psrc = dz_restore_col(dz_unwind_cap(dz_ccap(ff[0])));
		dz_reserve_col_t const r = dz_reserve_col(mem, ccol, range, &lk);

		dz_col_copy_swgv(ccol,  psrc, range);
		dz_col_copy_swgv(r.col, csrc, range);
		return((dz_prev_col_t){
			.first  = r.col,
			.second = ccol,
			.range  = *range,
		});

	#endif
}

#endif


static __dz_vectorize
void dz_col_rebase_core(dz_swgv_t *col, dz_range_t const *range, __m128i adjv)
{
	/* subtract adjustment */
	for(size_t p = range->sblk; p < range->eblk; p++) {
		/* load-compensate-store */
		dz_swgv_t const v  = dz_load_swgv(&col[p]);
		dz_swgv_t const cv = dz_subs_swgv(v, adjv);
		dz_store_swgv(&col[p], cv);
	}
	return;
}

static __dz_vectorize
void dz_col_rebase(dz_prev_col_t prev_col, dz_range_t const *range, uint32_t adj)
{
	debugblock({
		assert(adj < UINT16_MAX, "adj(%u)", adj);
	});

	/* broadcast */
	__m128i const adjv = _mm_set1_epi16(adj);

	/* rebase the second last if needed */
	#if defined(DZ_NUCL_DIAG)
		dz_col_rebase_core(prev_col.second, &prev_col.range, adjv);
	#endif

	/* rebase */
	dz_col_rebase_core(prev_col.first, range, adjv);
	return;
}



/*
 * fill-in core functions
 * are also flattened with macros
 */
#if defined(DZ_NUCL_DIAG)
#define dz_col_fill_load(w, prev_col, p) { \
	dz_update_load_gvec(&(prev_col).first[(p)]); \
	dz_update_load_svec(&(prev_col).second[(p)]); \
	if(dz_unlikely((p) >= (prev_col).range.eblk)) { \
		dz_update_reset_svec(); \
	} \
}
#else
#define dz_col_fill_load(w, prev_col, p) { \
	dz_update_load_gvec(&(prev_col).first[(p)]); \
	dz_update_load_svec(&(prev_col).first[(p)]); \
}
#endif

#define dz_col_fill_core(w, fetcher, prev_col, dst, p) ({ \
	/* load-fill-xdrop all-in-one */ \
	dz_col_fill_load(w, prev_col, p); \
	__m128i const sv = dz_update_body(); \
	dz_update_load_sc(&w->fetch, fetcher, (p) + 1); \
	dz_update_store(dst); \
	dz_maxv_is_drop(sv); \
})

typedef struct {
	uint64_t cont;
	size_t fwd;		/* forwarded length of range.sblk */
} dz_col_fill_head_t;

#define dz_col_fill_head(w, fetcher, prev_col, col) ({ \
	uint64_t cont = 1; \
	uint64_t const sblk = range.sblk; \
	dz_update_load_sc(&(w)->fetch, fetcher, sblk); \
	switch(0) do { \
		/* the first vector dropped; remove the vector and forward pointer */ \
		dz_range_shrink(&range);		/* forward sblk */ \
		/* accumulate max tracker offset */ \
		dz_maxv_forward(); \
	default: \
		/* first test if we can continue to DP fill-in; terminate when out-of-range */ \
		if(dz_range_is_term(&range)) { cont = 0; break; } \
		/* reset F vector before filling-in DP cells */ \
		dz_update_reset_head(); \
		/* DP calculation; continue the loop when xdrop */ \
	} while(dz_col_fill_core(w, fetcher, prev_col, &col[sblk], range.sblk)); \
	/* done */ \
	((dz_col_fill_head_t){ \
		.cont = cont, \
		.fwd  = range.sblk - sblk \
	}); \
})

#define dz_col_fill_body(w, fetcher, prev_col, col, blen) ({ \
	uint64_t cont = 0; \
	/* initialize counter */ \
	uint64_t p = range.sblk; \
	/* load-update-save loop; fallthrough from the fill_head */ \
	do { \
		/* update max tracker; we need bonus vector */ \
		__m128i const curr_bv = dz_fetch_get_bonus(&(w)->fetch, p + 1 == (blen)); \
		__m128i const curr_sv = dz_update_get_svec(); \
		dz_maxv_update(curr_sv, curr_bv);	/* (somewhat lazy) max update; see above */ \
		/* accumulate max tracker offset before testing range */ \
		dz_maxv_forward(); \
		/* continue to fill_tail loop */ \
		if(++p >= range.eblk) { cont = 1; break; }		/* eblk may be max_eblk */ \
		/* DP calculation; continue the loop when not xdrop */ \
	} while(!dz_col_fill_core(w, fetcher, prev_col, &col[p], p)); \
	/* write back tail when xdrop */ \
	range.eblk = p; \
	/* returns is_cont flag; do not continue fill-in when xdrop */ \
	cont; \
})

#define dz_col_fill_tail(w, fetcher, col, max_eblk) { \
	/* continued from fill_body; first test if we can fill further */ \
	if(!dz_range_is_bottom(&range, max_eblk)) { \
		/* forefront extension; clip the column length if too long */ \
		dz_update_reset_tail();		/* clear e vector */ \
		/* update the last vector */ \
		__m128i cv = dz_update_body(); \
		dz_update_store(&col[range.eblk]); \
		/* until all the cells drop; here E vector is still clipped to min */ \
		while(!dz_maxv_is_drop(cv)) { \
			/* save and update max tracker */ \
			dz_update_store(&(col)[range.eblk]); \
			dz_maxv_forward(); \
			/* forward eblk */ \
			if(dz_range_extend(&range, max_eblk)) { break; } \
			/* update vertical chain */ \
			cv = dz_update_tail(&(w)->cvec); \
			dz_update_store(&col[range.eblk]); \
		} \
	} \
}

#define dz_col_fill_init(w) ({ \
	/* init vectors */ \
	dz_update_reset_prev(&(w)->cvec, dz_range_is_head(&range)); \
	range; \
})

#define dz_col_fill_is_overflow(inc) ({ \
	(inc) > (UINT16_MAX - DZ_OVERFLOW_MARGIN); \
})

#define dz_col_fill_finalize(w, prev_col, max, xth) ({ \
	dz_bench_acc(w); \
	/* create cap object that contains [sblk, eblk) range (for use in the traceback routine) */ \
	dz_cap_t const *cap = dz_slice_cap((prev_col).first, range.eblk); \
	/* update max tracker then check overflow */ \
	dz_maxv_inc_t const inc = dz_maxv_finalize();	/* inc with offset */ \
	uint32_t inc_adj = dz_max_update(max, cap, xth, inc.diff); \
	if(dz_likely(!dz_col_fill_is_overflow(inc.peak))) { \
		inc_adj = 0; \
	} else { \
		debug("rebase needed, range(%u, %u), peak(%u, %u), inc(%u, %u)", range.sblk, range.eblk, inc.peak, dz_rm_ofs(inc.peak), inc.diff, dz_rm_ofs(inc.diff)); \
		/* rebase DP vectors */ \
		dz_max_rebase(max, inc_adj); \
		dz_col_rebase(prev_col, &range, inc_adj); \
	} \
	/* save adj (non-zero) when rebase took place */ \
	inc_adj; \
})


typedef struct {
	dz_prev_col_t prev_col;
	uint32_t adj;
} dz_col_fill_t;

static __dz_vectorize
dz_col_fill_t dz_col_fill(dz_work_t *w, dz_fetcher_t const *fetcher, dz_prev_col_t prev_col, dz_swgv_t *col, dz_ctx_t *ctx)
{
	dz_update_context(w, ctx);	/* declare vectors */
	dz_maxv_context(&(w)->cvec, &ctx->max, ctx->xth);

	/* init DP and maximum-tracking vectors */
	dz_range_t const prev_range = dz_col_fill_init(w);

	/* if reached the forefront of the query sequence, finish the extension */
	dz_col_fill_head_t const h = dz_col_fill_head(w, fetcher, prev_col, col);		/* poll head vectors and determine first survived one */
	debug("cont(%u), fwd(%zu)", (uint32_t)h.cont, h.fwd);

	/* update score offset for skipped blocks */
	dz_max_forward(&ctx->max, h.fwd);
	col -= h.fwd; 			/* dz_col_forward */

	/* then start fill-in */
	if(h.cont && dz_col_fill_body(w, fetcher, prev_col, col, ctx->blen)) {
		uint32_t const max_eblk = dz_fetch_clip_blen(&w->fetch, ctx->blen);
		dz_col_fill_tail(w, fetcher, col, max_eblk);
	}

	/* compose new prev_col */
	dz_prev_col_t const new_prev_col = dz_prev_col_push(prev_col, col, prev_range);

	/* calc and save adj in link */
	uint32_t const adj = dz_col_fill_finalize(w, new_prev_col, &ctx->max, ctx->xth);

	dz_update_finalize(w, ctx);
	return((dz_col_fill_t){
		.prev_col = new_prev_col,
		.adj      = adj
	});
}




static __dz_vectorize
dz_state_t *dz_extend_core_loop(
	dz_arena_t *mem,
	dz_ctx_t *ctx,
	dz_prev_col_t prev_col,
	dz_profile_t const *profile,
	dz_query_t const *query,
	dz_fetcher_t const *fetcher)
{
	/* working variables */
	dz_work_t w __attribute__(( aligned(16) )) = {
		.fetch = dz_fetch_init(profile, query),
		.cvec  = dz_cvec_init(profile)
	};

	/* until X-drop; try next base, break if it's the end */
	while(!dz_range_is_term(&ctx->range) && !dz_fetch_next(&w.fetch, fetcher, &ctx->link)) {
		debug("start, col(%p)", prev_col.first);

		/*
		 * create cap column for the previous column and reserve room for the next;
		 * get pointer to adj bin so that we can save adj afterward.
		 */
		dz_reserve_col_t const r = dz_reserve_col(mem, prev_col.first, &ctx->range, &ctx->link);

		/* forward one */
		dz_col_fill_t const f = dz_col_fill(&w, fetcher, prev_col, r.col, ctx);

		/* save adj and alias pointer for the next loop */
		*r.adj   = f.adj;
		prev_col = f.prev_col;
		debug("end, col(%p), adj(%d)", f.prev_col.first, f.adj);
	}

	/* place tail cap */
	dz_state_t *tail = dz_slice_tail(mem, prev_col.first, ctx->range.eblk);
	#if defined(DZ_BENCH_CNT)
		tail->vcnt = w.cnt;
	#endif
	return(tail);
}

static __dz_vectorize
dz_state_t const *dz_extend_core(
	dz_arena_t *mem,
	dz_profile_t const *profile,
	dz_query_t const *query,
	dz_fetcher_t const *fetcher,
	dz_state_t const **ff,
	size_t fcnt)
{
	if(fcnt == 0) { return(NULL); }		/* invalid input */

	#if defined(DZ_NUCL_DIAG)
		if(fcnt > 1) { return(NULL); }
	#endif

	/* merge metadata (range and max score) then vectors */
	dz_ctx_t ctx __attribute__(( aligned(16) ));
	dz_ctx_init(&ctx,
		dz_state_merge(query, ff, fcnt, dz_profile_is_scan(profile)),
		profile
	);

	/* extend */
	dz_state_t *tail = dz_extend_core_loop(mem, &ctx,
		dz_col_merge(mem, &ctx.range, &ctx.max, ff, fcnt),		/* merge DP columns */
		profile, query, fetcher
	);
	debug("done, cap(%p)", tail);

	/* save them in tail object */
	dz_ctx_finalize(&ctx, tail);
	return(tail);
}




/*
 * profile constructor
 */
typedef struct dz_alloc_profile_t {
	dz_profile_t *profile;
	dz_arena_t *mem;
} dz_alloc_profile_t;

static __dz_vectorize
size_t dz_estimate_arena_size(dz_score_conf_t const *conf)
{
	size_t const col_cnt = 2;

	/* estimate max root vector size */
	dz_range_t const range = dz_range_estimate_max(conf->max_gap_len);
	size_t const size = col_cnt * dz_calc_column_size(&range);
	return(size);
}

static __dz_vectorize
size_t dz_estimate_profile_size(dz_score_conf_t const *conf)
{
	dz_unused(conf);

	/* estimate (profile + arena) size */
	size_t const size = (
		  sizeof(dz_profile_t)
		+ DZ_REF_MAT_SIZE * DZ_QUERY_MAT_SIZE
		+ 2 * sizeof(__m128i)
	);
	return(size);
}

static __dz_vectorize
dz_alloc_profile_t dz_alloc_profile(dz_allocator_t *alloc, dz_score_conf_t const *conf)
{
	size_t const profile_size = dz_estimate_profile_size(conf);
	size_t const mem_size = dz_estimate_arena_size(conf);

	/* allocate together */
	uint8_t *ptr = alloc->fp(alloc->ctx, profile_size + mem_size);

	/* init profile */
	dz_profile_t *profile = (dz_profile_t *)ptr;
	profile->size = profile_size + mem_size;

	/* init mem */
	dz_arena_t *mem = (dz_arena_t *)&ptr[profile_size];
	dz_arena_init(mem, mem_size);
	return((dz_alloc_profile_t){
		.profile = profile,
		.mem = mem
	});
}

static __dz_vectorize
void dz_destroy_profile(dz_destructor_t *dtor, dz_profile_t *profile)
{
	dtor->fp(dtor->ctx, profile);
	return;
}


/* init contents */

#if defined(DZ_OFFSET)
static __dz_vectorize
int8_t dz_init_get_score_ofs(dz_score_conf_t const *conf)
{
	return(DZ_SCORE_OFS + conf->ins_extend);
}

#else		/* DZ_OFFSET */
static __dz_vectorize
int8_t dz_init_get_score_ofs(dz_score_conf_t const *conf)
{
	dz_unused(conf);
	return(DZ_SCORE_OFS);
}

#endif

static __dz_vectorize
void dz_init_score_matrix(dz_profile_t *profile, dz_score_conf_t const *conf)
{
	int8_t const *score_matrix = conf->score_matrix;
	int8_t const ofs = dz_init_get_score_ofs(conf);
	debug("ofs(%d)", ofs);

	/* transpose if needed */
	for(size_t i = 0; i < DZ_QUERY_MAT_SIZE; i++) {
		for(size_t j = 0; j < DZ_REF_MAT_SIZE; j++) {

			#if defined(DZ_TRANSPOSE_MATRIX) && (DZ_TRANSPOSE_MATRIX != 0)
				int8_t *dst = &profile->matrix[j * DZ_QUERY_MAT_SIZE + i];
			#else
				int8_t *dst = &profile->matrix[i * DZ_REF_MAT_SIZE + j];
			#endif

			/* copy */
			int8_t const *src = &score_matrix[i * DZ_REF_MAT_SIZE + j];
			*dst = *src + ofs;
		}
	}

	/* padding invalid reference base */ {
		int8_t *pad = &profile->matrix[DZ_REF_MAT_SIZE * DZ_QUERY_MAT_SIZE];
		_mm_store_si128((__m128i *)pad, _mm_setzero_si128());
	}
	return;
}

static __dz_vectorize
void dz_init_gap_penalties(dz_profile_t *profile, dz_score_conf_t const *conf)
{
	/* insertion penalties */
	__m128i const ii = _mm_set1_epi16(conf->ins_open);
	__m128i const ie = _mm_set1_epi16(conf->ins_extend);
	_mm_store_si128((__m128i *)profile->iiv, ii);
	_mm_store_si128((__m128i *)profile->iev, ie);

	/* deletion penalties */
	__m128i const di = _mm_set1_epi16(conf->del_open);
	__m128i const de = _mm_set1_epi16(conf->del_extend);
	_mm_store_si128((__m128i *)profile->div, di);
	_mm_store_si128((__m128i *)profile->dev, de);
	return;
}

static __dz_vectorize
void dz_init_misc_constants(dz_profile_t *profile, dz_score_conf_t const *conf)
{
	/* gap lengths */
	profile->max_gap_len = conf->max_gap_len;	/* save raw value */
	profile->bonus       = conf->full_length_bonus;

	/* X-drop threshold */
	uint16_t const gi = dz_max2(conf->ins_open,   conf->del_open);
	uint16_t const ge = dz_max2(conf->ins_extend, conf->del_extend);

	/* (ofs, init, xth, 1) */
	profile->ofs  = dz_add_ofs(0);
	profile->init = dz_add_ofs(DZ_CELL_MIN);
	profile->xth  = gi + ge * conf->max_gap_len;	/* FIXME: X-drop threshold */
	profile->one  = 1;

	debug("gi(%u), ge(%u), xdrop_threshold(%u), full_length_bonus(%u), max_gap_len(%u)", gi, ge, profile->xth, profile->bonus, profile->max_gap_len);
	return;
}


#if !defined(DZ_NUCL_DIAG)

static __dz_vectorize
dz_range_t dz_init_fill_root(dz_profile_t const *profile, dz_prev_col_t col)
{
	/*
	 * constant vectors; always use ofs(0) for the root cell
	 * load vectors; gap vectors must be initialized in advance
	 */
	dz_cvec_t cvec __attribute__(( aligned(16) )) = dz_cvec_init_intl(profile, dz_add_ofs(0));

	/* working vectors */
	dz_maxv_context_intl(&cvec, dz_add_ofs(0), profile->xth);
	dz_unused(maxv);

	dz_update_context_uvec();
	dz_update_context_cvec(&cvec);
	dz_update_load_empty(&cvec, 1);	/* is_head == 1; root cell loaded here */

	/* the first vector always survives because the origin has score zero */
	next_sc = _mm_set1_epi16(DZ_SCORE_OFS);
	dz_update_body();	/* score == 0; maxv updated */
	dz_update_store(&col.first[0]);	/* X-drop never be true */

	/* p = 1 */ {
		dz_update_reset_tail();

		dz_maxv_forward();		/* for DZ_OFFSET */
		dz_update_body();

		__m128i const curr_sv = dz_update_store(&col.first[1]);
		if(dz_maxv_is_drop(curr_sv)) {
			return((dz_range_t){
				.sblk = 0,
				.eblk = 1
			});
		}
	}

	/* until X-drop */
	for(size_t p = 2;; p++) {
		dz_maxv_forward();		/* for DZ_OFFSET */
		dz_update_tail(&cvec);

		__m128i const curr_sv = dz_update_store(&col.first[p]);
		if(dz_maxv_is_drop(curr_sv)) {
			return((dz_range_t){
				.sblk = 0,
				.eblk = p
			});
		}
	}

	/* never reach here */
	dz_trap();
	return((dz_range_t){ .sblk = 0 });
}
#else	/* DZ_NUCL_DIAG */

static __dz_vectorize
dz_range_t dz_init_fill_root(dz_profile_t const *profile, dz_prev_col_t col)
{
	/*
	 * constant vectors; always use ofs(0) for the root cell
	 * load vectors; gap vectors must be initialized in advance
	 */
	dz_cvec_t cvec __attribute__(( aligned(16) )) = dz_cvec_init_intl(profile, dz_add_ofs(0));

	/* working vectors */
	dz_update_context_uvec();
	dz_update_context_cvec(&cvec);

	/* first */ {
		next_sc = _mm_set1_epi16(DZ_SCORE_OFS);

		dz_update_load_empty(&cvec, 1);
		dz_update_body();
		dz_update_store(&col.second[col.range.sblk]);
	}

	/* create link */ {
		dz_link_t const lk = { .idx = 0 };
		dz_cap_t *cap = dz_slice_cap(col.second, col.range.eblk);
		dz_init_cap(cap, &col.range, &lk);
	}

	/* second */ {
		next_sc = _mm_setzero_si128();
		__m128i const minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN));
		dz_swgv_t const empty = { minv, minv, minv };

		dz_update_reset_head();
		dz_update_reset_prev(&cvec, 0);
		dz_update_load_gvec(&col.second[col.range.sblk]);
		dz_update_load_svec(&empty);

		dz_update_body();
		dz_update_store(&col.first[col.range.sblk]);
	}
	return((dz_range_t){
		.sblk = 0,
		.eblk = 1
	});
}
#endif

static __dz_vectorize
dz_range_t dz_init_calc_range(size_t max_gap_len)
{
	#if defined(DZ_NUCL_DIAG)
		max_gap_len = 0;
	#endif

	/* roundup */
	max_gap_len = dz_roundup(max_gap_len, DZ_L);

	return((dz_range_t){
		.sblk = 0,
		.eblk = max_gap_len / DZ_L + 1
	});
}

static __dz_vectorize
dz_prev_col_t dz_init_slice_root(dz_arena_t *mem, dz_range_t const *range)
{
	/* no incoming vector for the root */
	dz_head_t *head = dz_slice_head(mem, range, NULL, 0);

	/* mark root */
	dz_cap_t *cap  = dz_init_root_head(head);
	dz_swgv_t *col = dz_slice_col(cap, range->sblk);		/* sblk == 0 */

	#if !defined(DZ_NUCL_DIAG)
		return((dz_prev_col_t){
			.first = col
		});
	#else
		dz_link_t const lk = { .idx = 0 };

		dz_reserve_col_t const r = dz_reserve_col(mem, col, range, &lk);
		return((dz_prev_col_t){
			/* internal cap placed */
			.first  = r.col,
			.second = col,
			.range  = *range
		});
	#endif
}

static __dz_vectorize
dz_forefront_t const *dz_init_root(dz_arena_t *mem, dz_profile_t const *profile, dz_score_conf_t const *conf)
{
	/* calc vector length; query = NULL for the first (root) column */
	dz_range_t const max_range = dz_init_calc_range(conf->max_gap_len);

	/* slice first (and second if needed) */
	dz_prev_col_t const col = dz_init_slice_root(mem, &max_range);
	dz_range_t const range = dz_init_fill_root(profile, col);

	/* tail object */
	dz_state_t *root = dz_state(dz_slice_cap(col.first, range.eblk));
	dz_init_tail(root, &range);
	debug("col(%p), range(%u, %u), tail(%p)", col.first, root->range.sblk, root->range.eblk, root);

	return(dz_cff(root));
}

static __dz_vectorize
dz_profile_t *dz_init_profile(dz_allocator_t *alloc, dz_score_conf_t const *conf)
{
	/* allocate score priflie object and fill the root column */
	dz_alloc_profile_t const p = dz_alloc_profile(alloc, conf);
	dz_init_score_matrix(p.profile, conf);
	dz_init_gap_penalties(p.profile, conf);
	dz_init_misc_constants(p.profile, conf);

	/* allocate first column and fill cells */
	p.profile->root = dz_init_root(p.mem, p.profile, conf);
	debug("profile(%p), root(%p)", p.profile, p.profile->root);
	return(p.profile);
}




/*
 * maximum-scoring position search
 */
typedef struct {
	uint32_t rpos, qpos;
	int32_t score, bonus;
} dz_max_pos_t;

static __dz_vectorize
uint32_t dz_calc_max_rpos_core(dz_profile_t const *profile, dz_query_t const *query, dz_state_t const *ff)
{
	dz_unused(profile);
	dz_unused(query);

	/* max-scoring column */
	dz_cap_t const *cap = ff->max.cap;
	if(cap == NULL) { return(0); }

	// dz_tail_t const *tail = dz_restore_tail(ff);
	dz_cap_t const *prev_cap = dz_unwind_cap(cap);
	uint32_t const rpos = (dz_is_head(prev_cap) ? 0 : prev_cap->link.idx);

	debug("tail(%p), prev_cap(%p, %u, %u, %u, %u), is_head(%u), rpos(%u)",
		ff, prev_cap, prev_cap->range.sblk, prev_cap->range.eblk, prev_cap->link.rch, prev_cap->link.idx,
		(uint32_t)dz_is_head(prev_cap), rpos
	);
	return(rpos);
}

typedef struct {
	uint32_t qpos;
	int32_t bonus;
} dz_max_qpos_t;

#if defined(DZ_OFFSET)
static __dz_vectorize
__m128i dz_qpos_load_adjv(dz_profile_t const *profile)
{
	__m128i const iev  = _mm_load_si128((__m128i const *)profile->iev);
	__m128i const iev8 = _mm_slli_epi16(iev, 3);


	print_vector_raw(iev8);
	return(iev8);
}

static __dz_vectorize
__m128i dz_qpos_load_maxv(dz_profile_t const *profile, dz_cap_t const *pcap, dz_state_t const *tail)
{
	/* offset already included in tail->max.inc */
	uint16_t const ie = profile->iev[0];
	uint32_t const diff_sblk = tail->range.sblk - pcap->range.sblk;
	__m128i const rawv = _mm_set1_epi16(
		  tail->max.score.inc
		- DZ_L * diff_sblk * (uint32_t)ie
	);

	debug("adjust sblk diff, inc(%u), sblk(%u, %u), ie(%u)", tail->max.score.inc, pcap->range.sblk, tail->range.sblk, profile->iev[0]);

	/* add grad */
	__m128i const iev  = _mm_load_si128((__m128i const *)profile->iev);
	__m128i const ievl = dz_cvec_fold_ievl(iev);
	__m128i const maxv = _mm_adds_epu16(rawv, ievl);

	print_vector(rawv);
	print_vector(maxv);
	return(maxv);
}

#else		/* DZ_OFFSET */
static __dz_vectorize
__m128i dz_qpos_load_adjv(dz_profile_t const *profile)
{
	dz_unused(profile);

	__m128i const z = _mm_setzero_si128();
	return(z);
}

static __dz_vectorize
__m128i dz_qpos_load_maxv(dz_profile_t const *profile, dz_cap_t const *pcap, dz_state_t const *tail)
{
	dz_unused(profile);
	dz_unused(pcap);

	__m128i const maxv = _mm_set1_epi16(tail->max.score.inc);	/* offset already included in tail->max.inc */
	return(maxv);
}

#endif		/* DZ_OFFSET */

typedef struct {
	uint64_t eq_mask;
	uint64_t is_tail;
} dz_qpos_match_t;

static __dz_vectorize
dz_qpos_match_t dz_qpos_match(dz_query_t const *query, dz_swgv_t const *col, uint32_t p, __m128i maxv)
{
	__m128i const v = _mm_load_si128(&col[p].s);

	/* adjust full-length bonus if needed */
	#ifdef DZ_FULL_LENGTH_BONUS
		uint64_t const is_tail = p + 1 == query->blen;
		__m128i const b = dz_query_get_bonus(query, is_tail);
		__m128i const s = _mm_add_epi16(v, b);

	#else
		dz_unused(query);

		uint64_t const is_tail = 0;
		__m128i const s = v;
	#endif

	print_vector(s);
	print_vector(maxv);

	__m128i const eqv = _mm_cmpeq_epi16(s, maxv);
	return((dz_qpos_match_t){
		.eq_mask = _mm_movemask_epi8(eqv),
		.is_tail = is_tail
	});
}

static __dz_vectorize
dz_max_qpos_t dz_qpos_finalize(dz_query_t const *query, size_t p, uint64_t eq, uint64_t is_tail)
{
	/* tzcntq is faster but avoid using it b/c it requires relatively newer archs */
	uint64_t zcnt = (eq - 1) & (~eq & 0x5555);	/* subq, andnq, andq; chain length == 2 */
	zcnt += zcnt>>2; zcnt &= 0x3333;			/* shrq, addq, andq; chain length might be shorter */
	zcnt += zcnt>>4; zcnt &= 0x0f0f;
	zcnt += zcnt>>8; zcnt &= 0x00ff;
	debug("found, p(%zu), eq(%zx), zcnt(%zu), idx(%zu)", p, (size_t)eq, (size_t)zcnt, (size_t)(p * DZ_L + zcnt));

	uint32_t const qpos = p * DZ_L + zcnt;
	return((dz_max_qpos_t){
		.qpos  = qpos,
		.bonus = query->bonus[zcnt + (is_tail ? DZ_L : 0)]
	});
}

static __dz_vectorize
dz_max_qpos_t dz_calc_max_qpos_core(dz_profile_t const *profile, dz_query_t const *query, dz_state_t const *tail)
{
	dz_max_qpos_t const fail = {
		.qpos  = 0,
		.bonus = 0
	};

	dz_cap_t const *pcap = tail->max.cap;
	if(pcap == NULL) { return(fail); }

	/* load max score */
	__m128i const adjv = dz_qpos_load_adjv(profile);
	__m128i maxv = dz_qpos_load_maxv(profile, pcap, tail);

	debug("range(%u, %u)", pcap->range.sblk, pcap->range.eblk);

	/* loac column pointer */
	dz_swgv_t const *col = dz_restore_col(pcap);
	for(size_t p = pcap->range.sblk; p < pcap->range.eblk; p++) {
		dz_qpos_match_t const m = dz_qpos_match(query, col, p, maxv);
		debug("p(%zu), m(%lx)", p, m.eq_mask);
		if(m.eq_mask != 0) {
			return(dz_qpos_finalize(query, p, m.eq_mask, m.is_tail));
		}

		/* forward maxv */
		maxv = _mm_adds_epu16(maxv, adjv);
	}
	return(fail);
}

static __dz_vectorize
dz_max_pos_t dz_calc_max_pos_core(dz_profile_t const *profile, dz_query_t const *query, dz_state_t const *tail)
{
	dz_cap_t const *pcap = tail->max.cap;
	if(pcap == NULL) {
		/* is root */
		debug("pcap == NULL");
		return((dz_max_pos_t){ 0, 0 });
	}

	debug("tail(%p), pcap(%p)", tail, pcap);
	dz_max_qpos_t const q = dz_calc_max_qpos_core(profile, query, tail);
	return((dz_max_pos_t){
		.rpos  = dz_calc_max_rpos_core(profile, query, tail),
		.qpos  = q.qpos,
		.bonus = q.bonus,
		.score = dz_score_abs(&tail->max.score)
	});
}




/*
 * forward declaration for id; used in traceback
 */
static uint32_t dz_extract_id(dz_state_t const *tail);




/*
 * traceback
 */

/* alignment path object */
typedef struct dz_path_span_s {
	uint32_t id;
	uint32_t offset;
} dz_path_span_t;

/* alignment object */
/* typedef */
struct dz_alignment_s {
	int32_t score;
	int32_t bonus;

	/* total alignment length on each sequence */
	uint32_t ref_length, query_length;

	/* path */
	dz_path_span_t const *span;
	uint8_t const *path;
	uint32_t span_length, path_length;

	/* gap counts */
	uint32_t mismatch_count, match_count, ins_count, del_count;
} /* dz_alignment_t */;

/* working buffer */
typedef struct {
	/* dst arrays */
	dz_alignment_t *aln;
	struct {
		dz_path_span_t *ptr, *base;
	} span;
	struct {
		uint8_t *ptr, *base;
	} path;

	/* lengths */
	size_t rlen;
	dz_query_t const *query;

	/* match / gap counters */
	uint32_t cnt[4];

	#if defined(DZ_NUCL_DIAG)
		dz_cap_t const *pcap, *gcap, *ccap;

	#else
		dz_cap_t const *pcap, *ccap;

	#endif

	size_t idx;
	uint16_t score, ie, de;
	// uint16_t rch;
} dz_trace_work_t;

static __dz_vectorize
size_t dz_trace_vector_idx(size_t idx)
{
	return(idx / DZ_L);
}

static __dz_vectorize
size_t dz_trace_cell_idx(size_t idx)
{
	return(idx & (DZ_L - 1));
}

static __dz_vectorize
uint16_t dz_trace_score(size_t layer, dz_cap_t const *cap, size_t idx)
{
	dz_swgv_t const *col = dz_restore_col(cap);

	/* vector */
	size_t const vidx = dz_trace_vector_idx(idx);
	uint16_t const *v = (uint16_t const *)&col[vidx];

	/* cell */
	size_t const cidx = dz_trace_cell_idx(idx);
	return(v[layer * DZ_L + cidx]);
}

#if !defined(DZ_NUCL_DIAG)
static __dz_vectorize
uint16_t dz_trace_score_adj(dz_trace_work_t *w)
{
	/* w->score is of w->ccap; corresponding adj is saved one before */
	debug("score(%d), adj(%d)", w->score, w->pcap->link.adj);
	return(w->score + w->pcap->link.adj);
}

#else
static __dz_vectorize
uint16_t dz_trace_score_adj(dz_trace_work_t *w)
{
	/* w->score is of w->ccap; corresponding adj is saved one before */
	debug("score(%d), adj(%d, %d)", w->score, w->pcap->link.adj, w->gcap->link.adj);
	return(w->score + w->gcap->link.adj + w->ccap->link.adj);
}

#endif

static __dz_vectorize
uint16_t dz_trace_score_raw(dz_trace_work_t *w)
{
	/* raw score for internal comparison (for vertical (F; ins) transition) */
	return(w->score);
}

static __dz_vectorize
void dz_trace_push_span(dz_trace_work_t *w, size_t id)
{
	/* finzlize current section */
	w->span.ptr->offset = w->path.base - w->path.ptr;
	debug("offset(%u)", w->span.ptr->offset);

	/* open new bin */
	w->span.ptr--;
	w->span.ptr->id = id;
	return;
}

static __dz_vectorize
int32_t dz_trace_merge_score(dz_state_t const *tail, size_t idx, size_t layer)
{
	/* restore offsetted score; adj == 0 for head */
	// int32_t const base = tail->max.score - dz_rm_ofs(tail->max.score.inc);
	int32_t const base = dz_score_base(&tail->max.score);
	int32_t const s = dz_trace_score(layer, dz_ccap(tail), idx) + base;
	return(s);
}

static __dz_vectorize
uint64_t dz_trace_rewind_merge(dz_trace_work_t *w, size_t layer)
{
	/* merging vector; load contents to find an edge */
	dz_head_t const *head = dz_chead(w->pcap);
	int32_t const prev_score = dz_trace_score(layer, w->ccap, w->idx) + dz_score_base(&head->base);
	debug("head(%p), prev_score(%u)", head, dz_rm_ofs(prev_score));

	/* load incoming vectors */
	size_t const vidx = dz_trace_vector_idx(w->idx);
	size_t const fcnt = dz_link_fcnt(&head->link);
	dz_state_t const **farr = dz_cpstate(w->pcap) - fcnt;
	for(size_t i = 0; i < fcnt; i++) {
		dz_state_t const *tail = farr[i];
		debug("i(%zu), tail(%p)", i, tail);

		/* skip if the cell is out of the vector */
		if(!dz_inside(tail->range.sblk, vidx, tail->range.eblk)) { continue; }

		/* adj[i] = w.max - (ffs[i]->max - ffs[i]->inc); base = max - inc */
		int32_t const score = dz_trace_merge_score(tail, w->idx, layer);
		debug("score(%u), prev_score(%u)", dz_rm_ofs(score), dz_rm_ofs(prev_score));
		if(prev_score == score) {
			w->pcap = dz_ccap(tail);
			return(1);
		}
	}
	return(0);
}

static __dz_vectorize
uint64_t dz_trace_reload_section(dz_trace_work_t *w, size_t layer)
{
	/* just load pcap if internal bridge */
	if(dz_is_ilink(w->pcap)) {
		w->pcap = dz_rewind_ilink(w->pcap);
		debug("rewind, pcap(%p)", w->pcap);
		return(0);
	}

	/* load merging vector */
	if(dz_trace_rewind_merge(w, layer)) {
		return(1);
	}

	/* something is wrong... */
	dz_trap();
	return(1);
}

#if !defined(DZ_NUCL_DIAG)
static __dz_vectorize
void dz_trace_unwind_cap(dz_trace_work_t *w)
{
	w->ccap = w->pcap;
	w->pcap = dz_unwind_cap(w->pcap);
	return;
}

static inline dz_cap_t const *dz_trace_mcap(dz_trace_work_t *w) { return(w->pcap); }		/* diag, cross-line */
static inline dz_cap_t const *dz_trace_icap(dz_trace_work_t *w) { return(w->ccap); }		/* up, in-lane */
static inline dz_cap_t const *dz_trace_dcap(dz_trace_work_t *w) { return(w->pcap); }		/* left, cross-lane */

#else		/* DZ_NUCL_DIAG */
static __dz_vectorize
void dz_trace_unwind_cap(dz_trace_work_t *w)
{
	w->ccap = w->gcap;
	w->gcap = w->pcap;
	w->pcap = dz_unwind_cap(w->pcap);
	return;
}

static inline dz_cap_t const *dz_trace_mcap(dz_trace_work_t *w) { return(w->pcap); }
static inline dz_cap_t const *dz_trace_icap(dz_trace_work_t *w) { return(w->gcap); }
static inline dz_cap_t const *dz_trace_dcap(dz_trace_work_t *w) { return(w->gcap); }

#endif		/* DZ_NUCL_DIAG */

static __dz_vectorize
uint64_t dz_trace_unwind_core(dz_trace_work_t *w, size_t layer)
{
	debug("pcap(%p), cap(%p)", w->pcap, w->ccap);
	dz_trace_unwind_cap(w);

	while(1) {
		debug("pcap(%p, %p), cap(%p), is_head(%u), is_root(%u)", w->pcap, dz_trace_dcap(w), w->ccap, (uint32_t)dz_is_head(w->pcap), (uint32_t)dz_is_root(w->pcap));
		if(dz_likely(!dz_is_head(w->pcap))) { break; }	/* escape if not head */
		if(dz_is_root(w->pcap)) { break; }				/* escape if root */

		/* pcap is valid merging node */
		if(!dz_trace_reload_section(w, layer)) { break; }

		/* push segment info */
		uint32_t const id = dz_extract_id(dz_cstate(w->pcap));
		dz_trace_push_span(w, id);

		/* reload score */
		w->score = dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx);
		w->pcap  = dz_unwind_cap(w->pcap);
	}
	return(0);
}

static __dz_vectorize
uint64_t dz_trace_unwind_h(dz_trace_work_t *w, size_t layer)
{
	w->rlen++;
	return(dz_trace_unwind_core(w, layer));
}

#if !defined(DZ_NUCL_DIAG)
static __dz_vectorize
uint64_t dz_trace_unwind_v(dz_trace_work_t *w, size_t layer)
{
	dz_unused(layer);

	w->idx--;
	return(0);
}

#else		/* DZ_NUCL_DIAG */
static __dz_vectorize
uint64_t dz_trace_unwind_v(dz_trace_work_t *w, size_t layer)
{
	/* decrement lane index first */
	w->idx--;
	return(dz_trace_unwind_core(w, layer));
}

#endif		/* DZ_NUCL_DIAG */

static __dz_vectorize
uint64_t dz_trace_test_idx(dz_trace_work_t *w, dz_cap_t const *cap, size_t ofs)
{
	size_t const vidx = dz_trace_vector_idx(w->idx - ofs);
	return(!dz_inside(cap->range.sblk, vidx, cap->range.eblk));
}



static __dz_vectorize
void dz_trace_push_op(dz_trace_work_t *w, uint64_t op, uint16_t next_score)
{
	*--w->path.ptr = DZ_CIGAR_INTL>>(op<<3);
	w->cnt[op]++;
	w->score = next_score;
	debug("score(%u, %d), op(%c), idx(%zu)", next_score, dz_rm_ofs(next_score), *w->path.ptr, w->idx);
	return;
}


static __dz_vectorize
uint64_t dz_trace_eat_match(dz_trace_work_t *w, dz_profile_t const *profile, dz_getter_t const *getter)
{
	if(dz_trace_test_idx(w, dz_trace_mcap(w), 1)) { debug("test match out of range, idx(%zu), range(%u, %u)", w->idx, dz_trace_mcap(w)->range.sblk, dz_trace_mcap(w)->range.eblk); return(0); }

	/* get diagonal score for this cell */
	uint16_t const s = dz_trace_score(DZ_S_MATRIX, dz_trace_mcap(w), w->idx - 1);
	dz_trace_match_t m = getter->get_match(getter->opaque, profile->matrix, w->query, w->idx, dz_ctlink(&w->pcap->link));

	/* skip if score does not match */
	debug("test match, rch(%x), adj_score(%d, %d), match(%d, %u)", w->pcap->link.rch, dz_rm_ofs(dz_trace_score_adj(w)), dz_rm_ofs(s + m.score - DZ_SCORE_OFS), m.score, m.match);
	if(dz_trace_score_adj(w) != (uint16_t)(s + m.score - DZ_SCORE_OFS)) { return(0); }

	/* determine match state; compensate adj before saving */
	dz_trace_push_op(w, DZ_S_MATRIX + m.match, s);

	/* unwind both */
	dz_trace_unwind_v(w, DZ_S_MATRIX);
	dz_trace_unwind_h(w, DZ_S_MATRIX);
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_eat_ins(dz_trace_work_t *w)
{
	if(dz_trace_test_idx(w, dz_trace_icap(w), 1)) { debug("test ins out of range, idx(%zu), range(%u, %u)", w->idx, w->ccap->range.sblk, w->ccap->range.eblk); return(0); }

	/* skip if score does not match */
	uint16_t const f = dz_trace_score(DZ_F_MATRIX, w->ccap, w->idx);
	if(dz_likely(dz_trace_score_raw(w) != f)) { debug("test ins score unmatch, idx(%zu), raw_score(%d), f(%u)", w->idx, dz_rm_ofs(dz_trace_score_raw(w)), f); return(0); }
	debug("ins, raw_score(%d), f(%u), idx(%zu)", dz_rm_ofs(dz_trace_score_raw(w)), dz_rm_ofs(f), w->idx);

	/* go vertical at least once; unwind row to test the next condition */
	dz_trace_unwind_v(w, DZ_F_MATRIX);

	/* back to S matrix if there is no room for vertical gap */
	while(!dz_trace_test_idx(w, dz_trace_icap(w), 1)) {
		uint16_t const x = dz_trace_score(DZ_F_MATRIX, w->ccap, w->idx - 1);
		debug("ins, raw_score(%d), x(%u), ie(%u), idx(%zu)", dz_rm_ofs(dz_trace_score_raw(w)), dz_rm_ofs(x), w->ie, w->idx);

		if(dz_trace_score_raw(w) != x - w->ie) { break; }

		/* test passed; we can go once more */
		dz_trace_push_op(w, DZ_F_MATRIX, x);		/* push op for the current row */
		dz_trace_unwind_v(w, DZ_F_MATRIX);			/* load the next row */
	}

	/* push op for the last column */
	dz_trace_push_op(w, DZ_F_MATRIX, dz_trace_score(DZ_S_MATRIX, w->ccap, w->idx));
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_eat_del(dz_trace_work_t *w)
{
	if(dz_trace_test_idx(w, dz_trace_dcap(w), 0)) { debug("test ins out of range, idx(%zu), range(%u, %u)", w->idx, w->ccap->range.sblk, w->ccap->range.eblk); return(0); }

	/* skip if score does not match */
	uint16_t const e = dz_trace_score(DZ_E_MATRIX, w->ccap, w->idx);
	if(dz_likely(dz_trace_score_raw(w) != e)) { debug("test ins score unmatch, idx(%zu), raw_score(%d), e(%u)", w->idx, dz_rm_ofs(dz_trace_score_raw(w)), e); return(0); }
	debug("ins, raw_score(%d), e(%u), idx(%zu)", dz_rm_ofs(dz_trace_score_raw(w)), dz_rm_ofs(e), w->idx);

	/* go vertical at least once; unwind row to test the next condition */
	dz_trace_unwind_h(w, DZ_E_MATRIX);

	/* back to S matrix if there is no room for vertical gap */
	while(!dz_trace_test_idx(w, dz_trace_dcap(w), 0)) {
		uint16_t const x = dz_trace_score(DZ_E_MATRIX, w->ccap, w->idx - 1);
		debug("ins, raw_score(%d), x(%u), ie(%u), idx(%zu)", dz_rm_ofs(dz_trace_score_raw(w)), dz_rm_ofs(x), w->ie, w->idx);

		if(dz_trace_score_raw(w) != x - w->ie) { break; }

		/* test passed; we can go once more */
		dz_trace_push_op(w, DZ_E_MATRIX, x);		/* push op for the current row */
		dz_trace_unwind_h(w, DZ_E_MATRIX);			/* load the next row */
	}

	/* push op for the last column */
	dz_trace_push_op(w, DZ_E_MATRIX, dz_trace_score(DZ_S_MATRIX, w->ccap, w->idx));
	return(1);
}
#if 0
static __dz_vectorize
uint64_t dz_trace_eat_del(dz_trace_work_t *w)
{
	if(dz_trace_test_idx(w, w->pcap, 0)) { debug("test del out of range, idx(%zu), range(%u, %u)", w->idx, w->pcap->range.sblk, w->pcap->range.eblk); return(0); }

	/* skip if score does not match */
	uint16_t const e = dz_trace_score(DZ_E_MATRIX, w->ccap, w->idx);
	if(dz_likely(dz_trace_score_raw(w) != e)) { debug("test del score unmatch, idx(%zu), raw_score(%d), e(%u)", w->idx, dz_rm_ofs(dz_trace_score_raw(w)), e); return(0); }
	debug("del, raw_score(%d), e(%u)", dz_rm_ofs(dz_trace_score_raw(w)), dz_rm_ofs(e));

	do {
		uint16_t const x = dz_trace_score(DZ_E_MATRIX, w->pcap, w->idx);
		debug("del, adj_score(%d), x(%u), de(%u)", dz_rm_ofs(dz_trace_score_adj(w)), dz_rm_ofs(x), w->de);
		if(dz_trace_score_adj(w) != x - w->de) { break; }

		dz_trace_push_op(w, DZ_E_MATRIX, x);
		dz_trace_unwind_h(w, DZ_E_MATRIX);
	} while(!dz_trace_test_idx(w, w->pcap, 0));

	/* eat last row */
	dz_trace_push_op(w, DZ_E_MATRIX, dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx));
	dz_trace_unwind_h(w, DZ_S_MATRIX);
	return(1);
}
#endif

static __dz_vectorize
void dz_trace_init_aln(dz_alignment_t *aln, dz_state_t const *tail, dz_max_qpos_t q)
{
	// aln->rrem  = dz_calc_max_rpos_core(tail);
	aln->score = dz_score_abs(&tail->max.score);
	aln->bonus = q.bonus;
	aln->query_length = q.qpos;

	debug("tail(%p), idx(%zu), score(%d)", tail, (size_t)q.qpos, tail->max.score.abs);
	return;
}

static __dz_vectorize
void dz_trace_allocate_aln(dz_trace_work_t *w, dz_arena_t *mem, dz_state_t const *tail, dz_max_qpos_t q)
{
	/* allocate aln object */
	size_t aln_size = (sizeof(dz_alignment_t)
		+ (tail->cnt.section + 6) * sizeof(dz_path_span_t)
		+ dz_roundup(tail->cnt.column + q.qpos + 16, 8) * sizeof(uint8_t)	/* +1 for tail '\0' */
	);
	w->aln = (dz_alignment_t *)dz_arena_malloc(mem, aln_size);
	dz_trace_init_aln(w->aln, tail, q);

	/* slice section and path */
	dz_path_span_t *span = (dz_path_span_t *)(w->aln + 1);
	w->span.ptr  = span + tail->cnt.section + 4;
	w->span.base = span + tail->cnt.section + 4;

	uint8_t *path = (uint8_t *)(w->span.base + 2);
	w->path.ptr  = path + tail->cnt.column + q.qpos;
	w->path.base = path + tail->cnt.column + q.qpos;

	/* push tail sentinels */
	uint32_t const id = dz_extract_id(tail);
	dz_trace_push_span(w, id);
	// *w->path.ptr = '\0';		/* make sure readable as C string */
	_mm_storeu_si128((__m128i *)w->path.ptr, _mm_setzero_si128());
	return;
}

static __dz_vectorize
void dz_trace_init_work(dz_trace_work_t *w, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const *tail, dz_max_qpos_t q)
{
	/* save lengths */
	w->rlen  = 1;		/* is unknown before traceback */
	w->query = query;
	w->ie = profile->iev[0];
	w->de = profile->dev[0];

	/* clear counters */
	memset(w->cnt, 0, 4 * sizeof(uint32_t));

	/* load max column pointers */
	#if !defined(DZ_NUCL_DIAG)
		w->pcap = tail->max.cap;
		w->ccap = NULL;
		w->idx  = q.qpos;

	#else
		w->pcap = tail->max.cap;
		w->gcap = NULL;
		w->ccap = NULL;
		w->idx  = q.qpos;

	#endif

	/* compensate adj when save; make sure ccap == NULL */
	w->score = dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx);

	/* init prev_col pointers */
	#if defined(DZ_NUCL_DIAG)
		dz_trace_unwind_core(w, DZ_S_MATRIX);
	#endif
	dz_trace_unwind_core(w, DZ_S_MATRIX);
	return;
}

static __dz_vectorize
void dz_trace_finalize_path(dz_alignment_t *aln, uint8_t *path, uint8_t *base)
{
	/* put head cap */
	path[-1] = '\0';

	/* save pointer and length */
	aln->path = path;
	aln->path_length = base - path;
	return;
}

static __dz_vectorize
void dz_trace_finalize_span(dz_alignment_t *aln, dz_path_span_t *span, dz_path_span_t *base, size_t plen)
{
	span++;				/* adjust root span */
	size_t const slen = base - span;

	/* save pointers */
	aln->span = span;
	aln->span_length = slen;

	/* offset fixup */
	span[0].offset = 0;
	for(size_t i = 1; i < slen; i++) {
		debug("i(%zu), plen(%zu), offset(%u, %u)", i, plen, span[i].offset, (uint32_t)(plen - span[i].offset));
		span[i].offset = plen - span[i].offset;
	}
	span[slen].offset = plen;
	return;
}

static __dz_vectorize
dz_alignment_t const *dz_trace_finalize(dz_trace_work_t *w)
{
	/* check buffer overrun */
	debugblock({
		assert((uintptr_t)(w->aln + 1) <= (uintptr_t)w->span.ptr + 1, "(%p, %p)", w->aln + 1, w->span.ptr);
		assert((uintptr_t)w->span.base <= (uintptr_t)w->path.ptr, "(%p, %p)", w->span.base, w->path.ptr);
	});

	dz_trace_finalize_path(w->aln, w->path.ptr, w->path.base);
	dz_trace_finalize_span(w->aln, w->span.ptr, w->span.base, w->aln->path_length);

	/* save length */
	w->aln->ref_length = w->rlen - 1;

	/* save match / gap counts */
	w->aln->del_count = w->cnt[DZ_E_MATRIX];
	w->aln->ins_count = w->cnt[DZ_F_MATRIX];
	w->aln->mismatch_count = w->cnt[DZ_S_MATRIX];
	w->aln->match_count = w->cnt[DZ_S_MATRIX + 1];
	return(w->aln);
}

static __dz_vectorize
dz_alignment_t const *dz_trace_core(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_getter_t const *getter, dz_state_t const *ff)
{
	/* no object is returned for root */
	if(ff->max.cap == NULL) { return(NULL); }

	/* get pos; uint32_t -> size_t promotion */
	dz_max_qpos_t const q = dz_calc_max_qpos_core(profile, query, ff);	/* vector index, cell index */
	debug("idx(%zu), root(%p)", (size_t)q.qpos, profile->root);

	/* init working buffer */
	dz_trace_work_t w;
	dz_trace_allocate_aln(&w, mem, ff, q);
	dz_trace_init_work(&w, profile, query, ff, q);

	/* core loop */
	while(!dz_is_root(w.pcap)) {
		if(dz_trace_eat_match(&w, profile, getter)) { continue; }
		if(dz_trace_eat_ins(&w)) { continue; }
		if(dz_trace_eat_del(&w)) { continue; }

		dz_trap();
	}

	/* done */
	return(dz_trace_finalize(&w));
}




/*
 * wrapped APIs
 */
#if !(defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0))

dz_static_assert(DZ_TAIL_MARGIN >= sizeof(dz_state_t));

static __dz_vectorize
uint32_t dz_extract_id(dz_state_t const *tail)
{
	dz_unused(tail);
	return(0);
}

#else	/* DZ_WRAPPED_API */


/*
 * metadata handling:
 * reference sequence metadata (will not be exposed outside)
 */
typedef struct dz_meta_s {
	dz_query_t const *query;

	char const *ref;
	size_t rlen;
	uint32_t rdir;
	uint32_t rid;
} dz_meta_t;

dz_static_assert(DZ_TAIL_MARGIN >= sizeof(dz_state_t) + sizeof(dz_meta_t));

#define dz_meta(_p)					( (dz_meta_t *)(_p) )
#define dz_cmeta(_p)				( (dz_meta_t const *)(_p) )


static __dz_vectorize
void dz_save_meta(dz_state_t const *tail, dz_meta_t meta)
{
	/* just save */
	dz_meta_t *p = dz_meta(tail + 1);
	*p = meta;
	return;
}

static __dz_vectorize
dz_meta_t dz_extract_meta(dz_state_t const *tail)
{
	dz_meta_t const *p = dz_cmeta(tail + 1);
	dz_meta_t meta = *p;
	return(meta);
}

static __dz_vectorize
dz_query_t const *dz_extract_query(dz_state_t const *tail)
{
	dz_meta_t const *p = dz_cmeta(tail + 1);
	return(p->query);
}

static __dz_vectorize
uint32_t dz_extract_id(dz_state_t const *tail)
{
	dz_meta_t const *p = dz_cmeta(tail + 1);
	return(p->rid);
}




/*
 * margined sequence
 */
typedef struct {
	uint8_t *ptr;
	size_t len;
} dz_rseq_wrap_t;

static __dz_vectorize
dz_rseq_wrap_t dz_rseq_wrap(dz_profile_t const *profile, uint8_t const *seq, size_t slen)
{
	#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT) || defined(DZ_NUCL_DIAG) || defined(DZ_NUCL_4BIT)
		__m128i const iv = _mm_set1_epi8(rN);
	#else
		__m128i const iv = _mm_set1_epi8(0);
	#endif

	/* head */
	size_t const margin_size = dz_roundup(profile->max_gap_len, 16);
	uint8_t *base = dz_malloc(slen + 2 * margin_size);
	for(size_t i = 0; i < margin_size; i += 16) {
		_mm_storeu_si128((__m128i *)&base[i], iv);
	}

	/* body */
	uint8_t *body = &base[margin_size];
	memcpy(body, seq, slen);

	/*
	for(size_t i = 0; i < slen; i++) {
		fprintf(stderr, "%c", "ACGT"[seq[i]]);
	}
	fprintf(stderr, "\n");
	*/

	/* tail */
	uint8_t *tail = &body[slen];
	for(size_t i = 0; i < margin_size; i += 16) {
		_mm_storeu_si128((__m128i *)&tail[i], iv);
	}

	debug("pack, base(%p), body(%p), tail(%p), slen(%zu)", base, body, tail, slen);

	return((dz_rseq_wrap_t){
		.ptr = body,
		.len = slen
	});
}

static __dz_vectorize
void dz_rseq_destroy(dz_profile_t const *profile, dz_rseq_wrap_t seq)
{
	size_t const margin_size = dz_roundup(profile->max_gap_len, 16);
	return(dz_free(seq.ptr - margin_size));
}




#if defined(DZ_NUCL_ASCII)
// #if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t const pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_ascii,
		.conv     = dz_query_conv_ascii
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}

static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_ascii_t w __attribute__(( aligned(16) ));
	dz_fetcher_init_ascii(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t const fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fetch_next_t)dz_fetch_next_ascii,
		.get_profile = (dz_get_profile_t)dz_get_profile_ascii,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, char const *ref, size_t rlen, uint32_t dir, dz_state_t const *ff)
{
	dz_unused(ref);
	dz_unused(rlen);
	dz_unused(dir);

	dz_getter_t getter = {
		.opaque    = NULL,
		.get_match = (dz_get_match_t)dz_get_match_ascii
	};
	return(dz_trace_core(mem, profile, query, &getter, ff));
}
// #endif



#elif defined(DZ_NUCL_2BIT)
// #if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t const pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_2bit,
		.conv     = dz_query_conv_2bit
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}

static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_2bit_t w __attribute__(( aligned(16) ));
	dz_fetcher_init_2bit(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t const fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fetch_next_t)dz_fetch_next_2bit,
		.get_profile = (dz_get_profile_t)dz_get_profile_2bit,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, char const *ref, size_t rlen, uint32_t dir, dz_state_t const *ff)
{
	dz_unused(ref);
	dz_unused(rlen);
	dz_unused(dir);

	dz_getter_t getter = {
		.opaque    = NULL,
		.get_match = (dz_get_match_t)dz_get_match_2bit
	};
	return(dz_trace_core(mem, profile, query, &getter, ff));
}
// #endif



#elif defined(DZ_NUCL_DIAG)
// #if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t const pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_diag,
		.conv     = dz_query_conv_diag
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}

static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	/* first add margins around reference-side */
	dz_rseq_wrap_t r = dz_rseq_wrap(profile, (uint8_t const *)ref, rlen);

	/* initialize fetcher */
	dz_fetcher_diag_t w __attribute__(( aligned(16) ));
	dz_fetcher_init_diag(&w, (uint8_t const *)r.ptr, r.len, dir);

	dz_fetcher_t const fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fetch_next_t)dz_fetch_next_diag,
		.get_profile = (dz_get_profile_t)dz_get_profile_diag,
		.get_bound   = NULL
	};
	dz_state_t const *tail = dz_extend_core(mem, profile, query, &fetcher, ff, fcnt);

	/* done */
	dz_rseq_destroy(profile, r);
	return(tail);
}

static __dz_vectorize
dz_alignment_t const *dz_trace_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, char const *ref, size_t rlen, uint32_t dir, dz_state_t const *ff)
{
	dz_getter_diag_t g __attribute__(( aligned(16) ));
	dz_getter_init_diag(&g, (uint8_t const *)ref, rlen, dir);

	dz_getter_t getter = {
		.opaque    = (void const *)&g,
		.get_match = (dz_get_match_t)dz_get_match_diag
	};
	return(dz_trace_core(mem, profile, query, &getter, ff));
}
// #endif



#elif defined(DZ_NUCL_4BIT)
// #if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t const pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_4bit,
		.conv     = dz_query_conv_4bit
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}

static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_4bit_t w __attribute__(( aligned(16) ));
	dz_fetcher_init_4bit(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t const fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fetch_next_t)dz_fetch_next_4bit,
		.get_profile = (dz_get_profile_t)dz_get_profile_4bit,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, char const *ref, size_t rlen, uint32_t dir, dz_state_t const *ff)
{
	dz_unused(ref);
	dz_unused(rlen);
	dz_unused(dir);

	dz_getter_t getter = {
		.opaque    = NULL,
		.get_match = (dz_get_match_t)dz_get_match_4bit
	};
	return(dz_trace_core(mem, profile, query, &getter, ff));
}
// #endif



#elif defined(DZ_PROTEIN)
// #if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t const pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_protein,
		.conv     = dz_query_conv_protein
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}

static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_protein_t w __attribute__(( aligned(16) ));
	dz_fetcher_init_protein(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t const fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fetch_next_t)dz_fetch_next_protein,
		.get_profile = (dz_get_profile_t)dz_get_profile_protein,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, char const *ref, size_t rlen, uint32_t dir, dz_state_t const *ff)
{
	dz_unused(ref);
	dz_unused(rlen);
	dz_unused(dir);

	dz_getter_t getter = {
		.opaque    = NULL,
		.get_match = (dz_get_match_t)dz_get_match_protein
	};
	return(dz_trace_core(mem, profile, query, &getter, ff));
}
// #endif



#elif defined(DZ_REALIGNER)
// #if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t const pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_realigner,
		.conv     = dz_query_conv_realigner
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}

static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_realigner_t w __attribute__(( aligned(16) ));
	dz_fetcher_init_realigner(&w, dir, (dz_realigner_node_t const *)ref);

	dz_fetcher_t const fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fetch_next_t)dz_fetch_next_realigner,
		.get_profile = (dz_get_profile_t)dz_get_profile_realigner,
		.get_bound   = (dz_get_bound_t)dz_get_bound_realigner
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, char const *ref, size_t rlen, uint32_t dir, dz_state_t const *ff)
{
	dz_unused(ref);
	dz_unused(rlen);
	dz_unused(dir);

	dz_getter_t getter = {
		.opaque    = NULL,
		.get_match = (dz_get_match_t)dz_get_match_realigner
	};
	return(dz_trace_core(mem, profile, query, &getter, ff));
}
// #endif
#endif




/*
 * API implementations
 */
static __dz_vectorize
dz_t *dz_initx(
	int8_t const *score_matrix,
	uint16_t ins_open,
	uint16_t ins_extend,
	uint16_t del_open,
	uint16_t del_extend,
	uint16_t max_gap_len,
	uint16_t full_length_bonus)
{
	dz_arena_t *mem = dz_arena_create(DZ_MEM_INIT_SIZE);
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* allocate wrapper context */
	dz_t *self = dz_arena_malloc(mem, sizeof(dz_t));
	self->mem = mem;

	/* init profile */
	dz_allocator_t alloc;
	alloc.ctx = NULL;	/* ctx */
	alloc.fp  = (dz_malloc_t)dz_malloc_wrap;	/* malloc */

	dz_score_conf_t conf;
	conf.score_matrix = score_matrix;
	conf.ins_open     = ins_open;
	conf.ins_extend   = ins_extend;
	conf.del_open     = del_open;
	conf.del_extend   = del_extend;
	conf.scan_ref     = 0;
	conf.full_length_bonus = full_length_bonus;
	conf.max_gap_len  = max_gap_len;

	self->profile = dz_init_profile(&alloc, &conf);
	return(self);
}

#ifdef DZ_FULL_LENGTH_BONUS
static __dz_vectorize
dz_t *dz_init(int8_t const *score_matrix, uint16_t gap_open, uint16_t gap_extend, uint64_t max_gap_len, uint16_t full_length_bonus)
{
	return(dz_initx(score_matrix,
		gap_open, gap_extend,
		gap_open, gap_extend,
		max_gap_len, full_length_bonus
	));
}
#else		/* DZ_FULL_LENGTH_BONUS */
static __dz_vectorize
dz_t *dz_init(int8_t const *score_matrix, uint16_t gap_open, uint16_t gap_extend, uint64_t max_gap_len)
{
	return(dz_initx(score_matrix,
		gap_open, gap_extend,
		gap_open, gap_extend,
		max_gap_len, 0
	));
}
#endif		/* DZ_FULL_LENGTH_BONUS */

static __dz_vectorize
void dz_destroy(dz_t *self)
{
	debug("self(%p)", self);
	if(self == NULL) { return; }

	/* free profile */
	dz_destructor_t dtor = {
		.ctx = NULL,
		.fp  = dz_free_wrap
	};
	dz_destroy_profile(&dtor, self->profile);

	/* arena */
	dz_arena_destroy(self->mem);
	return;
}

static __dz_vectorize
void dz_flush(dz_t *self)
{
	dz_arena_flush(self->mem);

	/* we have dz_t after dz_arena_t; make sure the space is reserved */
	self->mem->stack.top += dz_roundup(sizeof(dz_t), DZ_MEM_ALIGN_SIZE);
	return;
}

/* packing query */
static __dz_vectorize
dz_query_t *dz_pack_query(dz_t *self, char const *query, int64_t qlen)
{
	uint32_t const dir = qlen < 0;
	char const *q  = dir ? &query[qlen] : query;
	size_t const l = dir ? -qlen : qlen;
	debug("q(%p), l(%zu), dir(%u)", q, l, dir);

	return(dz_pack_query_intl(self, q, l, dir));	/* internally calls dz_pack_query_core */
}

static __dz_vectorize
dz_query_t *dz_pack_query_forward(dz_t *self, char const *query, size_t qlen)
{
	return(dz_pack_query_intl(self, query, qlen, 0));
}

static __dz_vectorize
dz_query_t *dz_pack_query_reverse(dz_t *self, char const *query, size_t qlen)
{
	return(dz_pack_query_intl(self, query, qlen, 1));
}

/**
 * @fn dz_extend, dz_scan
 * @brief dz_extend for semi-global alignment (head gaps are penalized for the both directions),
 * and dz_scan for searching anchor (head gap is only penalized for the query side).
 *
 * NOTE: DZ_N_AS_UNMATCHING_BASE is not recommended when dz_scan is used
 */
static __dz_vectorize
dz_forefront_t const *dz_extend(
	dz_t *self,
	dz_query_t const *query,
	dz_forefront_t const **forefronts, size_t fcnt,
	char const *ref, int32_t rlen, uint32_t rid)
{
	uint32_t const dir = rlen < 0 ? 1 : 0;		/* 1 for reverse */
	char const *p  = rlen < 0 ? &ref[rlen] : ref;
	size_t const l = rlen < 0 ? -rlen : rlen;	/* always positive */
	// debug("extend, self(%p), fcnt(%zu), query(%p), ref(%p), rlen(%d), rid(%u), p(%p), l(%zu), dir(%u)", self, fcnt, query, ref, rlen, rid, p, l, dir);

	dz_state_t const *tail = dz_extend_intl(self->mem, self->profile, query, dz_cpstate(forefronts), fcnt, p, l, dir);
	if(tail == NULL) { return(NULL); }

	/* save id and direction */
	dz_meta_t const meta = {
		.query = query,

		.ref  = p,
		.rlen = l,
		.rdir = dir,
		.rid  = rid
	};
	dz_save_meta(tail, meta);
	return(dz_cff(tail));
}

static __dz_vectorize
dz_forefront_t const *dz_scan(
	dz_t *self,
	dz_query_t const *query,
	dz_forefront_t const **forefronts, size_t fcnt,
	char const *ref, int32_t rlen, uint32_t rid)
{
	uint32_t const dir = rlen < 0 ? 1 : 0;		/* 1 for reverse */
	char const *p  = rlen < 0 ? &ref[rlen] : ref;
	size_t const l = rlen < 0 ? -rlen : rlen;	/* always positive */
	debug("scan, self(%p), fcnt(%zu), query(%p), ref(%p), rlen(%d), rid(%u), p(%p), l(%zu), dir(%u)", self, fcnt, query, ref, rlen, rid, p, l, dir);

	/* dup profile object */
	uint8_t buf[self->profile->size];
	dz_profile_t *pdup = (dz_profile_t *)buf;
	memcpy(pdup, self->profile, self->profile->size);
	pdup->init = dz_add_ofs(0);

	/* extend with modified profile */
	dz_state_t const *tail = dz_extend_intl(self->mem, pdup, query, dz_cpstate(forefronts), fcnt, p, l, dir);
	if(tail == NULL) { return(NULL); }

	/* save id and direction */
	dz_meta_t const meta = {
		.query = query,

		.ref  = p,
		.rlen = l,
		.rdir = dir,
		.rid  = rid
	};
	dz_save_meta(tail, meta);
	return(dz_cff(tail));
}

static __dz_vectorize
int64_t dz_calc_max_rpos(dz_t *self, dz_forefront_t const *ff)
{
	/* zero if invalid */
	if(ff == NULL) { return(0); }
	dz_state_t const *tail = dz_cstate(ff);
	dz_meta_t const meta = dz_extract_meta(tail);

	/* calculate tail rpos, convert to signed */
	int32_t const rpos = dz_calc_max_rpos_core(self->profile, meta.query, tail);

	/* fold direction in; int32_t -> int64_t promotion (sign expansion) */
	return(meta.rdir ? -rpos : rpos);	/* negate if reverse */
}

static __dz_vectorize
uint64_t dz_calc_max_qpos(dz_t *self, dz_forefront_t const *ff)
{
	dz_state_t const *tail = dz_cstate(ff);
	dz_meta_t const meta = dz_extract_meta(tail);
	dz_max_qpos_t const q = dz_calc_max_qpos_core(self->profile, meta.query, tail);
	return(q.qpos);
}

static __dz_vectorize
uint64_t dz_calc_max_pos(dz_t *self, dz_forefront_t const *ff)
{
	if(ff == NULL) { return(0); }
	dz_state_t const *tail = dz_cstate(ff);
	dz_meta_t const meta = dz_extract_meta(tail);

	dz_max_pos_t const pos = dz_calc_max_pos_core(self->profile, meta.query, tail);

	/* fold direction */
	int64_t const rpos  = pos.rpos;
	int64_t const rdpos = meta.rdir ? -rpos : rpos;
	return((uint64_t)(rdpos<<32) | (uint64_t)pos.qpos);
}

static __dz_vectorize
dz_alignment_t const *dz_trace(
	dz_t *self,
	dz_forefront_t const *ff)
{
	dz_state_t const *tail = dz_cstate(ff);
	dz_meta_t const meta = dz_extract_meta(tail);
	return(dz_trace_intl(self->mem, self->profile, meta.query, meta.ref, meta.rlen, meta.rdir, tail));
}




/*
 * unittests for wrapped APIs
 */
#if defined(UNITTEST) && UNITTEST != 0

/* unittest stuffs */
#if defined(DZ_NUCL_ASCII)
static size_t const dz_unittest_query_length = 560;
static char const *dz_unittest_query =
	"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
	"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
	"TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
	"ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
	"CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
	"GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
	"AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
	"AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT";

static int8_t const dz_unittest_score_matrix[16] = {
	2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2
};

#elif (defined(DZ_NUCL_2BIT) || defined(DZ_NUCL_DIAG))
static size_t const dz_unittest_query_length = 560;
static char const *dz_unittest_query =
	"\x0\x2\x1\x3\x3\x3\x3\x1\x0\x3\x3\x1\x3\x2\x0\x1\x3\x2\x1\x0\x0\x1\x2\x2\x2\x1\x0\x0\x3\x0\x3\x2\x3\x1\x3"
	"\x1\x3\x2\x3\x2\x3\x2\x2\x0\x3\x3\x0\x0\x0\x0\x0\x0\x0\x2\x0\x2\x3\x2\x3\x1\x3\x2\x0\x3\x0\x2\x1\x0\x2\x1"
	"\x3\x3\x1\x3\x2\x0\x0\x1\x3\x2\x2\x3\x3\x0\x1\x1\x3\x2\x1\x1\x2\x3\x2\x0\x2\x3\x0\x0\x0\x3\x3\x0\x0\x0\x0"
	"\x3\x3\x3\x3\x0\x3\x3\x2\x0\x1\x3\x3\x0\x2\x2\x3\x1\x0\x1\x3\x0\x0\x0\x3\x0\x1\x3\x3\x3\x0\x0\x1\x1\x0\x0"
	"\x3\x0\x3\x0\x2\x2\x1\x0\x3\x0\x2\x1\x2\x1\x0\x1\x0\x2\x0\x1\x0\x2\x0\x3\x0\x0\x0\x0\x0\x3\x3\x0\x1\x0\x2"
	"\x0\x2\x3\x0\x1\x0\x1\x0\x0\x1\x0\x3\x1\x1\x0\x3\x2\x0\x0\x0\x1\x2\x1\x0\x3\x3\x0\x2\x1\x0\x1\x1\x0\x1\x1"
	"\x0\x3\x3\x0\x1\x1\x0\x1\x1\x0\x1\x1\x0\x3\x1\x0\x1\x1\x0\x3\x3\x0\x1\x1\x0\x1\x0\x2\x2\x3\x0\x0\x1\x2\x2"
	"\x3\x2\x1\x2\x2\x2\x1\x3\x2\x0\x1\x2\x1\x2\x3\x0\x1\x0\x2\x2\x0\x0\x0\x1\x0\x1\x0\x2\x0\x0\x0\x0\x0\x0\x2"
	"\x1\x1\x1\x2\x1\x0\x1\x1\x3\x2\x0\x1\x0\x2\x3\x2\x1\x2\x2\x2\x1\x3\x3\x3\x3\x3\x3\x3\x3\x3\x1\x2\x0\x1\x1"
	"\x0\x0\x0\x2\x2\x3\x0\x0\x1\x2\x0\x2\x2\x3\x0\x0\x1\x0\x0\x1\x1\x0\x3\x2\x1\x2\x0\x2\x3\x2\x3\x3\x2\x0\x0"
	"\x2\x3\x3\x1\x2\x2\x1\x2\x2\x3\x0\x1\x0\x3\x1\x0\x2\x3\x2\x2\x1\x0\x0\x0\x3\x2\x1\x0\x2\x0\x0\x1\x2\x3\x3"
	"\x3\x3\x1\x3\x2\x1\x2\x3\x2\x3\x3\x2\x1\x1\x2\x0\x3\x0\x3\x3\x1\x3\x2\x2\x0\x0\x0\x2\x1\x0\x0\x3\x2\x1\x1"
	"\x0\x2\x2\x1\x0\x2\x2\x2\x2\x1\x0\x2\x2\x3\x2\x2\x1\x1\x0\x1\x1\x2\x3\x1\x1\x3\x1\x3\x1\x3\x2\x1\x1\x1\x1"
	"\x1\x2\x1\x1\x0\x0\x0\x0\x3\x1\x0\x1\x1\x0\x0\x1\x1\x0\x1\x1\x3\x2\x2\x3\x2\x2\x1\x2\x0\x3\x2\x0\x3\x3\x2"
	"\x0\x0\x0\x0\x0\x0\x1\x1\x0\x3\x3\x0\x2\x1\x2\x2\x1\x1\x0\x2\x2\x0\x3\x2\x1\x3\x3\x3\x0\x1\x1\x1\x0\x0\x3"
	"\x0\x3\x1\x0\x2\x1\x2\x0\x3\x2\x1\x1\x2\x0\x0\x1\x2\x3\x0\x3\x3\x3\x3\x3\x2\x1\x1\x2\x0\x0\x1\x3\x3\x3\x3";

static int8_t const dz_unittest_score_matrix[16] = {
	2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2
};

#elif defined(DZ_NUCL_4BIT)
static size_t const dz_unittest_query_length = 560;
static char const *dz_unittest_query =
	"\x1\x4\x2\x8\x8\x8\x8\x2\x1\x8\x8\x2\x8\x4\x1\x2\x8\x4\x2\x1\x1\x2\x4\x4\x4\x2\x1\x1\x8\x1\x8\x4\x8\x2\x8"
	"\x2\x8\x4\x8\x4\x8\x4\x4\x1\x8\x8\x1\x1\x1\x1\x1\x1\x1\x4\x1\x4\x8\x4\x8\x2\x8\x4\x1\x8\x1\x4\x2\x1\x4\x2"
	"\x8\x8\x2\x8\x4\x1\x1\x2\x8\x4\x4\x8\x8\x1\x2\x2\x8\x4\x2\x2\x4\x8\x4\x1\x4\x8\x1\x1\x1\x8\x8\x1\x1\x1\x1"
	"\x8\x8\x8\x8\x1\x8\x8\x4\x1\x2\x8\x8\x1\x4\x4\x8\x2\x1\x2\x8\x1\x1\x1\x8\x1\x2\x8\x8\x8\x1\x1\x2\x2\x1\x1"
	"\x8\x1\x8\x1\x4\x4\x2\x1\x8\x1\x4\x2\x4\x2\x1\x2\x1\x4\x1\x2\x1\x4\x1\x8\x1\x1\x1\x1\x1\x8\x8\x1\x2\x1\x4"
	"\x1\x4\x8\x1\x2\x1\x2\x1\x1\x2\x1\x8\x2\x2\x1\x8\x4\x1\x1\x1\x2\x4\x2\x1\x8\x8\x1\x4\x2\x1\x2\x2\x1\x2\x2"
	"\x1\x8\x8\x1\x2\x2\x1\x2\x2\x1\x2\x2\x1\x8\x2\x1\x2\x2\x1\x8\x8\x1\x2\x2\x1\x2\x1\x4\x4\x8\x1\x1\x2\x4\x4"
	"\x8\x4\x2\x4\x4\x4\x2\x8\x4\x1\x2\x4\x2\x4\x8\x1\x2\x1\x4\x4\x1\x1\x1\x2\x1\x2\x1\x4\x1\x1\x1\x1\x1\x1\x4"
	"\x2\x2\x2\x4\x2\x1\x2\x2\x8\x4\x1\x2\x1\x4\x8\x4\x2\x4\x4\x4\x2\x8\x8\x8\x8\x8\x8\x8\x8\x8\x2\x4\x1\x2\x2"
	"\x1\x1\x1\x4\x4\x8\x1\x1\x2\x4\x1\x4\x4\x8\x1\x1\x2\x1\x1\x2\x2\x1\x8\x4\x2\x4\x1\x4\x8\x4\x8\x8\x4\x1\x1"
	"\x4\x8\x8\x2\x4\x4\x2\x4\x4\x8\x1\x2\x1\x8\x2\x1\x4\x8\x4\x4\x2\x1\x1\x1\x8\x4\x2\x1\x4\x1\x1\x2\x4\x8\x8"
	"\x8\x8\x2\x8\x4\x2\x4\x8\x4\x8\x8\x4\x2\x2\x4\x1\x8\x1\x8\x8\x2\x8\x4\x4\x1\x1\x1\x4\x2\x1\x1\x8\x4\x2\x2"
	"\x1\x4\x4\x2\x1\x4\x4\x4\x4\x2\x1\x4\x4\x8\x4\x4\x2\x2\x1\x2\x2\x4\x8\x2\x2\x8\x2\x8\x2\x8\x4\x2\x2\x2\x2"
	"\x2\x4\x2\x2\x1\x1\x1\x1\x8\x2\x1\x2\x2\x1\x1\x2\x2\x1\x2\x2\x8\x4\x4\x8\x4\x4\x2\x4\x1\x8\x4\x1\x8\x8\x4"
	"\x1\x1\x1\x1\x1\x1\x2\x2\x1\x8\x8\x1\x4\x2\x4\x4\x2\x2\x1\x4\x4\x1\x8\x4\x2\x8\x8\x8\x1\x2\x2\x2\x1\x1\x8"
	"\x1\x8\x2\x1\x4\x2\x4\x1\x8\x4\x2\x2\x4\x1\x1\x2\x4\x8\x1\x8\x8\x8\x8\x8\x4\x2\x2\x4\x1\x1\x2\x8\x8\x8\x8";

static int8_t const dz_unittest_score_matrix[16 * 16] = {
	2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2
};

#elif defined(DZ_PROTEIN)
static size_t const dz_unittest_query_length = 511;
static char const *dz_unittest_query =
	"MATLVQTGKAKQLTLLGFFAITASMVMAVYEYPTFATSGFSLVFFLLLGGILWFIPVGLC"
	"AAEMATVDGWEEGGVFAWVSNTLGPRWGFAAISFGYLQIAIGFIPMLYFVLGALSYILKW"
	"PALNEDPITKTIAALIILWALALTQFGGTKYTARIAKVGFFAGILLPAFILIALAAIYLH"
	"SGAPVAIEMDSKTFFPDFSKVGTLVVFVAFILSYMGVEASATHVNEMSNPGRDYPLAMLL"
	"LMVAAICLSSVGGLSIAMVIPGNEINLSAGVMQTFTVLMSHVAPEIEWTVRVISALLLLG"
	"VLAEIASWIVGPSRGMYVTAQKNLLPAAFAKMNKNGVPVTLVISQLVITSIALIILTNTG"
	"GGNNMSFLIALALTVVIYLCAYFMLFIGYIVLVLKHPDLKRTFNIPGGKGVKLVVAIVGL"
	"LTSIMAFIVSFLPPDNIQGDSTDMYVELLVVSFLVVLALPFILYAVHDRKGKANTGVTLE"
	"PINSQNAPKGHFFLHPRARSPHYIVMNDKKH";

/* BLOSUM62 */
dz_static_assert(DZ_MAT_SIZE == 32);
static int8_t const dz_unittest_score_matrix[DZ_MAT_SIZE * DZ_MAT_SIZE] = {
	 1, -4, -4, -4, -4, -4, -4, -4, -4, -4,  0, -4, -4, -4, -4,  0, -4, -4, -4, -4, -4,  0, -4, -4, -4, -4, -4,  0,  0,  0,  0,  0,
	-4,  4, -2,  0, -2, -1, -2,  0, -2, -1,  0, -1, -1, -1, -2,  0, -1, -1, -1,  1,  0,  0,  0, -3,  0, -2, -1,  0,  0,  0,  0,  0,
	-4, -2,  4, -3,  4,  1, -3, -1,  0, -3,  0,  0, -4, -3,  3,  0, -2,  0, -1,  0, -1,  0, -3, -4, -1, -3,  1,  0,  0,  0,  0,  0,
	-4,  0, -3,  9, -3, -4, -2, -3, -3, -1,  0, -3, -1, -1, -3,  0, -3, -3, -3, -1, -1,  0, -1, -2, -2, -2, -3,  0,  0,  0,  0,  0,
	-4, -2,  4, -3,  6,  2, -3, -1, -1, -3,  0, -1, -4, -3,  1,  0, -1,  0, -2,  0, -1,  0, -3, -4, -1, -3,  1,  0,  0,  0,  0,  0,
	-4, -1,  1, -4,  2,  5, -3, -2,  0, -3,  0,  1, -3, -2,  0,  0, -1,  2,  0,  0, -1,  0, -2, -3, -1, -2,  4,  0,  0,  0,  0,  0,
	-4, -2, -3, -2, -3, -3,  6, -3, -1,  0,  0, -3,  0,  0, -3,  0, -4, -3, -3, -2, -2,  0, -1,  1, -1,  3, -3,  0,  0,  0,  0,  0,
	-4,  0, -1, -3, -1, -2, -3,  6, -2, -4,  0, -2, -4, -3,  0,  0, -2, -2, -2,  0, -2,  0, -3, -2, -1, -3, -2,  0,  0,  0,  0,  0,
	-4, -2,  0, -3, -1,  0, -1, -2,  8, -3,  0, -1, -3, -2,  1,  0, -2,  0,  0, -1, -2,  0, -3, -2, -1,  2,  0,  0,  0,  0,  0,  0,
	-4, -1, -3, -1, -3, -3,  0, -4, -3,  4,  0, -3,  2,  1, -3,  0, -3, -3, -3, -2, -1,  0,  3, -3, -1, -1, -3,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-4, -1,  0, -3, -1,  1, -3, -2, -1, -3,  0,  5, -2, -1,  0,  0, -1,  1,  2,  0, -1,  0, -2, -3, -1, -2,  1,  0,  0,  0,  0,  0,
	-4, -1, -4, -1, -4, -3,  0, -4, -3,  2,  0, -2,  4,  2, -3,  0, -3, -2, -2, -2, -1,  0,  1, -2, -1, -1, -3,  0,  0,  0,  0,  0,
	-4, -1, -3, -1, -3, -2,  0, -3, -2,  1,  0, -1,  2,  5, -2,  0, -2,  0, -1, -1, -1,  0,  1, -1, -1, -1, -1,  0,  0,  0,  0,  0,
	-4, -2,  3, -3,  1,  0, -3,  0,  1, -3,  0,  0, -3, -2,  6,  0, -2,  0,  0,  1,  0,  0, -3, -4, -1, -2,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-4, -1, -2, -3, -1, -1, -4, -2, -2, -3,  0, -1, -3, -2, -2,  0,  7, -1, -2, -1, -1,  0, -2, -4, -2, -3, -1,  0,  0,  0,  0,  0,
	-4, -1,  0, -3,  0,  2, -3, -2,  0, -3,  0,  1, -2,  0,  0,  0, -1,  5,  1,  0, -1,  0, -2, -2, -1, -1,  3,  0,  0,  0,  0,  0,
	-4, -1, -1, -3, -2,  0, -3, -2,  0, -3,  0,  2, -2, -1,  0,  0, -2,  1,  5, -1, -1,  0, -3, -3, -1, -2,  0,  0,  0,  0,  0,  0,
	-4,  1,  0, -1,  0,  0, -2,  0, -1, -2,  0,  0, -2, -1,  1,  0, -1,  0, -1,  4,  1,  0, -2, -3,  0, -2,  0,  0,  0,  0,  0,  0,
	-4,  0, -1, -1, -1, -1, -2, -2, -2, -1,  0, -1, -1, -1,  0,  0, -1, -1, -1,  1,  5,  0,  0, -2,  0, -2, -1,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-4,  0, -3, -1, -3, -2, -1, -3, -3,  3,  0, -2,  1,  1, -3,  0, -2, -2, -3, -2,  0,  0,  4, -3, -1, -1, -2,  0,  0,  0,  0,  0,
	-4, -3, -4, -2, -4, -3,  1, -2, -2, -3,  0, -3, -2, -1, -4,  0, -4, -2, -3, -3, -2,  0, -3, 11, -2,  2, -3,  0,  0,  0,  0,  0,
	-4,  0, -1, -2, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1,  0, -2, -1, -1,  0,  0,  0, -1, -2, -1, -1, -1,  0,  0,  0,  0,  0,
	-4, -2, -3, -2, -3, -2,  3, -3,  2, -1,  0, -2, -1, -1, -2,  0, -3, -1, -2, -2, -2,  0, -1,  2, -1,  7, -2,  0,  0,  0,  0,  0,
	-4, -1,  1, -3,  1,  4, -3, -2,  0, -3,  0,  1, -3, -1,  0,  0, -1,  3,  0,  0, -1,  0, -2, -3, -1, -2,  4,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

#endif

#ifdef DZ_FULL_LENGTH_BONUS
#define DZ_UNITTEST_SCORE_PARAMS	dz_unittest_score_matrix, 5, 1, 20, 0
#else
#define DZ_UNITTEST_SCORE_PARAMS	dz_unittest_score_matrix, 5, 1, 20
#endif

unittest() {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);
	ut_assert(dz_root(dz) != NULL);
	dz_destroy(dz);
}

unittest() {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_unused(q);
	dz_destroy(dz);
}

/* short, exact matching sequences */
#if !defined(DZ_NUCL_DIAG)
unittest( "extend.base" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 3; trial++) {
		dz_query_t const *q = dz_pack_query(dz, dz_unittest_query,
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			:              10
		);
		dz_forefront_t const *forefront = NULL;

		/* nothing happens */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(0, 0, 0, 0), "max(%u), cell(%u)", forefront->max, dz_ut_sel(0, 0, 0, 0));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 0, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  0);

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("A", "\x0", "\x1", "M"), 1, 4);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 2, 5), "max(%u), cell(%u)", forefront->max, dz_ut_sel(2, 2, 2, 5));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 1, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  1);

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 5);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 4, 9), "max(%u), cell(%u)", forefront->max, dz_ut_sel(4, 4, 4, 9));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 2, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  2);
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "\x1\x4\x1\x8\x8\x8\x8", "MASLVQT"), 7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 9, 28), "max(%u), cell(%u)", forefront->max, dz_ut_sel(9, 9, 9, 28));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 7, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  7);
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "\x1\x4\x1\x8\x8\x8\x8\x2", "MASLVQTG"), 8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 11, 34), "max(%u), cell(%u)", forefront->max, dz_ut_sel(11, 11, 11, 34));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 8, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  8);

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "\x1\x4\x1\x8\x8\x8\x8\x2\x1", "MASLVQTGK"), 9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 13, 39), "max(%u), cell(%u)", forefront->max, dz_ut_sel(13, 13, 13, 39));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 9, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  9);

		dz_unused(forefront);
	}
	dz_destroy(dz);
}

#else	/* DZ_NUCL_DIAG */
unittest( "extend.base" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 3; trial++) {
		dz_query_t const *q = dz_pack_query(dz, dz_unittest_query,
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			:              10
		);
		dz_forefront_t const *forefront = NULL;

		/* nothing happens */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront != NULL && forefront->max == 0, "max(%u), cell(%u)", forefront->max, 0);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0", 1, 4);
		ut_assert(forefront != NULL && forefront->max == 2, "max(%u), cell(%u)", forefront->max, 2);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2", 2, 5);
		ut_assert(forefront != NULL && forefront->max == 4, "max(%u), cell(%u)", forefront->max, 4);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2\x0\x3\x3\x3\x3", 7, 6);
		ut_assert(forefront != NULL && forefront->max == 9, "max(%u), cell(%u)", forefront->max, 9);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2\x0\x3\x3\x3\x3\x1", 8, 7);
		ut_assert(forefront != NULL && forefront->max == 11, "max(%u), cell(%u)", forefront->max, 11);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2\x0\x3\x3\x3\x3\x1\x0", 9, 8);
		ut_assert(forefront != NULL && forefront->max == 13, "max(%u), cell(%u)", forefront->max, 13);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		dz_unused(forefront);
	}
	dz_destroy(dz);
}

#endif

#if !defined(DZ_PROTEIN)
#if !defined(DZ_NUCL_DIAG)
unittest( "extend.base.revcomp" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	char revcomp[dz_unittest_query_length];
	for(size_t i = 0; i < dz_unittest_query_length; i++) {
		switch(dz_unittest_query[dz_unittest_query_length - i - 1]) {
			case 'A': revcomp[i] = 'T'; break;
			case 'C': revcomp[i] = 'G'; break;
			case 'G': revcomp[i] = 'C'; break;
			case 'T': revcomp[i] = 'A'; break;
			case rA: revcomp[i] = rT; break;
			case rC: revcomp[i] = rG; break;
			case rG: revcomp[i] = rC; break;
			case rT: revcomp[i] = rA; break;
		}
	}
	for(size_t trial = 0; trial < 3; trial++) {
		size_t length = trial == 0 ? dz_unittest_query_length
					  : trial == 1 ? 3
					  :              10;
		dz_query_t const *q = dz_pack_query_reverse(dz, &revcomp[dz_unittest_query_length - length], length);
		dz_forefront_t const *forefront = NULL;
#if 0
		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(0, 0, 0, 0), "max(%u), cell(%u)", forefront->max, dz_ut_sel(0, 0, 0, 0));

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("A", "\x0", "\x1", "M"), 1, 4);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 2, 5));
#endif
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 5);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 4, 9));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 2, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  2);
		if(trial == 1) { continue; }


		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "\x1\x4\x1\x8\x8\x8\x8", "MASLVQT"), 7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 9, 28));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 7, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  7);

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"AAAATCT"[7], &"\x0\x0\x0\x0\x3\x1\x3"[7], &"\x1\x1\x1\x1\x8\x2\x8"[7], "MASLVQT"), -7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 9, 28));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 7, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  7);
		if(trial == 2) { continue; }


		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "\x1\x4\x1\x8\x8\x8\x8\x2", "MASLVQTG"), 8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 11, 34));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 8, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  8);

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"GAAAATCT"[8], &"\x2\x0\x0\x0\x0\x3\x1\x3"[8], &"\x4\x1\x1\x1\x1\x8\x2\x8"[8], "MASLVQTG"), -8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 11, 28));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 8, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  8);


		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "\x1\x4\x1\x8\x8\x8\x8\x2\x1", "MASLVQTGK"), 9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 13, 39));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 9, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  9);

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"TGAAAATCT"[9], &"\x3\x2\x0\x0\x0\x0\x3\x1\x3"[9], &"\x8\x4\x1\x1\x1\x1\x8\x2\x8"[9], "MASLVQTGK"), -9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 13, 28));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		ut_assert(forefront != NULL && forefront->cnt.column  == 9, "cnt.column(%u), cell(%u)",  forefront->cnt.column,  9);

		dz_unused(forefront);
	}
	dz_destroy(dz);
}

#else		/* DZ_NUCL_DIAG */
unittest( "extend.base.revcomp" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	char revcomp[dz_unittest_query_length];
	for(size_t i = 0; i < dz_unittest_query_length; i++) {
		switch(dz_unittest_query[dz_unittest_query_length - i - 1]) {
			case 'A': revcomp[i] = 'T'; break;
			case 'C': revcomp[i] = 'G'; break;
			case 'G': revcomp[i] = 'C'; break;
			case 'T': revcomp[i] = 'A'; break;
			case rA: revcomp[i] = rT; break;
			case rC: revcomp[i] = rG; break;
			case rG: revcomp[i] = rC; break;
			case rT: revcomp[i] = rA; break;
		}
	}
	for(size_t trial = 0; trial < 3; trial++) {
		size_t length = trial == 0 ? dz_unittest_query_length
					  : trial == 1 ? 3
					  :              10;
		dz_query_t const *q = dz_pack_query_reverse(dz, &revcomp[dz_unittest_query_length - length], length);
		dz_forefront_t const *forefront = NULL;

		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2", 2, 5);
		ut_assert(forefront != NULL && forefront->max == 4);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		if(trial == 1) { continue; }


		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2\x0\x3\x3\x3\x3", 7, 6);
		ut_assert(forefront != NULL && forefront->max == 9);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		forefront = dz_extend(dz, q, dz_root(dz), 1, &"\x0\x0\x0\x0\x3\x1\x3"[7], -7, 6);
		ut_assert(forefront != NULL && forefront->max == 9);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);
		if(trial == 2) { continue; }


		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2\x0\x3\x3\x3\x3\x1", 8, 7);
		ut_assert(forefront != NULL && forefront->max == 11);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		forefront = dz_extend(dz, q, dz_root(dz), 1, &"\x2\x0\x0\x0\x0\x3\x1\x3"[8], -8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 11, 28));
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);


		forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2\x0\x3\x3\x3\x3\x1\x0", 9, 8);
		ut_assert(forefront != NULL && forefront->max == 13);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		forefront = dz_extend(dz, q, dz_root(dz), 1, &"\x3\x2\x0\x0\x0\x0\x3\x1\x3"[9], -9, 8);
		ut_assert(forefront != NULL && forefront->max == 13);
		ut_assert(forefront != NULL && forefront->cnt.section == 1, "cnt.section(%u), cell(%u)", forefront->cnt.section, 1);

		dz_unused(forefront);
	}
	dz_destroy(dz);
}

#endif		/* DZ_NUCL_DIAG */
#endif		/* DZ_PROTEIN */

#if !defined(DZ_NONGRAPH)		/* skip graphs for diag */
/* a small graph */
unittest( "extend.small" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 4; trial++) {
		dz_query_t const *q = dz_pack_query(dz, dz_unittest_query,
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			: trial == 2 ? 4
			:              7
		);
		dz_forefront_t const *forefronts[5] = { NULL };

		/*
		 * AG---TTTT------CTGA
		 *   \ /    \    /
		 *    C      CATT
		 */
		forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 1);
		ut_assert(forefronts[0] != NULL && forefronts[0]->max == dz_ut_sel(4, 4, 4, 9));
		ut_assert(forefronts[0] != NULL && forefronts[0]->cnt.section == 1);
		ut_assert(forefronts[0] != NULL && forefronts[0]->cnt.column  == 2);

		forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "\x2", "T"), 1, 2);
		ut_assert(forefronts[1] != NULL && forefronts[1]->max == dz_ut_sel(6, 6, 6, 14));
		ut_assert(forefronts[1] != NULL && forefronts[1]->cnt.section == 2);
		ut_assert(forefronts[1] != NULL && forefronts[1]->cnt.column  == 3);
		if(trial == 1) { continue; }

		forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "\x8\x8\x8\x8", "LVQT"), 4, 3);
		if(trial == 2) {
			ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(8, 8, 8, 18));
			continue;
		}
		ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(14, 14, 14, 32));
		ut_assert(forefronts[2] != NULL && forefronts[2]->cnt.section == 3);
		ut_assert(forefronts[2] != NULL && forefronts[2]->cnt.column  == 7);
		if(trial == 3) { continue; }

		forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATT", "\x1\x0\x3\x3", "\x2\x1\x8\x8", "CKAK"), 4, 4);
		ut_assert(forefronts[3] != NULL && forefronts[3]->max == dz_ut_sel(22, 22, 22, 43));
		ut_assert(forefronts[3] != NULL && forefronts[3]->cnt.section == 4);
		ut_assert(forefronts[3] != NULL && forefronts[3]->cnt.column  == 11);

		forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "\x2\x8\x4\x1", "QLTL"), 4, 5);
		ut_assert(forefronts[4] != NULL && forefronts[4]->max == dz_ut_sel(30, 30, 30, 61));
		ut_assert(forefronts[4] != NULL && forefronts[4]->cnt.section == 5);
		ut_assert(forefronts[4] != NULL && forefronts[4]->cnt.column  == 15);
	}
	dz_destroy(dz);
}

unittest( "calc_max.small" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	struct dz_forefront_s const *forefronts[5] = { NULL };

	/*
	 * AG---TTTT------CTGA
	 *   \ /    \    /
	 *    C      CATT
	 */
	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 1);
	ut_assert(forefronts[0] != NULL && forefronts[0]->max == dz_ut_sel(4, 4, 4, 9));
	ut_assert(dz_calc_max_rpos(dz, forefronts[0]) == 2);
	ut_assert(dz_calc_max_qpos(dz, forefronts[0]) == 2);
	ut_assert(dz_calc_max_pos(dz, forefronts[0]) == ((2ULL<<32) | 2ULL));

	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "\x2", "T"), 1, 2);
	ut_assert(forefronts[1] != NULL && forefronts[1]->max == dz_ut_sel(6, 6, 6, 14));
	ut_assert(dz_calc_max_rpos(dz, forefronts[1]) == 1);
	ut_assert(dz_calc_max_qpos(dz, forefronts[1]) == 3);
	ut_assert(dz_calc_max_pos(dz, forefronts[1]) == ((1ULL<<32) | 3ULL));

	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "\x8\x8\x8\x8", "LVQT"), 4, 3);
	ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(14, 14, 14, 32));
	ut_assert(dz_calc_max_rpos(dz, forefronts[2]) == 4);
	ut_assert(dz_calc_max_qpos(dz, forefronts[2]) == 7);
	ut_assert(dz_calc_max_pos(dz, forefronts[2]) == ((4ULL<<32) | 7ULL));

	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATG", "\x1\x0\x3\x2", "\x2\x1\x8\x4", "CKAM"), 4, 4);
	ut_assert(forefronts[3] != NULL && forefronts[3]->max == dz_ut_sel(20, 20, 20, 38), "%d", forefronts[3]->max);
	ut_assert(dz_calc_max_rpos(dz, forefronts[3]) == 3);
	ut_assert(dz_calc_max_qpos(dz, forefronts[3]) == 10);
	ut_assert(dz_calc_max_pos(dz, forefronts[3]) == ((3ULL<<32) | 10ULL));

	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "\x2\x8\x4\x1", "QLTL"), 4, 5);
	ut_assert(forefronts[4] != NULL && forefronts[4]->max == dz_ut_sel(25, 25, 25, 55), "%d", forefronts[4]->max);
	ut_assert(dz_calc_max_rpos(dz, forefronts[4]) == 4);
	ut_assert(dz_calc_max_qpos(dz, forefronts[4]) == 15);
	ut_assert(dz_calc_max_pos(dz, forefronts[4]) == ((4ULL<<32) | 15ULL));

	dz_destroy(dz);
}

unittest( "calc_max.small.revcomp" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	struct dz_forefront_s const *forefronts[5] = { NULL };

	/*
	 * AG---TTTT------CTGA
	 *   \ /    \    /
	 *    C      CATT
	 */
	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"CT"[2], &"\x1\x3"[2], &"\x2\x8"[2], &"AM"[2]), -2, 1);
	ut_assert(dz_calc_max_rpos(dz, forefronts[0]) == -2, "%ld", dz_calc_max_rpos(dz, forefronts[0]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[0]) == 2);
	ut_assert(dz_calc_max_pos(dz, forefronts[0]) == ((((uint64_t)-2LL)<<32) | 2LL));

	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel(&"G"[1], &"\x2"[1], &"\x4"[1], &"T"[1]), -1, 2);
	ut_assert(dz_calc_max_rpos(dz, forefronts[1]) == -1, "%ld", dz_calc_max_rpos(dz, forefronts[1]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[1]) == 3);
	ut_assert(dz_calc_max_pos(dz, forefronts[1]) == ((((uint64_t)-1LL)<<32) | 3LL));

	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel(&"AAAA"[4], &"\x0\x0\x0\x0"[4], &"\x1\x1\x1\x1"[4], &"TQVL"[4]), -4, 3);
	ut_assert(dz_calc_max_rpos(dz, forefronts[2]) == -4, "%ld", dz_calc_max_rpos(dz, forefronts[2]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[2]) == 7);
	ut_assert(dz_calc_max_pos(dz, forefronts[2]) == ((((uint64_t)-4LL)<<32) | 7LL));

	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel(&"CATG"[4], &"\x1\x0\x3\x2"[4], &"\x2\x1\x8\x4"[4], &"MAKC"[4]), -4, 4);
	ut_assert(dz_calc_max_rpos(dz, forefronts[3]) == -3, "%lu", dz_calc_max_rpos(dz, forefronts[3]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[3]) == 10);
	ut_assert(dz_calc_max_pos(dz, forefronts[3]) == ((((uint64_t)-3LL)<<32) | 10LL));

	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel(&"TCAG"[4], &"\x3\x1\x0\x2"[4], &"\x8\x2\x1\x4"[4], &"LTLQ"[4]), -4, 5);
	ut_assert(dz_calc_max_rpos(dz, forefronts[4]) == -4);
	ut_assert(dz_calc_max_qpos(dz, forefronts[4]) == 15);
	ut_assert(dz_calc_max_pos(dz, forefronts[4]) == ((((uint64_t)-4LL)<<32) | 15LL));

	dz_destroy(dz);
}

#endif		/* DZ_NONGRAPH */

unittest( "calc_max.null" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	struct dz_forefront_s const *forefront = NULL;

	forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("TT", "\x3\x3", "\x8\x8", "KC"), 2, 1);
	ut_assert(dz_calc_max_rpos(dz, forefront) == 0);
	ut_assert(dz_calc_max_qpos(dz, forefront) == 0);
	ut_assert(dz_calc_max_pos(dz, forefront) == 0);

	forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"AA"[2], &"\x0\x0"[2], &"\x1\x1"[2], &"KC"[2]), -2, 1);
	ut_assert(dz_calc_max_rpos(dz, forefront) == 0, "%ld", dz_calc_max_rpos(dz, forefront));
	ut_assert(dz_calc_max_qpos(dz, forefront) == 0);
	ut_assert(dz_calc_max_pos(dz, forefront) == 0);

	dz_destroy(dz);
}

/* short, exact matching sequences */
#if !defined(DZ_NONGRAPH)
unittest( "trace" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_forefront_t const *forefronts[5] = { NULL };
	dz_alignment_t const *aln = NULL;

	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 1);
	aln = dz_trace(dz, forefronts[0]);

	ut_assert(aln->score == dz_ut_sel(4, 4, 4, 9));
	ut_assert(aln->ref_length   == 2);
	ut_assert(aln->query_length == 2);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 1);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 2);
	ut_assert(strcmp((char const *)aln->path, "==") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 2);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "\x2", "T"), 1, 2);
	aln = dz_trace(dz, forefronts[1]);

	ut_assert(aln->score == dz_ut_sel(6, 6, 6, 14));
	ut_assert(aln->ref_length   == 3);
	ut_assert(aln->query_length == 3);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 2);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 3);
	ut_assert(strcmp((char const *)aln->path, "===") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 3);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "\x8\x8\x8\x8", "LVQT"), 4, 3);
	aln = dz_trace(dz, forefronts[2]);

	ut_assert(aln->score == dz_ut_sel(14, 14, 14, 32));
	ut_assert(aln->ref_length   == 7);
	ut_assert(aln->query_length == 7);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 3);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 7);
	ut_assert(strcmp((char const *)aln->path, "=======") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 7);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATG", "\x1\x0\x3\x2", "\x2\x1\x8\x4", "GKAM"), 4, 4);
	aln = dz_trace(dz, forefronts[3]);

	ut_assert(aln->score == dz_ut_sel(20, 20, 20, 47));
	ut_assert(aln->ref_length   == 10);
	ut_assert(aln->query_length == 10);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 4);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 10);
	ut_assert(strcmp((char const *)aln->path, "==========") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 10);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "\x2\x8\x4\x1", "QLTL"), 4, 5);
	aln = dz_trace(dz, forefronts[4]);

	ut_assert(aln->score == dz_ut_sel(25, 25, 25, 64));
	ut_assert(aln->ref_length   == 15);
	ut_assert(aln->query_length == 15);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 5);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);
	ut_assert(aln->span[4].id     == 5);
	ut_assert(aln->span[4].offset == 11);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 15);
	ut_assert(strcmp((char const *)aln->path, "==========X====") == 0);

	ut_assert(aln->mismatch_count == 1);
	ut_assert(aln->match_count    == 14);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);

	(void)aln;
	dz_destroy(dz);
}

unittest( "trace.revcomp" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_forefront_t const *forefronts[5] = { NULL };
	dz_alignment_t const *aln = NULL;

	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"CT"[2], &"\x1\x3"[2], &"\x2\x8"[2], &"AM"[2]), -2, 1);
	aln = dz_trace(dz, forefronts[0]);

	ut_assert(aln->score == dz_ut_sel(4, 4, 4, 9));
	ut_assert(aln->ref_length   == 2);
	ut_assert(aln->query_length == 2);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 1);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 2);
	ut_assert(strcmp((char const *)aln->path, "==") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 2);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel(&"G"[1], &"\x2"[1], &"\x4"[1], &"T"[1]), -1, 2);
	aln = dz_trace(dz, forefronts[1]);

	ut_assert(aln->score == dz_ut_sel(6, 6, 6, 14));
	ut_assert(aln->ref_length   == 3);
	ut_assert(aln->query_length == 3);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 2);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 3);
	ut_assert(strcmp((char const *)aln->path, "===") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 3);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel(&"AAAA"[4], &"\x0\x0\x0\x0"[4], &"\x1\x1\x1\x1"[4], &"TQVL"[4]), -4, 3);
	aln = dz_trace(dz, forefronts[2]);

	ut_assert(aln->score == dz_ut_sel(14, 14, 14, 32));
	ut_assert(aln->ref_length   == 7);
	ut_assert(aln->query_length == 7);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 3);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 7);
	ut_assert(strcmp((char const *)aln->path, "=======") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 7);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel(&"CATG"[4], &"\x1\x0\x3\x2"[4], &"\x2\x1\x8\x4"[4], &"MAKG"[4]), -4, 4);
	aln = dz_trace(dz, forefronts[3]);

	ut_assert(aln->score == dz_ut_sel(20, 20, 20, 47));
	ut_assert(aln->ref_length   == 10);
	ut_assert(aln->query_length == 10);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 4);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 10);
	ut_assert(strcmp((char const *)aln->path, "==========") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 10);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel(&"TCAG"[4], &"\x3\x1\x0\x2"[4], &"\x8\x2\x1\x4"[4], &"LTLQ"[4]), -4, 5);
	aln = dz_trace(dz, forefronts[4]);

	ut_assert(aln->score == dz_ut_sel(25, 25, 25, 64));
	ut_assert(aln->ref_length   == 15);
	ut_assert(aln->query_length == 15);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 5);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);
	ut_assert(aln->span[4].id     == 5);
	ut_assert(aln->span[4].offset == 11);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 15);
	ut_assert(strcmp((char const *)aln->path, "==========X====") == 0);

	ut_assert(aln->mismatch_count == 1);
	ut_assert(aln->match_count    == 14);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);

	(void)aln;
	dz_destroy(dz);
}

/* full-length bonus */
#if defined(DZ_FULL_LENGTH_BONUS)
unittest( "trace.bonus" ) {
	struct dz_s *dz = dz_init(dz_unittest_score_matrix, 5, 1, 20, 10);
	ut_assert(dz != NULL);

	char const *query = dz_ut_sel(
		"ACACTTCTAGACTTTACCACTA",
		"\x0\x1\x0\x1\x3\x3\x1\x3\x0\x2\x0\x1\x3\x3\x3\x0\x1\x1\x0\x1\x3\x0",
		"\x1\x2\x1\x2\x8\x8\x2\x8\x1\x4\x1\x2\x8\x8\x8\x1\x2\x2\x1\x2\x8\x1",
		"MSALLILALVGAAVAPLEDDDK"
	);
	struct dz_query_s const *q = dz_pack_query_forward(dz, query, 22);


	struct dz_forefront_s const *ff[10] = { 0 };
	ff[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("ACAC", "\x0\x1\x0\x1", "\x1\x2\x1\x2", "MSAL"), 4, 0);

	ff[1] = dz_extend(dz, q, &ff[0], 1, dz_ut_sel("TTGT", "\x3\x3\x2\x3", "\x8\x8\x4\x8", "LILA"), 4, 1);
	ff[2] = dz_extend(dz, q, &ff[0], 1, dz_ut_sel("ATCC", "\x0\x3\x1\x1", "\x1\x8\x2\x2", "KKLG"), 4, 2);
	ff[3] = dz_extend(dz, q, &ff[1], 2, dz_ut_sel("AGAC", "\x0\x2\x0\x1", "\x1\x4\x1\x2", "LVGA"), 4, 3);

	ff[4] = dz_extend(dz, q, &ff[3], 1, dz_ut_sel("T", "\x3", "\x8", "A"), 1, 4);
	ff[5] = dz_extend(dz, q, &ff[3], 2, dz_ut_sel("TTCTA", "\x3\x3\x1\x3\x0", "\x8\x8\x2\x8\x1", "AVAFP"), 5, 5);

	ff[6] = dz_extend(dz, q, &ff[5], 1, dz_ut_sel("A", "\x0", "\x1", "A"), 1, 6);
	ff[7] = dz_extend(dz, q, &ff[5], 1, dz_ut_sel("C", "\x1", "\x2", "E"), 1, 6);
	ff[8] = dz_extend(dz, q, &ff[5], 1, dz_ut_sel("G", "\x2", "\x4", "M"), 1, 6);
	ff[9] = dz_extend(dz, q, &ff[6], 3, dz_ut_sel("CACGG", "\x1\x0\x1\x2\x2", "\x2\x1\x2\x4\x4", "EDDCL"), 5, 9);

	/* detect max */
	struct dz_forefront_s const *max = NULL;
	for(size_t i = 0; i < 10; i++) {
		if(max == NULL || ff[i]->max > max->max) { max = ff[i]; }
	}
	ut_assert(max == ff[9]);

	/* traceback */
	struct dz_alignment_s const *aln = dz_trace(dz, max);
	ut_assert(aln->ref_length     == dz_ut_sel(23, 23, 23, 23), "%u", aln->ref_length);
	ut_assert(aln->query_length   == dz_ut_sel(22, 22, 22, 22), "%u", aln->query_length);
	ut_assert(aln->score == dz_ut_sel(33, 33, 33, 88), "%d", aln->score);
	ut_assert(strcmp((char const *)aln->path, dz_ut_sel("======X=======I======XX", "======X=======I======XX", "======X=======I======XX", "===============I=X===XX")) == 0, "%s", aln->path);

	ut_assert(aln->path_length    == dz_ut_sel(23, 23, 23, 23), "%u", aln->path_length);
	ut_assert(aln->span_length    == 6, "%u", aln->span_length);

	ut_assert(aln->span[0].id     == 0, "%u", aln->span[0].id);
	ut_assert(aln->span[0].offset == 0, "%u", aln->span[0].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[0].offset], dz_ut_sel("====", "====", "====", "===="), 4) == 0);

	ut_assert(aln->span[1].id     == 1, "%u", aln->span[1].id);
	ut_assert(aln->span[1].offset == 4, "%u", aln->span[1].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[1].offset], dz_ut_sel("==X=", "==X=", "==X=", "===="), 4) == 0);

	ut_assert(aln->span[2].id     == 3, "%u", aln->span[2].id);
	ut_assert(aln->span[2].offset == 8, "%u", aln->span[2].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[2].offset], dz_ut_sel("====", "====", "====", "===="), 4) == 0);

	ut_assert(aln->span[3].id     == 5, "%u", aln->span[3].id);
	ut_assert(aln->span[3].offset == 12, "%u", aln->span[3].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[3].offset], dz_ut_sel("==I==", "==I==", "==I==", "===I="), 5) == 0);

	ut_assert(aln->span[4].id     == 6, "%u", aln->span[4].id);
	ut_assert(aln->span[4].offset == 17, "%u", aln->span[4].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[4].offset], dz_ut_sel("=", "=", "=", "X"), 1) == 0);

	ut_assert(aln->span[5].id     == 9, "%u", aln->span[5].id);
	ut_assert(aln->span[5].offset == 18, "%u", aln->span[5].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[5].offset], dz_ut_sel("===XX", "===XX", "===XX", "===XX"), 5) == 0);

	dz_destroy(dz);
}
#endif		/* DZ_FULL_LENGTH_BONUS */

#else		/* DZ_NONGRAPH */
unittest( "trace" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_forefront_t const *forefront = NULL;
	dz_alignment_t const *aln = NULL;

	forefront = dz_extend(dz, q, dz_root(dz), 1, "\x0\x2\x1\x3\x3\x3\x3\x1\x0\x3\x2\x1\x3\x2\x0", 15, 1);
	aln = dz_trace(dz, forefront);

	ut_assert(aln->score == 25);
	ut_assert(aln->ref_length   == 15, "%u", aln->ref_length);
	ut_assert(aln->query_length == 15, "%u", aln->query_length);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 1);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 15);
	ut_assert(strcmp((char const *)aln->path, "==========X====") == 0);

	ut_assert(aln->mismatch_count == 1);
	ut_assert(aln->match_count    == 14);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);

	(void)aln;
	dz_destroy(dz);
}

unittest( "trace.revcomp" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_forefront_t const *forefront = NULL;
	dz_alignment_t const *aln = NULL;

	forefront = dz_extend(dz, q, dz_root(dz), 1, &"\x3\x1\x0\x2\x1\x0\x3\x2\x0\x0\x0\x0\x2\x1\x3"[15], -15, 1);
	aln = dz_trace(dz, forefront);

	ut_assert(aln->score == 25);
	ut_assert(aln->ref_length   == 15);
	ut_assert(aln->query_length == 15);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 1);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 15);
	ut_assert(strcmp((char const *)aln->path, "==========X====") == 0);

	ut_assert(aln->mismatch_count == 1);
	ut_assert(aln->match_count    == 14);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);

	(void)aln;
	dz_destroy(dz);
}

/* gap counts */
unittest( "trace.gaps" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	char const *query = "\x0\x1\x0\x1\x3\x3\x1\x3\x0\x2\x0\x1\x3\x3\x3\x0\x1\x1\x0\x1\x3\x0\x0\x1\x0\x1\x3\x3";
	struct dz_query_s const *q = dz_pack_query_forward(dz, query, 28);


	debug("start start start");
	struct dz_forefront_s const *ff = NULL;
	ff = dz_extend(dz, q, dz_root(dz), 1, "\x0\x1\x0\x1\x3\x3\x2\x3\x0\x2\x0\x1\x3\x3\x1\x3\x0\x1\x1\x0\x1\x2\x2\x0\x1\x0\x1\x3\x3", 29, 0);

	/* traceback */
	struct dz_alignment_s const *aln = dz_trace(dz, ff);
	ut_assert(aln->ref_length     == 29, "%u", aln->ref_length);
	ut_assert(aln->query_length   == 28, "%u", aln->query_length);
	ut_assert(aln->score == 35, "%d", aln->score);
	ut_assert(strcmp((char const *)aln->path, "======X=======I======XX======") == 0, "%s", aln->path);

	ut_assert(aln->path_length    == 29, "%u", aln->path_length);
	ut_assert(aln->span_length    == 1, "%u", aln->span_length);

	ut_assert(aln->span[0].id     == 0, "%u", aln->span[0].id);
	ut_assert(aln->span[0].offset == 0, "%u", aln->span[0].offset);

	dz_destroy(dz);
}

#endif		/* DZ_NONGRAPH */


static __dz_force_inline
uint8_t dz_unittest_random_char(void)
{
	#if !defined(DZ_PROTEIN)
		switch(rand() % 4) {
			case 0:  return(dz_ut_sel('A', 0, 1, 0));
			case 1:  return(dz_ut_sel('C', 1, 2, 0));
			case 2:  return(dz_ut_sel('G', 2, 4, 0));
			case 3:  return(dz_ut_sel('T', 3, 8, 0));
			default: return(dz_ut_sel('A', 0, 1, 0));
		}
	#else
		return(rand() % 20);
	#endif
}

typedef struct {
	char *seq;
	size_t len;
} dz_unittest_seq_t;

static __dz_force_inline
dz_unittest_seq_t dz_unittest_rand_seq(size_t len)
{
	char *seq = malloc(sizeof(char) * (len + 1));
	for(size_t i = 0; i < len; i++) {
		seq[i] = dz_unittest_random_char();
	}
	seq[len] = '\0';

	return((dz_unittest_seq_t){
		.seq = seq,
		.len = len
	});
}

static __dz_force_inline
dz_unittest_seq_t dz_unittest_mod_seq(char const *seq, size_t len)
{
	char *mod = malloc(sizeof(char) * (4 * len + 1024));
	size_t j = 0;
	for(size_t i = 0; i < len; i++) {
		switch(rand() % 20) {
			case 0: break;	/* del */
			case 1:			/* ins */
				mod[j++] = dz_unittest_random_char();
				mod[j++] = seq[i];
				break;
			// case 2:
				// mod[j++] = dz_unittest_random_char();
				// break;
			default:
				mod[j++] = seq[i];
				break;
		}
	}
	mod[j] = '\0';
	return((dz_unittest_seq_t){
		.seq = mod,
		.len = j
	});
}


unittest( "long" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	size_t const len[3] = { 32, 16384 };
	for(size_t i = 0; i < 2; i++) {
		dz_unittest_seq_t query = dz_unittest_rand_seq(len[i]);
		struct dz_query_s const *q = dz_pack_query_forward(dz, query.seq, query.len);
		ut_assert(q != NULL);

		struct dz_forefront_s const *ff = dz_extend(dz, q, dz_root(dz), 1, query.seq, query.len, 0);
		ut_assert(ff != NULL);
		#if !defined(DZ_PROTEIN)
			/* FIXME: accumulate score for protein */
			ut_assert(ff->max == (int32_t)(dz_unittest_score_matrix[0] * query.len), "%d, %d", ff->max, (int32_t)(dz_unittest_score_matrix[0] * query.len));
		#endif

		struct dz_alignment_s const *aln = dz_trace(dz, ff);
		ut_assert(aln != NULL);

		free(query.seq);
	}
	dz_destroy(dz);
}

unittest( "long.gaps" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	size_t const cnt = 32 * 1024;
	size_t const len = 16384;

	for(size_t i = 0; i < cnt; i++) {
		dz_flush(dz);

		dz_unittest_seq_t query = dz_unittest_rand_seq(len);
		struct dz_query_s const *q = dz_pack_query_forward(dz, query.seq, query.len);
		ut_assert(q != NULL);

		dz_unittest_seq_t ref = dz_unittest_mod_seq(query.seq, query.len);
		struct dz_forefront_s const *ff = dz_extend(dz, q, dz_root(dz), 1, ref.seq, ref.len, 0);
		ut_assert(ff != NULL);
		/* FIXME? we don't exactly know how large the score is */

		uint64_t pos = dz_calc_max_pos(dz, ff);
		ut_assert(pos != 0);
		/* FIXME? we don't exactly know where the maxpos is */

		struct dz_alignment_s const *aln = dz_trace(dz, ff);
		ut_assert(aln != NULL);		/* we only expect the trace back is done without failure */

		free(query.seq);
		free(ref.seq);
	}
	dz_destroy(dz);
}

unittest( "bench" ) {
	#if defined(DZ_FULL_LENGTH_BONUS)
		struct dz_s *dz = dz_init(dz_unittest_score_matrix, 5, 1, 100, 0);
	#else
		struct dz_s *dz = dz_init(dz_unittest_score_matrix, 5, 1, 100);
	#endif
	ut_assert(dz != NULL);

	size_t const cnt = 64 * 1024;
	size_t const len = 4096;

	dz_unittest_seq_t query = dz_unittest_rand_seq(len);
	dz_unittest_seq_t ref = dz_unittest_mod_seq(query.seq, query.len);

	struct dz_query_s const *q = dz_pack_query_forward(dz, query.seq, query.len);
	ut_assert(q != NULL);

	for(size_t i = 0; i < cnt; i++) {
		dz_freeze_t const *fz = dz_arena_freeze(dz->mem);

		struct dz_forefront_s const *ff = dz_extend(dz, q, dz_root(dz), 1, ref.seq, ref.len, 0);
		ut_assert(ff != NULL);

		dz_arena_restore(dz->mem, fz);
	}

	free(query.seq);
	free(ref.seq);
	dz_destroy(dz);
}

static __dz_force_inline
dz_unittest_seq_t dz_unittest_sample_seq(dz_unittest_seq_t *seq, size_t len, size_t cnt)
{
	size_t const rlen = 1 * 1024 * 1024;
	size_t const blen = cnt * len + 2 * rlen;
	char *buf = malloc(sizeof(char) * blen);

	/*
	 * compose graph as following:
	 *
	 *   B---D---F---H     X
	 *  / \ / \ / \ /       \
	 * A   X   X   X   ...   Z
	 *  \ / \ / \ / \       /
	 *   C---E---G---I     Y
	 *
	 * and build a seq by sampling either top or bottom nodes.
	 */
	size_t bpos = 0;
	memcpy(&buf[bpos], seq[0].seq, seq[0].len); bpos += seq[0].len;
	for(size_t i = 1; i < cnt - 1; i += 2) {
		size_t j = rand() & 0x01ULL;
		memcpy(&buf[bpos], seq[i + j].seq, seq[i + j].len); bpos += seq[i + j].len;
	}

	/* clip tail */
	size_t const tlen = seq[cnt - 1].len - dz_min2(seq[cnt - 1].len - 1, 8192);
	memcpy(&buf[bpos], seq[cnt - 1].seq, tlen); bpos += tlen;

	/* add random bases */
	dz_unittest_seq_t rseq = dz_unittest_rand_seq(rlen);
	memcpy(&buf[bpos], rseq.seq, rseq.len); bpos += rseq.len;
	free(rseq.seq);

	/* cap */
	buf[bpos] = '\0';

	#if 0
		/* put gaps and substitutions */
		dz_unittest_seq_t m = dz_unittest_mod_seq(buf, bpos);

		free(buf);
		return(m);
	#else
		return((dz_unittest_seq_t){
			.seq = buf,
			.len = bpos
		});
	#endif
}

#if !defined(DZ_NONGRAPH)
unittest( "long.graph" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	size_t const cnt = 128;
	dz_unittest_seq_t ref[cnt];
	memset(ref, 0, cnt * sizeof(dz_unittest_seq_t));

	size_t const len = 2 * 1024 * 1024;
	for(size_t i = 0; i < cnt; i++) {
		ref[i] = dz_unittest_rand_seq((i % 3) == 0 ? 16 : len);
	}

	dz_unittest_seq_t query = dz_unittest_sample_seq(ref, len, cnt);
	struct dz_query_s const *q = dz_pack_query_forward(dz, query.seq, query.len);
	ut_assert(q != NULL);

	/* root */
	struct dz_forefront_s const *ff[cnt];
	ff[0] = dz_extend(dz, q, dz_root(dz), 1, ref[0].seq, ref[0].len, 0);
	ut_assert(ff[0] != NULL);
	ut_assert(ff[0]->max > 0, "%d", ff[0]->max);		/* score positive */

	/* extend over net graph */
	int32_t prev_max = ff[0]->max;
	for(size_t i = 1; i < cnt - 1; i += 2) {
		struct dz_forefront_s const *incoming[2];
		size_t icnt = 0;
		if(i == 1) {
			incoming[icnt++] = ff[0];
		} else {
			if(ff[i - 2] != NULL) { incoming[icnt++] = ff[i - 2]; }
			if(ff[i - 1] != NULL) { incoming[icnt++] = ff[i - 1]; }
		}

		debug("first");
		ff[i] = dz_extend(dz, q, incoming, icnt, ref[i].seq, ref[i].len, i);
		ut_assert(ff[i] != NULL);

		debug("second");
		ff[i + 1] = dz_extend(dz, q, incoming, icnt, ref[i + 1].seq, ref[i + 1].len, i + 1);
		ut_assert(ff[i + 1] != NULL);

		/* either is unterminated */
		ut_assert(!dz_is_terminated(ff[i]) || !dz_is_terminated(ff[i + 1]));

		/* score increasing */
		int32_t max = prev_max;
		if(!dz_is_terminated(ff[i    ])) { max = dz_max2(max, ff[i    ]->max); }
		if(!dz_is_terminated(ff[i + 1])) { max = dz_max2(max, ff[i + 1]->max); }

		ut_assert(max > prev_max, "i(%zu), max(%d), pair(%d, %d), prev_max(%d)", i, max, ff[i]->max, ff[i + 1]->max, prev_max);
		debugblock({ if(max <= prev_max) { trap(); } });
		prev_max = max;
	}

	/* tail */
	ff[cnt - 1] = dz_extend(dz, q,
		&ff[cnt - 3], 2,
		ref[cnt - 1].seq, ref[cnt - 1].len, cnt - 1
	);
	ut_assert(ff[cnt - 1] != NULL);
	ut_assert(ff[cnt - 1]->max > prev_max);
	ut_assert(ff[cnt - 1]->cap != NULL);
	/* we don't exactly know how large the score is */

	struct dz_alignment_s const *aln = dz_trace(dz, ff[cnt - 1]);
	ut_assert(aln != NULL);		/* we only expect the trace back is done without failure */


	free(query.seq);
	for(size_t i = 0; i < cnt; i++) {
		free(ref[i].seq);
	}
	dz_destroy(dz);
}

#endif	/* DZ_NONGRAPH */



#endif	/* defined(UNITTEST) && UNITTEST != 0 */
#endif	/* DZ_WRAPPED_API */



#ifdef __cplusplus
};	/* extern "C" { */
#endif

#endif		/* _DZ_H_INCLUDED */
/**
 * end of dozeu.h
 */
