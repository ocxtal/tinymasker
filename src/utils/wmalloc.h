
/**
 * @file wmalloc.h
 * @brief malloc wrapper
 */

#ifndef _WMALLOC_H_INCLUDED
#define _WMALLOC_H_INCLUDED


/* include global header *before* we include individual dependencies */
#include "common.h"
#include "debug.h"

#include <sys/resource.h>		/* getrusage */
#ifdef __linux__
#  include <sys/mman.h>			/* madvise */
#endif


/* maxrss coef */
#if defined(__linux__)
#  define WM_MAXRSS_COEF	( 1000.0 )
#elif defined(__APPLE__)
#  define WM_MAXRSS_COEF	( 1000.0 * 1000.0 )
#else
#  define WM_MAXRSS_COEF	( 1000.0 )
#endif


/* madvise configurations */
#define WM_PAGE_SIZE					ARCH_PAGE_SIZE
#define WM_HUGEPAGE_SIZE				ARCH_HUGEPAGE_SIZE
#define WM_ROUNDUP_THRESH				( 128ULL * 1024 * 1024 )


/**
 * @fn oom_abort
 * @brief called when malloc / realloc failed. dump thread information and exit (noreturn).
 */
static
void oom_abort(
	char const *name,
	size_t req)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	fprintf(stderr, "[E::%s] Out of memory. (required: %.1f MB, maxrss: %.1f MB)\n", name, (double)req / (1024.0 * 1024.0), (double)r.ru_maxrss / WM_MAXRSS_COEF);
	trap();								/* segv trap for debugging; see log.h */
	exit(128);							/* 128 reserved for out of memory */
}


/**
 * @macro wrapped_malloc
 */
#define wrapped_malloc(_x) ({ \
	assert((_x) < 20ULL * 1024 * 1024 * 1024, "size(%zu)", (size_t)(_x)); \
	size_t _req_size = (size_t)(_x); \
	void *_ptr = malloc(_req_size); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, _req_size); \
	} \
	_ptr; \
})
#define wrapped_new(type_t, _content) ({ \
	type_t *_ptr = wrapped_malloc(sizeof(type_t)); \
	*_ptr = (type_t)_content; \
	_ptr; \
})
#define malloc(_x)		wrapped_malloc(_x)

/**
 * @macro wrapped_realloc
 */
#define wrapped_realloc(_x, _y) ({ \
	assert((_y) < 20ULL * 1024 * 1024 * 1024, "size(%zu)", (size_t)(_y)); \
	void *_prev_ptr = (void *)(_x); \
	size_t _req_size = (size_t)(_y); \
	void *_ptr = realloc(_prev_ptr, _req_size); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, _req_size); \
	} \
	_ptr; \
})
#define realloc(_x, _y)	wrapped_realloc(_x, _y)

/**
 * @macro wrapped_calloc
 */
#define wrapped_calloc(_x, _y) ({ \
	size_t _req_cnt = (size_t)(_x), _req_size = (size_t)(_y); \
	void *_ptr = calloc(_req_cnt, _req_size); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, (_req_cnt) * (_req_size)); \
	} \
	_ptr; \
})
#define calloc(_x, _y) wrapped_calloc(_x, _y)



/* aligned_malloc; madvise to use transparent hugepage */
typedef struct {
	size_t adjusted_size, align_size;
	uint64_t use_madvise;
} wm_aligned_malloc_hint_t;

static _force_inline
wm_aligned_malloc_hint_t wm_calc_size(size_t size)
{
	size_t const small_align_size = 64;

	#ifdef __linux__
		size_t const hugepage_size = WM_HUGEPAGE_SIZE;
		size_t const roundup_thresh = WM_ROUNDUP_THRESH;

		/* check if size is multiple of 2M */
		uint64_t use_madvise = (size & (uint64_t)(hugepage_size - 1)) == 0;
		return((wm_aligned_malloc_hint_t){
			.adjusted_size = size < roundup_thresh ? size : _roundup(size, hugepage_size),
			.align_size    = use_madvise ? hugepage_size : small_align_size,
			.use_madvise   = use_madvise
		});
	#else
		return((wm_aligned_malloc_hint_t){
			.adjusted_size = size,
			.align_size    = small_align_size,
			.use_madvise   = 0
		});
	#endif
}

static _force_inline
void wm_madvice(wm_aligned_malloc_hint_t hint, void *ptr, char const *fn)
{
	#ifdef __linux__
		if(hint.use_madvise != 0) {
			int ret = posix_madvise(ptr, hint.adjusted_size, MADV_HUGEPAGE);

			if(_unlikely(ret != 0)) {
				fprintf(stderr, "madvise failed for ptr(%p), size(%zu, %zx), ret(%d), in %s\n", ptr, hint.adjusted_size, hint.adjusted_size, ret, fn);
			}
		}
	#else
		_unused(hint);
		_unused(ptr);
		_unused(fn);
	#endif
	return;
}

/**
 * @macro wrapped_aligned_malloc
 */
#define wrapped_aligned_malloc(_x) ({ \
	wm_aligned_malloc_hint_t _hint = wm_calc_size((size_t)(_x)); \
	void *_ptr; \
	if(_unlikely(posix_memalign(&_ptr, _hint.align_size, _hint.adjusted_size) != 0)) { \
		oom_abort(__func__, _hint.adjusted_size); \
	} \
	wm_madvice(_hint, _ptr, __func__); \
	_ptr; \
})
#define aligned_malloc(_x)				wrapped_aligned_malloc(_x)

/*
#define wrapped_posix_memalign(_ptr, _align, _size) ({ \
	void *_p; \
	if(_unlikely(posix_memalign(&_p, (size_t)(_align), (size_t)(_size)) != 0)) { \
		oom_abort(__func__, (_size)); \
	} \
	*(_ptr) = _p; \
	0; \
})
#define posix_memalign(_x, _y, _z)		wrapped_posix_memalign(_x, _y, _z)
*/

#if 0
unittest( .name = "wrapped_malloc" ) {
	uint64_t const size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size);
	ut_assert(p != NULL);

	memset(p, 0, size);							/* make sure we can touch this area */

	p = realloc(p, 2*size);
	ut_assert(p != NULL);

	p = realloc(p, size/2);
	ut_assert(p != NULL);
	free(p);
}
#endif


#endif
/**
 * end of wmalloc.h
 */
