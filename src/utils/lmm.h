
/**
 * @fn lmm.h
 * @brief malloc with local context
 */
#ifndef _LMM_H_INCLUDED
#define _LMM_H_INCLUDED


/* include global header *before* we include individual dependencies */
#include "common.h"


/* linear-probing local malloc */
#define LMM_ALIGN_SIZE				( 16 )
#define LMM_HEADER_SIZE				( 16 )
#define LMM_MIN_BASE_SIZE			( 128 )
#define LMM_DEFAULT_BASE_SIZE		( 1024 )

/**
 * @struct lmm_s
 */
struct lmm_s {
	uint8_t need_free;
	uint8_t pad1[7];
	uint32_t head_margin, tot_margin;
	void *ptr, *lim;
};
typedef struct lmm_s lmm_t;

static _force_inline
lmm_t *lmm_init_margin(void *base, size_t base_size, uint32_t head_margin, uint32_t tail_margin)
{
	struct lmm_s *lmm = NULL;
	uint8_t need_free = 0;
	if(base == NULL || base_size < LMM_MIN_BASE_SIZE) {
		base_size = MAX2(base_size, LMM_DEFAULT_BASE_SIZE);
		base = malloc(base_size);
		need_free = 1;
	}

	lmm = (struct lmm_s *)base;
	lmm->need_free = need_free;
	lmm->ptr = (void *)((uintptr_t)base + sizeof(struct lmm_s));
	lmm->lim = (void *)((uintptr_t)base + _cutdown(base_size, LMM_ALIGN_SIZE));
	lmm->head_margin = LMM_HEADER_SIZE + head_margin;
	lmm->tot_margin = LMM_HEADER_SIZE + head_margin + tail_margin;
	return(lmm);
}

static _force_inline
lmm_t *lmm_init(void *base, size_t base_size)
{
	return(lmm_init_margin(base, base_size, 0, 0));
}

static _force_inline
void *lmm_clean(lmm_t *lmm)
{
	if(lmm != NULL && lmm->need_free == 1) {
		free((void *)lmm); return(NULL);
	}
	return((void *)lmm);
}
static _force_inline
void *lmm_reserve_mem(lmm_t *lmm, void *ptr, uint64_t tot_size)
{
	uint64_t *sp = (uint64_t *)ptr;
	*sp = tot_size;
	lmm->ptr = _add_offset(sp, tot_size);
	return(sp);
}

static _force_inline
void *lmm_malloc(lmm_t *lmm, size_t size)
{
	if(lmm == NULL) { return(malloc(size)); }

	size += lmm->tot_margin;
	size = _roundup(size, LMM_ALIGN_SIZE);
	#ifndef LMM_DEBUG
		if((uintptr_t)_add_offset(lmm->ptr, size) < (uintptr_t)lmm->lim) {
			return(_add_offset(lmm_reserve_mem(lmm, lmm->ptr, size), lmm->head_margin));
		}
	#endif

	void *ptr = malloc(size);
	if(ptr == NULL) { return(NULL); }
	return(_add_offset(ptr, lmm->head_margin));
}

static _force_inline
void *lmm_realloc(lmm_t *lmm, void *ptr, size_t size)
{
	if(ptr == NULL) { lmm_malloc(lmm, size); }
	if(lmm == NULL) { return(realloc(ptr, size)); }

	size += lmm->tot_margin;
	ptr = _sub_offset(ptr, lmm->head_margin);
	#ifndef LMM_DEBUG
		/* check if prev mem (ptr) is inside mm */
		if((void *)lmm < ptr && ptr < lmm->lim) {

			uint64_t prev_size = *((uint64_t *)ptr);
			if((uintptr_t)ptr + prev_size == (uintptr_t)lmm->ptr	/* the last block */
			&& (uintptr_t)ptr + size < (uintptr_t)lmm->lim) {		/* and room for expansion */
				return(_add_offset(lmm_reserve_mem(lmm, ptr, size), lmm->head_margin));
			}

			/* no room realloc to outside */
			void *np = malloc(size);
			if(np == NULL) { return(NULL); }
			np = _add_offset(np, lmm->head_margin);
			memcpy(np, ptr, prev_size);
			lmm->ptr = ptr;
			return(np);
		}
	#endif

	/* pass to library realloc */
	void *p = realloc(ptr, size);
	if(p == NULL) { p = ptr; }
	return(_add_offset(p, lmm->head_margin));
}

static _force_inline
void lmm_free(lmm_t *lmm, void *ptr)
{
	if(ptr == NULL) { return; }
	if(lmm == NULL) { free(ptr); }

	ptr = _sub_offset(ptr, lmm->head_margin);
	#ifndef LMM_DEBUG
		if((void *)lmm < ptr && ptr < lmm->lim) {
			/* no need to free */
			uint64_t prev_size = *((uint64_t *)ptr);
			if((uintptr_t)ptr + prev_size == (uintptr_t)lmm->ptr) {
				lmm->ptr = ptr;
			}
			return;
		}
	#endif
	free(ptr);
	return;
}

static inline
char *lmm_strdup(lmm_t *lmm, char const *str)
{
	int64_t len = strlen(str);
	char *s = (char *)lmm_malloc(lmm, len + 1);
	memcpy(s, str, len + 1);
	return(s);
}

#endif /* _LMM_H_INCLUDED */
/**
 * end of lmm.h
 */
