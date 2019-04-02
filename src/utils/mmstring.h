
/**
 * @file mmstring.h
 * @brief string handling
 */

#ifndef _MMSTRING_H_INCLUDED
#define _MMSTRING_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"
#include "log.h"

#include <math.h>

static _force_inline
char *mm_strndup(char const *p, size_t l)
{
	if(!p) { return(NULL); }
	char *q = malloc(sizeof(char) * (l + 1));
	memcpy(q, p, l); q[l] = '\0';
	return(q);
}

static _force_inline
char *mm_strdup(char const *p)
{
	if(!p) { return(NULL); }
	kvec_t(char) b = { 0 };
	while(*p != '\0') { kv_push(char, b, *p); p++; }
	kv_push(char, b, '\0');
	return(kv_ptr(b));
}

static _force_inline
char *mm_join(char const *const *p, char c)
{
	kvec_t(char) b = { 0 };
	while(*p != NULL) {
		char const *q = *p;
		while(*q != '\0') { kv_push(char, b, *q); q++; }
		if(*++p != NULL && c != '\0') {
			kv_push(char, b, c);
		}
	}
	kv_push(char, b, '\0');
	return(kv_ptr(b));
}

static _force_inline
size_t mm_strlen(char const *p)
{
	if(p == NULL) { return(0); }
	return(strlen(p));
}

static _force_inline
size_t mm_strcmp(char const *p, char const *q)
{
	if(q == NULL) { return(p != NULL); }
	if(p == NULL) { return(-1); }
	return(strcmp(p, q));
}

static _force_inline
size_t mm_strncmp(char const *p, char const *q, size_t l)
{
	if(q == NULL) { return(p != NULL); }
	if(p == NULL) { return(-1); }
	return(strncmp(p, q, l));
}

static _force_inline
int mm_startswith(char const *p, char const *prf)
{
	size_t l = strlen(p), r = strlen(prf);
	return(l >= r && strncmp(p, prf, r) == 0);
}

static _force_inline
int mm_endswith(char const *p, char const *suf)
{
	size_t l = strlen(p), r = strlen(suf);
	return(l >= r && strncmp(p + l - r, suf, r) == 0);
}

static _force_inline
char *mm_append(char *p, char const *suf)
{
	size_t l = strlen(p), r = strlen(suf);
	p = realloc(p, sizeof(char) * (l + r + 1));
	memcpy(p + l, suf, r); p[l + r] = '\0';
	return(p);
}

/**
 * @fn mm_atoi, mm_atof
 */
static _force_inline
int64_t mm_atoi(char const *arg, size_t len)
{
	if(arg == NULL) { return(0); }
	if(len == 0) { len = strlen(arg); }

	/* FIXME: use pcmp** instruction for better performance? */
	for(char const *p = arg, *t = &arg[len]; p < t && *p != '\0'; p++) {
		if(_unlikely(!isdigit(*p))) {
			error("unparsable number `%.*s'.", (int)len, arg);
			return(INT64_MIN);
		}
	}
	return(atoi(arg));
}
static _force_inline
double mm_atof(char const *arg, size_t len)
{
	if(arg == NULL) { return(0); }
	if(len == 0) { len = strlen(arg); }
	static char const allowed[16] __attribute__(( aligned(16) )) = {
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '+', '.', ',', 'e', 'E'
	};
	v16i8_t v = _load_v16i8(allowed);

	/* FIXME: use pcmp** instruction for better performance? */
	for(char const *p = arg, *t = &arg[len]; p < t && *p != '\0'; p++) {
		if(((v16i8_masku_t){ .mask = _mask_v16i8(_eq_v16i8(_set_v16i8(*p), v)) }).all == 0) {
			error("unparsable number `%.*s'.", (int)len, arg);
			return(NAN);
		}
	}
	return(atof(arg));
}


typedef struct {
	size_t size, used;
	char *ptr;
} mm_strbin_blk_t;
typedef struct {
	kvec_t(mm_strbin_blk_t) blk;
} mm_strbin_t;

static _force_inline
void mm_strbin_init_static(mm_strbin_t *bin)
{
	kv_clear(bin->blk);
	return;
}

static _force_inline
void mm_strbin_destroy_static(mm_strbin_t *bin)
{
	/* memory arena */
	kv_foreach(mm_strbin_blk_t, bin->blk, { free(p->ptr); });
	kv_destroy(bin->blk);
	return;
}

static _force_inline
mm_strbin_blk_t *mm_strbin_allocate(mm_strbin_t *bin, size_t len)
{
	size_t msize = kv_size(bin->blk);
	mm_strbin_blk_t const *tail = kv_tail(bin->blk);

	if(msize == 0 || tail->size - tail->used < len + 1) {
		size_t const min_blk_size = 4 * 1024;
		size_t blk_size = MAX2(len + 1, min_blk_size);

		char *ptr = malloc(sizeof(char) * blk_size);
		kv_push(mm_strbin_blk_t, bin->blk, ((mm_strbin_blk_t){
			.ptr  = ptr,
			.size = blk_size,
			.used = 0
		}));
	}
	return(kv_tail(bin->blk));
}

static _force_inline
char const *mm_bin_strdup(mm_strbin_t *bin, char const *str, size_t len)
{
	/* calc len if needed */
	if(len == 0) { len = strlen(str); }

	/* allocate memory bin */
	mm_strbin_blk_t *blk = mm_strbin_allocate(bin, len);

	/* slice space from bin */
	char *base = &blk->ptr[blk->used];
	blk->used += len + 1;

	/* copy */
	memcpy(base, str, len);
	base[len] = '\0';
	return(base);
}

static _force_inline
char const *mm_bin_join(mm_strbin_t *bin, char const *const *p, char c)
{
	char *str = mm_join(p, c);
	char const *dup = mm_bin_strdup(bin, str, 0);
	free(str);
	return(dup);
}

static _force_inline
char const *mm_bin_append(mm_strbin_t *bin, char const *p, char const *q)
{
	char const *b[3] = { p, q, NULL };
	char *str = mm_join(b, '\0');
	char const *dup = mm_bin_strdup(bin, str, 0);
	free(str);
	return(dup);
}


#endif
/**
 * end of mmstring.h
 */
