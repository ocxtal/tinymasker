
/**
 * @file mmstring.h
 * @brief string handling
 */

#ifndef _MMSTRING_H_INCLUDED
#define _MMSTRING_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"
#include "log.h"


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
int64_t mm_atoi(char const *arg, uint64_t len)
{
	if(arg == NULL) { return(0); }
	if(len == 0) { len = strlen(arg); }
	for(char const *p = arg, *t = &arg[len]; p < t && *p != '\0'; p++) {
		if(!isdigit(*p)) { error("unparsable number `%.*s'.", (int)len, arg); return(0); }
	}
	return(atoi(arg));
}
static _force_inline
double mm_atof(char const *arg, uint64_t len)
{
	if(arg == NULL) { return(0); }
	if(len == 0) { len = strlen(arg); }
	static char const allowed[16] __attribute__(( aligned(16) )) = {
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '+', '.', ',', 'e', 'E'
	};
	v16i8_t v = _load_v16i8(allowed);
	for(char const *p = arg, *t = &arg[len]; p < t && *p != '\0'; p++) {
		if(((v16i8_masku_t){ .mask = _mask_v16i8(_eq_v16i8(_set_v16i8(*p), v)) }).all == 0) {
			error("unparsable number `%.*s'.", (int)len, arg); return(0);
		}
	}
	return(atof(arg));
}


#endif
/**
 * end of mmstring.h
 */
