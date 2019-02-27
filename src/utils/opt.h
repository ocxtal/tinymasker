
/**
 * @file opt.h
 * @brief GNU-style option parser (and logger) (no need to maintain large switch-case block)
 */

#ifndef _OPT_H_INCLUDED
#define _OPT_H_INCLUDED


/* include global header *before* we include individual dependencies */
#include "common.h"
#include "log.h"


enum opt_type_t { OPT_BOOL = 1, OPT_REQ = 2, OPT_OPT = 3 };
typedef struct opt_s opt_t;
typedef void (*opt_callback_t)(opt_t *o, void *opaque, char const *optarg);
typedef struct { uint8_t type; opt_callback_t fn; } opt_parser_t;
typedef struct { char *p; size_t size, used; } opt_mem_t;

/**
 * @struct opt_t
 * @brief parameter, logger, and parser container
 */
typedef int (*opt_log_t)(opt_t const *o, char level, char const *func, char const *fmt, ...);
struct opt_s {
	ptr_v parg;								/* positional arguments */
	kvec_t(opt_mem_t) mem;					/* string bin (for positional arguments and opt_strdup) */
	uint32_t ecnt;							/* error counter */
	// uint32_t verbose;
	// double inittime;
	// opt_log_t log;
	void *log;								/* output file pointer */
	opt_parser_t t[256];					/* parser functions */
};

#if 0
/**
 * @fn opt_log_printer
 * @brief 0, 1, 2,... for normal message, 8, 9, 10,... for message with timestamp, 16, 17, 18, ... for without header.
 */
static
int opt_log_printer(
	opt_t const *o,			/* option object */
	char level,				/* 'E' and 'W' for error and warning, 0, 1,... for message */
	char const *func,		/* __func__ must be passed */
	char const *fmt,		/* format string */
	...)
{
	if (level < ' ' && (level & 0x07) > o->verbose) {
		return(0);
	}

	va_list l;
	va_start(l, fmt);

	FILE *fp = (FILE *)o->fp;
	int r = 0;
	if(level >= ' ' || (level & 0x10) == 0) {
		if(level >= ' ' || (level & 0x08) == 0) {
			r += fprintf(fp, "[%c::%s] ", level < ' ' ? 'M' : level, func);
		} else {
			r += fprintf(fp, "[%c::%s::%.3f*%.2f] ",
				level < ' ' ? 'M' : level,					/* 'E' for error */
				func,										/* function name */
				realtime() - o->inittime,					/* realtime */
				cputime() / (realtime() - o->inittime));	/* average cpu usage */
		}
	}
	r += vfprintf(fp, fmt, l);								/* body */
	r += fprintf(fp, "\n");
	va_end(l);
	return(r);
}
#endif

/**
 * @macro oassert, olog
 */
#define oassert(_o, _cond, ...)		{ if(!(_cond)) { error("" __VA_ARGS__); (_o)->ecnt++; } }
#define olog(_o, ...)				{ (_o)->log(_o, __VA_ARGS__); }

/**
 * @macro split_foreach
 * @brief split string into tokens and pass each to _body.
 * p, len, and i are reserved for pointer, length, and #parsed.
 * break / continue can be used in the _body to terminate / skip the current element.
 */
#define split_foreach(_ptr, _len, _delims, _body) ({ \
	char const *_q = (_ptr), *_t = &(_ptr)[(_len) == 0 ? UINT32_MAX : (_len)]; \
	size_t i = 0; \
	v16i8_t _dv = _loadu_v16i8(_delims); \
	uint16_t _m, _mask = 0x02<<_tzc_u64(((v16i8_masku_t){		/* reserve space for '\0' */ \
		.mask = _mask_v16i8(_eq_v16i8(_set_v16i8('\0'), _dv)) \
	}).all); \
	_dv = _bsl_v16i8(_dv, 1); _mask--;						/* push '\0' at the head of the vector */ \
	do { \
		char const *_p = _q; \
		/* test char one by one until dilimiter found */ \
		while(((_m = ((v16i8_masku_t){ .mask = _mask_v16i8(_eq_v16i8(_set_v16i8(*_q), _dv)) }).all) & _mask) == 0) { _q++; } \
		/* delimiter found, pass to _body */ \
		char const *p = _p; \
		uint64_t l = _q++ - _p; \
		if(l > 0) { _body; i++; } \
	} while((_m & 0x01) == 0 && _q < _t); \
	i; \
})

#define opt_inittime(_o)			( (_o)->inittime )
#define opt_pargv(_o)				( (_o)->parg )
#define opt_parg(_o)				( (_o)->parg.a )
#define opt_parg_cnt(_o)			( (_o)->parg.n )
// #define opt_verbose(_o)				( (_o)->verbose )
// #define opt_fp(_o)					( (_o)->fp )
#define opt_ecnt(_o)				( (_o)->ecnt )
#define opt_log(_o)					( (_o)->log )

/**
 * @fn opt_init, opt_destroy
 */
static _force_inline
void opt_init_static(opt_t *o, void *fp)
{
	// o->verbose = 2;
	o->log = fp;

	/* make sure kv_ptr(o->mem) is always available */
	kv_reserve(void *, o->parg, 16);
	return;
}
static _force_inline
void opt_destroy_static(opt_t *o)
{
	/* memory arena */
	kv_foreach(opt_mem_t, o->mem, { free(p->p); });
	kv_destroy(o->mem);

	/* positional arguments */
	kv_destroy(o->parg);
	return;
}
#define opt_init(_fp, _log)			{ opt_t *o = calloc(sizeof(opt_t)); opt_init_static(o, _fp, _log); }
#define opt_destroy(_o)				{ opt_destroy_static(_o); free(_o); }

/**
 * @fn opt_strdup
 */
static _force_inline
opt_mem_t *opt_strdup_allocate_mem(opt_t *o, size_t len)
{
	size_t msize = kv_size(o->mem);
	opt_mem_t const *tail = kv_tail(o->mem);
	if(msize == 0 || tail->size - tail->used < len + 1) {
		size_t const min_blk_size = 4 * 1024;
		size_t blk_size = MAX2(len + 1, min_blk_size);

		char *p = malloc(sizeof(char) * blk_size);
		kv_push(opt_mem_t, o->mem, ((opt_mem_t){ .p = p, .size = blk_size, .used = 0 }));
	}
	return(kv_tail(o->mem));	
}
static _force_inline
char const *opt_strdup(opt_t *o, char const *str, size_t len)
{
	/* calc len if needed */
	if(len == 0) { len = strlen(str); }

	/* allocate memory bin */
	opt_mem_t *mem = opt_strdup_allocate_mem(o, len);

	/* copy */
	char *base = &mem->p[mem->used];
	memcpy(base, str, len);
	base[len] = '\0';
	mem->used += len + 1;
	return(base);
}
static _force_inline
char const *opt_join(opt_t *o, char const *const *p, char c)
{
	char *str = mm_join(p, c);
	char const *dup = opt_strdup(o, str, 0);
	free(str);
	return(dup);
}
static _force_inline
char const *opt_append(opt_t *o, char const *p, char const *q)
{
	char const *b[3] = { p, q, NULL };
	char *str = mm_join(b, '\0');
	char const *dup = opt_strdup(o, str, 0);
	free(str);
	return(dup);
}

/**
 * @fn opt_parse_argv
 * @brief parse (argc, argv)-style option arrays. only argv is required here (MUST be NULL-terminated).
 */
static _force_inline
void opt_push_parg(opt_t *o, char const *arg)
{
	/* save positional argument */
	char const *ptr = opt_strdup(o, arg, 0);
	kv_push(void *, o->parg, (void *)ptr);

	/* push and pop; always keep NULL-terminated */
	kv_push(void *, o->parg, NULL);
	(void)kv_pop(o->parg);
	return;
}
static _force_inline
int opt_parse_argv(opt_t *o, void *opaque, char const *const *argv)
{
	char const *const *p = argv - 1;
	char const *q;
	#define _isarg(_q)	( (_q)[0] != '-' || (_q)[1] == '\0' )
	#define _x(x)		( (size_t)(x) )
	while((q = *++p)) {											/* p is a jagged array, must be NULL-terminated */
		if(_isarg(q)) {	opt_push_parg(o, q); continue; }		/* option starts with '-' and longer than 2 letters, such as "-a" and "-ab" */
		while(o->t[_x(*++q)].type == OPT_BOOL) { o->t[_x(*q)].fn(o, opaque, NULL); }/* eat boolean options */
		if(*q == '\0') { continue; }													/* end of positional argument */
		if(!o->t[_x(*q)].fn) { error("unknown option `-%c'.", *q); continue; }	/* argument option not found */
		char const *r = q[1] ? q+1 : (p[1] && _isarg(p[1]) ? *++p : NULL);				/* if the option ends without argument, inspect the next element in the jagged array (originally placed after space(s)) */
		oassert(o, o->t[_x(*q)].type != OPT_REQ || r, "missing argument for option `-%c'.", *q);
		if(o->t[_x(*q)].type != OPT_REQ || r) { o->t[_x(*q)].fn(o, opaque, r); }	/* option with argument would be found at the tail */
	}
	#undef _isarg
	#undef _x
	return(o->ecnt);
}

/**
 * @fn opt_parse_line
 * @brief parse tab/space-delimited command line arguments
 */
static void opt_parse_line(opt_t *o, void *opaque, char const *arg)
{
	kvec_t(char) str = { 0 };
	kvec_t(char const *) ptr = { 0 };
	split_foreach(arg, 0, " \t\r\n", {
		kv_push(char const *, ptr, (char const *)kv_size(str));
		kv_pushm(char, str, p, l);
		kv_push(char, str, '\0');
	});
	kv_foreach(char const *, ptr, { *p += (ptrdiff_t)str.a; });
	kv_push(char const *, ptr, NULL);
	opt_parse_argv(o, opaque, ptr.a);
	free(str.a); free(ptr.a);
	return;
}

/**
 * @fn opt_load_conf
 */
static int opt_load_conf(opt_t *o, void *opaque, char const *arg)
{
	FILE *fp = fopen(arg, "r");
	if(fp == NULL) { error("failed to find configuration file `%s'.", arg); return(0); }

	kvec_t(char) str = { 0 };
	kv_reserve(char, str, 1024);
	while(1) {
		if((str.n += fread(str.a, sizeof(char), str.m - str.n, fp)) < str.m) { break; }
		kv_reserve(char, str, 2 * str.n);
	}
	fclose(fp); kv_push(char, str, '\0');
	for(uint64_t i = 0; i < str.n; i++) {
		if(str.a[i] == '\n' || str.a[i] == '\t') { str.a[i] = ' '; }
	}

	/* parse */
	opt_parse_line(o, opaque, str.a);

	/* cleanup */
	message(o->log, "loading preset params from `%s': `%s'", arg, str.a);
	free(str.a);
	return(1);
}

#endif
/**
 * end of opt.h
 */

