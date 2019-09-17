
/**
 * @file opt.h
 * @brief GNU-style option parser (and logger) (no need to maintain large switch-case block)
 */

#ifndef _OPT_H_INCLUDED
#define _OPT_H_INCLUDED


/* include global header *before* we include individual dependencies */
#include "common.h"
#include "log.h"
#include "mmstring.h"


enum opt_type_t { OPT_BOOL = 1, OPT_REQ = 2, OPT_OPT = 3 };
typedef struct opt_s opt_t;
typedef int (*opt_callback_t)(void *opaque, char const *optarg);		/* non-zero if fail */
typedef struct { uint8_t type; opt_callback_t fn; } opt_parser_t;
typedef struct { char *p; size_t size, used; } opt_mem_t;

/**
 * @struct opt_t
 * @brief parameter, logger, and parser container
 */
typedef int (*opt_log_t)(opt_t const *o, char level, char const *func, char const *fmt, ...);
struct opt_s {
	mm_bin_t bin;
	ptr_v parg;								/* positional arguments */
	size_t ecnt;							/* error counter */
	void *log;								/* output file pointer */
	opt_parser_t t[256];					/* parser functions */
};

#define opt_inittime(_o)			( (_o)->inittime )
#define opt_pargv(_o)				( (_o)->parg )
#define opt_parg(_o)				( (_o)->parg.a )
#define opt_parg_cnt(_o)			( (_o)->parg.n )
#define opt_ecnt(_o)				( (_o)->ecnt )
#define opt_log(_o)					( (_o)->log )

/**
 * @fn opt_init, opt_destroy
 */
static _force_inline
void opt_init_static(opt_t *opt, void *fp)
{
	opt->log = fp;

	/* make sure kv_ptr(opt->mem) is always available */
	mm_bin_init_static(&opt->bin);
	kv_reserve(void *, opt->parg, 16);
	return;
}
static _force_inline
void opt_destroy_static(opt_t *opt)
{
	/* positional arguments */
	mm_bin_destroy_static(&opt->bin);
	kv_destroy(opt->parg);
	return;
}
#define opt_init(_fp, _log)			{ opt_t *_opt = calloc(sizeof(opt_t)); opt_init_static(_opt, _fp, _log); }
#define opt_destroy(_o)				{ opt_destroy_static(_o); free(_o); }


static _force_inline
uint64_t opt_is_argument(char const *q)
{
	return(q[0] != '-' || q[1] == '\0');
}

static _force_inline
void opt_push_parg(opt_t *opt, char const *arg)
{
	/* save positional argument */
	char const *ptr = mm_bin_strdup(&opt->bin, arg, 0);
	kv_push(void *, opt->parg, (void *)ptr);

	/* push and pop; always keep NULL-terminated */
	kv_push(void *, opt->parg, NULL);
	(void)kv_pop(opt->parg);
	return;
}

/**
 * @fn opt_parse_argv
 * @brief parse (argc, argv)-style option arrays. only argv is required here (MUST be NULL-terminated).
 */
static _force_inline
int opt_parse_argv(opt_t *opt, void *opaque, char const *const *argv)
{
	char const *const *p = argv - 1;
	char const *q = NULL;

	/* p is a jagged array, must be NULL-terminated */
	while((q = *++p) != NULL) {
		/* if the option does not start with '-' and has length, it's positional */
		if(opt_is_argument(q)) {
			opt_push_parg(opt, q);
			continue;
		}

		/* option starts with '-' and longer than 2 letters, such as "-a" and "-ab"; eat boolean options other than the last one */
		while(opt->t[(size_t)*++q].type == OPT_BOOL) {
			if(opt->t[(size_t)*q].fn(opaque, NULL)) { opt->ecnt++; }	/* error if nonzero */
		}
		if(*q == '\0') { continue; }		/* end of positional argument */

		/* skip if argument option not found */
		if(opt->t[(size_t)*q].fn == NULL) {
			error("unknown option `-%c'.", *q);
			continue;
		}

		/* if the option ends without argument, inspect the next element in the jagged array (originally placed after space(s)) */
		char const *r = q[1] != '\0' ? q + 1 : (p[1] && opt_is_argument(p[1]) ? *++p : NULL);
		if(opt->t[(size_t)*q].type == OPT_REQ && r == NULL) {
			error("missing argument for option `-%c'.", *q);
		} else {
			if(opt->t[(size_t)*q].fn(opaque, r)) { opt->ecnt++; }
		}
	}
	return(opt->ecnt);
}

/**
 * @fn opt_parse_line
 * @brief parse tab/space-delimited command line arguments
 */
static void opt_parse_line(opt_t *o, void *opaque, char const *arg)
{
	kvec_t(char) str = { 0 };
	kvec_t(char const *) ptr = { 0 };
	mm_split_foreach(arg, 0, " \t\r\n", {
		kv_push(char const *, ptr, (char const *)kv_size(str));
		kv_pushm(char, str, p, l);
		kv_push(char, str, '\0');
	});
	kv_foreach(char const *, ptr, { *p += (ptrdiff_t)str.a; });
	kv_push(char const *, ptr, NULL);
	opt_parse_argv(o, opaque, ptr.a);

	/* done */
	free(str.a);
	free(ptr.a);
	return;
}

/**
 * @fn opt_load_conf
 */
static int opt_load_conf(opt_t *o, void *opaque, char const *arg)
{
	/* open file */
	FILE *fp = fopen(arg, "r");
	if(fp == NULL) {
		error("failed to find configuration file `%s'.", arg);
		return(0);
	}

	/* file opend successfully; init buffer */
	kvecm_t(char) str;
	kvm_init(str, 0, 32);
	kvm_reserve(char, str, 1024);

	/* dump */
	while(1) {
		if((str.n += fread(str.a, sizeof(char), str.m - str.n, fp)) < str.m) { break; }
		kv_reserve(char, str, 2 * str.n);
	}
	fclose(fp);

	/* escape \n and \t */
	kv_push(char, str, '\0');
	v32i8_t const nv = _set_v32i8('\n'), tv = _set_v32i8('\t'), sv = _set_v32i8(' ');
	for(size_t i = 0; i < str.n; i += 32) {
		v32i8_t const v = _loadu_v32i8(&str.a[i]);
		v32i8_t const nm = _eq_v32i8(nv, v);
		v32i8_t const tm = _eq_v32i8(tv, v);
		v32i8_t const w = _sel_v32i8(_or_v32i8(nm, tm), sv, v);
		_storeu_v32i8(&str.a[i], w);

		// if(str.a[i] == '\n' || str.a[i] == '\t') { str.a[i] = ' '; }
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

