
/**
 * @file tinymasker.c
 * @brief command line argument parser and main function
 *
 * @author Hajime Suzuki
 * @license MIT
 */


#ifndef UNITTEST
#  define UNITTEST				( 1 )
#endif
#define UNITTEST_UNIQUE_ID		1


#include "utils/utils.h"		/* include all */
#include "tinymasker.h"
unittest_config( .name = "main" );


#include "dozeu.h"
#include "baln.h"
#include "dbg.h"		/* de Bruijn hash */
#include "index.h"		/* index wrapper */
#include "align.h"
#include "masker.h"
#include "patch.h"


/* return codes */
enum main_error_codes {
	ERROR_INTERNAL = 255,
	ERROR_NO_ARG = 1,

	ERROR_OPEN_IDX = 2,
	ERROR_LOAD_IDX = 3,

	ERROR_OPEN_RSEQ = 4,
	ERROR_OPEN_QSEQ = 5
};

typedef struct {
	/* global args */
	char const *args;
	uint32_t verbose, help;
	size_t nth;

	/* index dump mode if not NULL */
	char const *idxdump;

	/* patch mode if not NULL */
	char const *pafload;
	tm_patch_conf_t patch;

	/* scan-and-mask params */
	char const *profile;
	tm_idx_conf_t fallback;			/* tm_conf_t: tm_idx_conf_t */
	// tm_print_conf_t print;

	/* option parser */
	FILE *log;
	opt_t opt;
} tm_conf_t;

#define tm_conf_assert(_cond, ...) ({ \
	if(!(_cond)) { \
		error("" __VA_ARGS__); \
		return(1); \
	} \
})


static int tm_conf_verbose(tm_conf_t *conf, char const *arg) {
	if(arg == NULL || *arg == '\0') {
		conf->verbose = 1;
	} else if(isdigit(*arg)) {
		conf->verbose = (size_t)mm_atoi(arg, 0);
	} else if(*arg == 'v') {
		tm_conf_verbose(conf, arg + 1);
		conf->verbose++;
	} else {
		return(1);
	}
	return(0);
}
static int tm_conf_threads(tm_conf_t *conf, char const *arg) {
	int64_t const nth = mm_atoi(arg, 0);
	tm_conf_assert(nth < TM_MAX_THREADS, "#threads must be less than %d.", TM_MAX_THREADS);

	conf->nth = nth;
	return(0);
}
static int tm_conf_help(tm_conf_t *conf, char const *arg) {
	_unused(arg);

	conf->verbose++;
	conf->help++;
	return(0);
}

/* index filename */
static int tm_conf_idxdump(tm_conf_t *conf, char const *arg) {
	/* FIXME: check sanity */
	conf->idxdump = mm_strdup(arg);
	return(0);
}
static int tm_conf_profile(tm_conf_t *conf, char const *arg) {
	/* FIXME: check sanity */
	conf->profile = mm_strdup(arg);
	return(0);
}

/* indexing behavior */
static int tm_conf_ign_fail(tm_conf_t *conf, char const *arg) {
	_unused(arg);
	conf->fallback.ignore_overflow = 1;
	return(0);
}

/* indexing params */
static int tm_conf_kmer(tm_conf_t *conf, char const *arg) {
	conf->fallback.kmer = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_window(tm_conf_t *conf, char const *arg) {
	conf->fallback.window = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_ccnt(tm_conf_t *conf, char const *arg) {
	conf->fallback.min_scnt = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_qspan(tm_conf_t *conf, char const *arg) {
	conf->fallback.span_thresh = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_match(tm_conf_t *conf, char const *arg) {
	conf->fallback.match = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_mismatch(tm_conf_t *conf, char const *arg) {
	conf->fallback.mismatch = tm_idx_wrap(mm_atoi(arg, 0));		/* positive number expected */
	return(0);
}
static int tm_conf_gap_open(tm_conf_t *conf, char const *arg) {
	conf->fallback.gap_open = tm_idx_wrap(mm_atoi(arg, 0));		/* positive number expected */
	return(0);
}
static int tm_conf_gap_extend(tm_conf_t *conf, char const *arg) {
	conf->fallback.gap_extend = tm_idx_wrap(mm_atoi(arg, 0));	/* positive number expected */
	return(0);
}
static int tm_conf_max_gap(tm_conf_t *conf, char const *arg) {
	conf->fallback.max_gap_len = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_anc_bonus(tm_conf_t *conf, char const *arg) {
	conf->fallback.full_length_bonus = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_min_score(tm_conf_t *conf, char const *arg) {
	conf->fallback.min_score = tm_idx_wrap(mm_atoi(arg, 0));
	return(0);
}
static int tm_conf_use_raw(tm_conf_t *conf, char const *arg) {
	_unused(arg);
	conf->fallback.use_raw_score = tm_idx_wrap(1);
	return(0);
}
/*
static int tm_conf_flip(tm_conf_t *conf, char const *arg) {
	_unused(conf);
	_unused(arg);
	// conf->print.flip = 1;
	return(0);
}
*/


/* patch */
static int tm_conf_pafload(tm_conf_t *conf, char const *arg) {
	conf->pafload = mm_strdup(arg);
	return(0);
}

static int tm_conf_patch_gap(tm_conf_t *conf, char const *arg) {
	conf->patch.max_gap_len = tm_patch_wrap(mm_atoi(arg, 0));
	return(0);
}

static int tm_conf_patch_margin(tm_conf_t *conf, char const *arg) {
	conf->patch.margin_len = tm_patch_wrap(mm_atoi(arg, 0));
	return(0);
}


static int tm_conf_preset(tm_conf_t *conf, char const *arg)
{
	struct tm_conf_preset_s {
		/* key-value pair */
		char const *key;
		char const *val;

		/* child */
		struct tm_conf_preset_s const *children[6];
	};

	/* preset is defined recursively */
	#define _n(_k, ...)		&((struct tm_conf_preset_s const){ (_k), __VA_ARGS__ })
	struct tm_conf_preset_s const *presets[] = { 0 };		/* FIXME */
	#undef _n

	struct tm_conf_preset_s const *const *q = presets;
	mm_split_foreach(arg, 0, ".:", {	/* traverse preset param tree along with parsing */
		while(*q != NULL && strncmp(p, (*q)->key, l) != 0) { q++; }
		if(*q == NULL) {				/* terminate if not matched, not loaded from file */
			int ret = opt_load_conf(&conf->opt, conf, p);
			if(ret == 0) {
				error("no preset params found for `%.*s'.", (int)l, p);
				return(1);
			}
			break;
		}
		opt_parse_line(&conf->opt, conf, (*q)->val);/* apply recursively */
		q = (*q)->children;				/* visit child nodes */
	});
	return(0);
}


static _force_inline
uint64_t tm_conf_check_sanity(tm_conf_t *conf)
{
	return(opt_ecnt(&conf->opt));
}

#if 0
/* invalidated for now; equivalent is found in tm_idx_fill_default */
static _force_inline
uint64_t tm_conf_restore_default(tm_conf_t *conf)
{
	tm_conf_t defaults = {
		.fallback = {
			.kmer   = 4,
			.window = 32,
			.min_scnt    = 4,
			.skip_filter = 0,
			.match       = 2,
			.mismatch    = 3,
			.gap_open    = 5,
			.gap_extend  = 1,
			.max_gap_len = 16,
			.min_score   = 0
		},
		.print = {
			.flip = 0
		}
	};

	/* we assume all element is sized 64bit */
	uint64_t *q = (uint64_t *)conf;
	uint64_t const *p = (uint64_t const *)&defaults;
	for(size_t i = 0; i < offsetof(tm_conf_t, opt) / sizeof(uint64_t); i++) {
		if(q[i] != 0) { continue; }		/* parameter set */
		if(p[i] == 0) { continue; }		/* default value not found */
		q[i] = p[i];					/* load default */
	}
	return(0);
}
#endif

static _force_inline
void tm_conf_destroy_static(tm_conf_t *conf)
{
	if(conf == NULL) { return; }

	free((void *)conf->args);
	free((void *)conf->idxdump);
	free((void *)conf->profile);

	opt_destroy_static(&conf->opt);
	return;
}

static _force_inline
uint64_t tm_conf_init_static(tm_conf_t *conf, char const *const *argv, FILE *fp)
{
	/* NOTE: all the other members are cleared */
	*conf = (tm_conf_t){
		/* logger */
		.log = stderr,

		#define _c(_x)			( (opt_callback_t)(_x) )
		.opt.t = {
			['\0'] = { 0, NULL },

			['h'] = { OPT_BOOL, _c(tm_conf_help) },
			['v'] = { OPT_OPT,  _c(tm_conf_verbose) },
			['t'] = { OPT_REQ,  _c(tm_conf_threads) },

			/* mapping and output */
			// ['F'] = { OPT_BOOL, _c(tm_conf_flip) },

			/* preset and configuration file */
			['d'] = { OPT_REQ,  _c(tm_conf_idxdump) },
			['x'] = { OPT_REQ,  _c(tm_conf_preset) },
			['z'] = { OPT_REQ,  _c(tm_conf_profile) },
			['L'] = { OPT_BOOL, _c(tm_conf_ign_fail) },

			/* fallback parameters */
			['k'] = { OPT_REQ,  _c(tm_conf_kmer) },
			['w'] = { OPT_REQ,  _c(tm_conf_window) },
			['c'] = { OPT_REQ,  _c(tm_conf_ccnt) },
			['S'] = { OPT_REQ,  _c(tm_conf_qspan) },
			['a'] = { OPT_REQ,  _c(tm_conf_match) },
			['b'] = { OPT_REQ,  _c(tm_conf_mismatch) },
			['p'] = { OPT_REQ,  _c(tm_conf_gap_open) },
			['q'] = { OPT_REQ,  _c(tm_conf_gap_extend) },
			['m'] = { OPT_REQ,  _c(tm_conf_min_score) },
			['R'] = { OPT_BOOL, _c(tm_conf_use_raw) },
			['g'] = { OPT_REQ,  _c(tm_conf_max_gap) },
			['l'] = { OPT_REQ,  _c(tm_conf_anc_bonus) },

			/* patch */
			['F'] = { OPT_REQ,  _c(tm_conf_pafload) },
			['G'] = { OPT_REQ,  _c(tm_conf_patch_gap) },
			['M'] = { OPT_REQ,  _c(tm_conf_patch_margin) }
		}
		#undef _c
	};

	/* parse */
	opt_init_static(&conf->opt, fp);
	if(opt_parse_argv(&conf->opt, conf, argv + 1)) { goto _tm_conf_init_fail; }

	/* make sure the result is correct */
	if(tm_conf_check_sanity(conf)) { goto _tm_conf_init_fail; }

	/* parsed without error */
	conf->args = mm_join(argv, ' ');
	return(0);


_tm_conf_init_fail:;
	opt_destroy_static(&conf->opt);
	return(1);
}

/* determine help and verbose level */
typedef struct {
	FILE *fp;
	uint64_t help, quit;
} tm_conf_outfp_t;

static _force_inline
tm_conf_outfp_t tm_conf_get_outfp(tm_conf_t *conf)
{
	/* use stdout for explicit version (-v) and help (-h) options */
	if(conf->verbose == 1 || conf->help > 0) {
		return((tm_conf_outfp_t){
			.fp = stdout,
			.help = conf->help,
			.quit = 1
		});
	}

	/* restore default verbose level */
	if(conf->verbose == 0) {
		conf->verbose = 1;
	}

	/* implicit help invoked when no input file found; always redirected to stderr */
	return((tm_conf_outfp_t){
		.fp = conf->log,
		.help = opt_parg_cnt(&conf->opt) == 0,
		.quit = 0
	});
}

/* print help */
static _force_inline
void tm_conf_print_help(tm_conf_t const *conf, FILE *lfp)
{
	if(conf->verbose == 0) { return; }

	/* compose temporal profile for default params */
	tm_idx_profile_t *profile = tm_idx_build_profile(&conf->fallback);

	/* get match and mismatch */
	tm_idx_score_t const s = tm_idx_get_score(profile);

	#define _level(_conf) ( (_conf)->verbose + 1 )
	#define _msg_impl(_lv, _fmt, ...) { logger_printf(lfp, _fmt "%s\n", __VA_ARGS__); }
	#define _msg(_lv, ...) { if((_lv) <= _level(conf)) { _msg_impl(_lv, __VA_ARGS__, ""); } }

	_msg(2, "\n"
			"  tinymasker - fast repeat masking tool\n"
			"");
	_msg(2, "Usage:\n"
			"  index construction:\n"
			"    $ tinymasker -t4 -d index.tmi repeats.fa\n"
			"\n"
			"  mask and patch:\n"
			"    $ tinymasker -t4 index.tmi reads.fa > mask.paf\n"
			"    $ tinymasker -t4 -f mask.paf reads.fa > patched_reads.fa\n"
			"   -- or --\n"
			"    $ tinymasker -t4 -f repeats.fa reads.fa > patched_reads.fa\n"
			"");
	_msg(2, "General options:");
	_msg(2, "  -t INT       number of threads [%zu]", conf->nth);
	_msg(2, "  -v [INT]     show version number or set verbose level (when number passed)");
	_msg(2, "");
	_msg(2, "Indexing options:");
	_msg(2, "  -d FILE      dump index to FILE (index construction mode)");
	_msg(3, "  -x STR       load preset params for a specific setting []");
	_msg(3, "  -z STR       custom profile configuration (in TOML format) []");
	_msg(3, "  -L           ignore index construction failure (loose mode)");
	_msg(2, "");
	_msg(2, "Indexing and mapping fallbacks:");
	_msg(2, "  NOTE: ignored for prebuilt indices; only applicable when index is built");
	_msg(2, "  -k INT       k-mer length [%u]", profile->kbits>>1);
	_msg(2, "  -w INT       chaining window size [%u]", profile->chain.window.sep.u);
	_msg(3, "  -c INT       minimum seed count for chain [%u]", profile->chain.min_scnt + 1);
	_msg(3, "  -S INT       minimum q-side span for skipping filter [%u]", profile->filter.span_thresh);
	_msg(2, "  -a INT       match award [%u]", s.match);
	_msg(2, "  -b INT       mismatch penalty [%u]", s.mismatch);
	_msg(2, "  -p INT       gap-open penalty [%u]", profile->extend.giv);
	_msg(2, "  -q INT       gap-extension penalty [%u]", profile->extend.gev);
	_msg(2, "  -m INT       minimum score threshold [%d]", profile->extend.min_score);
	_msg(2, "  -R           use raw score instead of complexity adjusted score [%s]", profile->extend.use_raw ? "yes" : "no");
	_msg(3, "  -g INT       max gap length allowed (X-drop threshold) [%u]", profile->extend.vlim);
	_msg(3, "  -l INT       full-length bonus per end [%u]", profile->extend.bonus);
	_msg(2, "");
	_msg(2, "Patching options:");
	_msg(2, "  -F FILE      load alignments for patching; PAF, FASTA, or tmi index");
	_msg(3, "  -G INT       connect gaps less than this length []");
	_msg(3, "  -M INT       margin around masked regions []");
	_msg(2, "");

	if(_level(conf) <= 2) {
		_msg(1, "Some optinos are hidden. Type `tinymasker -hh' to show all.");
		_msg(1, "");
	}

	#undef _level
	#undef _msg_impl
	#undef _msg

	tm_idx_destroy_profile(profile);
	return;
}

/* return 1 if overridden */
static _force_inline
uint64_t tm_conf_is_overridden(tm_conf_t *conf)
{
	/* we assume all element is sized 64bit */
	uint64_t *q = (uint64_t *)&conf->fallback;
	for(size_t i = 0; i < sizeof(tm_idx_conf_t) / sizeof(uint64_t); i++) {
		if(!tm_idx_is_default(q[i])) { return(1); }
	}
	return(0);
}



/* index construction */

static _force_inline
int main_index_error(tm_conf_t const *conf, int error_code, char const *filename)
{
	_unused(conf);
	switch(error_code) {
	/* argument missing */
	case ERROR_NO_ARG: error("argument is not enough. at least one reference file is required."); break;

	/* opening files */
	case ERROR_OPEN_IDX: error("failed to open index file `%s' in write mode. Please check file path and its permission.", filename); break;
	case ERROR_OPEN_RSEQ: error("failed to build index `%s'.", filename); break;
	}
	return(error_code);
}

static _force_inline
int main_index_intl(tm_conf_t *conf, pg_t *pg, pt_t *pt)
{
	/* iterate over index blocks */
	kv_foreach(void *, opt_pargv(&conf->opt), ({
		debug("p(%s)", *p);

		/* generate index */
		tm_idx_t *mi = tm_idx_gen(&conf->fallback, pt, *p, conf->profile, stderr);
		if(mi == NULL) {
			debug("mi == NULL");
			/* failed to open file */
			return(main_index_error(conf, ERROR_OPEN_RSEQ, conf->idxdump));
		}

		/* dump index */
		size_t size = tm_idx_dump(mi, pg, (write_t)pgwrite);

		/* flush output for next batch */
		pg_flush(pg);

		tm_idx_destroy(mi);
		message(conf->log, "built and dumped index for `%s', on-memory size of this chunk: %.1f MB", (char const *)*p, (double)size / (1024ULL * 1024));
	}));
	return(0);
}

static _force_inline
int main_index(tm_conf_t *conf, pt_t *pt)
{
	if(opt_parg_cnt(&conf->opt) == 0) {
		return(main_index_error(conf, ERROR_NO_ARG, NULL));
	}

	/* add suffix if missing and if not /dev/xxx */
	if(!mm_startswith(conf->idxdump, "/dev") && !mm_endswith(conf->idxdump, ".tmi")) {
		message(conf->log, "index filename does not end with `.tmi' (added).");
		conf->idxdump = mm_append((char *)conf->idxdump, ".tmi");
	}

	/* open file in write mode */
	FILE *fp = fopen(conf->idxdump, "wb");
	if(fp == NULL) {
		/* failed open file (locked?) */
		return(main_index_error(conf, ERROR_OPEN_IDX, conf->idxdump));
	}

	/* initialize compressor */
	pg_t pg;
	pg_init_static(&pg, fp, pt);

	/* index_intl does everything */
	int error_code = main_index_intl(conf, &pg, pt);

	/* done */
	pg_stat_t stat = pg_destroy_static(&pg);
	fclose(fp);
	if(error_code) { return(error_code); }

	message(conf->log, "done. total index size (compressed) on disk: %.1f MB.", (double)stat.out / (1024ULL * 1024));
	return(0);
}


/* scan-and-mask */

static _force_inline
int main_scan_error(tm_conf_t const *conf, int error_code, char const *file)
{
	_unused(conf);
	switch(error_code) {
	/* unknown */
	case ERROR_INTERNAL: error("failed to instanciate alignment context."); break;

	/* in mapping */
	case ERROR_OPEN_QSEQ: error("failed to open sequence file `%s'. Please check file path and format.", file); break;

	/* argument missing */
	case ERROR_NO_ARG: error("argument is not enough. a reference and at least one query file are required."); break;

	/* index loading */
	case ERROR_OPEN_IDX: error("failed to open index file `%s'. Please check file path and permission.", file); break;
	case ERROR_LOAD_IDX: error("failed to load index block from `%s'. Please check file path and version, or rebuild the index.", file); break;

	/* index construction */
	case ERROR_OPEN_RSEQ: error("failed to build index from `%s'.", file); break;
	}
	return(error_code);
}

/* working buffer */
typedef struct {
	/* always available */
	tm_conf_t *conf;
	size_t fcnt;				/* fetched count */

	pt_t *pt;
	// tm_print_t *printer;

	char const *ref;			/* reference filename */
	char const *const *query;	/* query filename */
	size_t qcnt;


	/* for prebuilt index */
	FILE *prebuilt;				/* != NULL if prebuilt index is available */
	pg_t pg;
} main_scan_tbuf_t;

static _force_inline
int main_scan_tbuf_destroy_static(main_scan_tbuf_t *w)
{
	if(w->prebuilt != NULL) {
		fclose(w->prebuilt);
		pg_destroy_static(&w->pg);
	}
	return(0);
}

static _force_inline
int main_scan_tbuf_init_static(main_scan_tbuf_t *w, tm_conf_t *conf, char const *const *parg, size_t pcnt/*, tm_print_t *printer*/, pt_t *pt)
{
	/* save args */
	*w = (main_scan_tbuf_t){
		.conf = conf,
		.fcnt = 0,

		.pt = pt,
		// .printer = printer,

		/* inputs */
		.ref = parg[0],
		.query = &parg[1],
		.qcnt = pcnt - 1
	};

	/* check suffix to determine if the file is prebuilt index */
	if(!mm_endswith(parg[0], ".tmi")) {
		return(0);			/* not a prebuilt index */
	}

	/* open; if fails, it would be invalid path or permission */
	if((w->prebuilt = fopen(parg[0], "rb")) == NULL) {
		return(main_scan_error(conf, ERROR_OPEN_IDX, parg[0]));
	}

	/* done; create threads */
	pg_init_static(&w->pg, w->prebuilt, pt);
	return(0);
}


static _force_inline
int main_scan_idx_gen(main_scan_tbuf_t *w, tm_idx_t **pmi)
{
	if(++w->fcnt > 1) { return(0); }

	/* first call */
	tm_idx_t *mi = tm_idx_gen(&w->conf->fallback, w->pt, w->ref, w->conf->profile, stderr);
	if(mi == NULL) {
		/* failed to open file */
		return(main_scan_error(w->conf, ERROR_OPEN_RSEQ, w->ref));
	}
	message(stderr, "built index for `%s'.", w->ref);
	*pmi = mi;
	return(0);
}

static _force_inline
int main_scan_idx_load(main_scan_tbuf_t *w, tm_idx_t **pmi)
{
	if(pg_eof(&w->pg)) { return(0); }

	/* prebuilt index available; try to fetch next block */
	tm_idx_t *mi = tm_idx_load(&w->pg, (read_t)pgread);
	if(++w->fcnt == 1 && mi == NULL) {
		/* would be broken */
		return(main_scan_error(w->conf, ERROR_LOAD_IDX, w->ref));
	}

	/* NULL for tail */
	if(mi != NULL) { message(stderr, "loaded index block from `%s'.", w->ref); }
	*pmi = mi;
	return(0);
}

/* for each query file */
static _force_inline
int main_scan_foreach_qfile(main_scan_tbuf_t *w, tm_mtscan_t *mt)
{
	for(size_t i = 0; i < w->qcnt; i++) {
		if(tm_mtscan_file(mt, w->query[i])) {
			return(main_scan_error(w->conf, ERROR_OPEN_QSEQ, w->query[i]));
		}
	}
	return(0);
}

/* for each index chunk */
static _force_inline
int main_scan_foreach_idx(main_scan_tbuf_t *w)
{
	while(1) {
		tm_idx_t *mi = NULL;

		int const fetcher_error_code = (w->prebuilt == NULL
			? main_scan_idx_gen(w, &mi)
			: main_scan_idx_load(w, &mi)
		);
		if(fetcher_error_code != 0 || mi == NULL) {
			return(fetcher_error_code);		/* error should be handled inside */
		}

		/* instanciate multithreading context for this index chunk */
		tm_mtscan_t *mt = tm_mtscan_init(mi/*, w->printer*/, w->pt);
		if(mt == NULL) {
			return(main_scan_error(w->conf, ERROR_INTERNAL, NULL));
		}

		/* do the task; for each query file */
		int const scan_error_code = main_scan_foreach_qfile(w, mt);

		/* done for this index chunk */
		tm_mtscan_destroy(mt);
		tm_idx_destroy(mi);
		if(scan_error_code != 0) { return(scan_error_code); }
	}
	return(0);
}

/* for each index filename */
static _force_inline
int main_scan(tm_conf_t *conf/*, tm_print_t *printer*/, pt_t *pt)
{
	char const *const *parg = (char const *const *)opt_parg(&conf->opt);
	size_t const pcnt = opt_parg_cnt(&conf->opt);

	/* error if no argument is given */
	if(pcnt < 2) {
		return(main_scan_error(conf, ERROR_NO_ARG, NULL));
	}

	/* instanciate thread-local working buffers */
	main_scan_tbuf_t w;
	int const init_error_code = main_scan_tbuf_init_static(&w, conf, parg, pcnt/*, printer*/, pt);
	if(init_error_code != 0) {
		return(init_error_code);
	}

	/* print warning if conf is overridden */
	if(w.prebuilt != NULL && tm_conf_is_overridden(conf)) {
		warn("Indexing and mapping parameters are ignored for prebuilt indices.");
	}

	/* we consider the first argument as reference */
	int const scan_error_code = main_scan_foreach_idx(&w);

	/* done */
	main_scan_tbuf_destroy_static(&w);
	return(scan_error_code);
}

#if 0
/* printer */
static _force_inline
int main_scan(tm_conf_t *conf, pt_t *pt)
{
	/* instanciate alignment formatter */
	// tm_print_t printer;
	// tm_print_init_static(&printer, &conf->print, conf->args);

	/* dispatch */
	int const error_code = main_scan_intl(conf/*, &printer*/, pt);

	// tm_print_destory_static(&printer);
	return(error_code);
}
#endif


/* patch */
static _force_inline
int main_patch_error(tm_conf_t const *conf, int error_code, char const *filename)
{
	_unused(conf);
	switch(error_code) {
	/* argument missing */
	case ERROR_NO_ARG: error("argument is not enough. at least one reference file is required."); break;

	/* opening files */
	case ERROR_OPEN_PAF: error("failed to open alignment file `%s'. Please check file path and its format.", filename); break;
	case ERROR_OPEN_QSEQ: error("failed to open sequence file `%s'. Please check file path and format.", filename); break;
	}
	return(error_code);
}

static _force_inline
int main_patch_intl(tm_conf_t *conf, rbread_t *fp, pt_t *pt)
{
	tm_patch_tbuf_t w;
	tm_patch_tbuf_init_static(&w, conf->patch, fp);

	kv_foreach(void *, opt_pargv(&conf->opt)) ({
		if(tm_patch_file(&w, *p, pt) != 0) {
			return(main_patch_error(conf, ERROR_OPEN_QSEQ, *p));
		}
	});

	tm_patch_tbuf_destroy_static(&w);
	return(0);
}

static _force_inline
int main_patch(tm_conf_t *conf, pt_t *pt)
{
	if(opt_parg_cnt(&conf->opt) == 0) {
		return(main_patch_error(conf, ERROR_NO_ARG, NULL));
	}

	/* open file in read-binary mode; (because it might be compressed) */
	rbread_t *fp = rbopen(conf->pafload);
	if(fp == NULL) {
		return(main_patch_error(conf, ERROR_OPEN_PAF, conf->pafload));
	}

	int error_code = main_patch_intl(conf, fp, pt);
	rbclose(fp);
	if(error_code) { return(error_code); }

	message(conf->log, "done.");
	return(0);
}


/* create worker threads for indexing and mapping */
static _force_inline
int main_dispatch(tm_conf_t *conf)
{
	pt_t *pt = pt_init(conf->nth);
	if(pt == NULL) {
		return(main_scan_error(conf, ERROR_INTERNAL, NULL));
	}

	/* dispatch either index or scan */
	int const error_code = (conf->idxdump ? main_index : main_scan)(conf, pt);

	/* done */
	pt_destroy(pt);
	return(error_code);
}


/* entry */
int main(int argc, char *argv[])
{
	int error_code = ERROR_NO_ARG;

	/* unittest dispatcher */
	#if defined(UNITTEST) && UNITTEST != 0
		if(argc > 1 && strcmp(argv[1], "unittest") == 0) {
			return(unittest_main(argc, argv));
		}
	#else
		_unused(argc);
	#endif

	/* instanciate option object */
	logger_init();
	tm_conf_t conf;
	if(tm_conf_init_static(&conf, (char const *const *)argv, stderr)) {
		error("error in parsing arguments. abort.");
		goto _main_final;
	}

	/* always print version */
	tm_conf_outfp_t const out = tm_conf_get_outfp(&conf);
	message(out.fp, "Version: %s (%s), Build: %s", tm_version(), tm_commit(), tm_arch_name());

	/* when -h is passed or no input file is given, print help message */
	if(out.help) {
		if(conf.help > 0) { error_code = 0; }	/* also exit status is 0 (not an error) */
		tm_conf_print_help(&conf, out.fp);		/* we use stdout when invoked by -h option */
	}
	if(out.quit) { goto _main_final; }

	/* dispatch tasks and get return code */
	if((error_code = main_dispatch(&conf)) == 0) {
		message(conf.log, "Command: %s", conf.args);	/* print log when succeeded */
	}
	logger_destroy();

_main_final:;
	tm_conf_destroy_static(&conf);
	return(error_code);
}


/**
 * end of tinymasker.c
 */
