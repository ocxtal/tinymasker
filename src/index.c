
/**
 * @file index.c
 * @brief index data structure construction APIs; wrapper of de Bruijn hash
 *
 * @author Hajime Suzuki
 * @license MIT
 */


#ifndef UNITTEST
#  define UNITTEST				( 1 )
#endif
#define UNITTEST_UNIQUE_ID		4


#include "utils/utils.h"		/* include all */
#include "tinymasker.h"
unittest_config( .name = "index" );


#include "dozeu.h"
#include "dbg.h"		/* de Bruijn hash */
#include "bseq.h"		/* FASTA/Q parser */
#include "toml.h"		/* toml parser */
#include "re.h"			/* regex */


/* forward declarations */
#include "index.h"


static _force_inline
void *tm_idx_malloc(void *unused, size_t size)
{
	_unused(unused);
	return(malloc(size));
}

static _force_inline
void tm_idx_free(void *unused, void *ptr)
{
	_unused(unused);

	free(ptr);
	return;
}

// static _force_inline
void tm_idx_destroy_profile(tm_idx_profile_t *profile)
{
	dz_destructor_t f = {
		.ctx = NULL,
		.fp  = (dz_free_t)tm_idx_free
	};
	dz_destroy_profile(&f, profile->extend.dz);
	free(profile);
	return;
}

static _force_inline
void tm_idx_destroy_parr(tm_idx_profile_t **arr, size_t cnt)
{
	for(size_t i = 0; i < cnt; i++) {
		tm_idx_destroy_profile(arr[i]);
		arr[i] = NULL;
	}
	return;
}

static _force_inline
void tm_idx_destroy_sketch(tm_idx_sketch_t **arr, size_t cnt)
{
	for(size_t i = 0; i < cnt; i++) {
		tm_ref_destroy(&arr[i]->ref);
	}
	return;
}

static _force_inline
void tm_idx_destroy_mono(tm_idx_t *mi)
{
	dz_destructor_t f = {
		.ctx = NULL,
		.fp  = (dz_free_t)tm_idx_free
	};

	/* still we need to destroy profile */
	for(size_t i = 0; i < mi->profile.cnt; i++) {
		dz_destroy_profile(&f, mi->profile.arr[i]->extend.dz);
	}
	free(mi->base);
	return;
}

static _force_inline
void tm_idx_destroy_normal(tm_idx_t *mi)
{
	if(mi->filename.ref) {
		free((void *)mi->filename.ref);
	}
	if(mi->filename.profile) {
		free((void *)mi->filename.profile);
	}

	tm_idx_destroy_parr(mi->profile.arr, mi->profile.cnt);
	tm_idx_destroy_sketch(mi->sketch.arr, mi->sketch.cnt);

	free(mi->profile.arr);
	free(mi->sketch.arr);

	free(mi);
	return;
}

// static _force_inline
void tm_idx_destroy(tm_idx_t *mi)
{
	if(mi == NULL) { return; }

	if(mi->base != NULL) {
		/* monolithic index loaded by tm_idx_load */
		tm_idx_destroy_mono(mi);
		return;
	}

	/* built by tm_idx_gen; not monolithic */
	tm_idx_destroy_normal(mi);
	return;
}


/* index builder */
typedef struct {
	re_pattern_t *name;
	re_pattern_t *comment;
} tm_idx_matcher_t;

typedef struct {
	tm_idx_t mi;

	struct {
		kvec_t(tm_idx_profile_t *) profile;
		kvec_t(tm_idx_matcher_t) matcher;
	} score;

	struct {
		bseq_file_t *fp;
		kvec_t(tm_idx_sketch_t *) bin;
	} col;

	/* working buffers */
	tm_ref_tbuf_t *ref[];
} tm_idx_gen_t;

typedef struct {
	size_t base_rid, _pad[3];
	bseq_batch_t bin;
} tm_idx_batch_t;


static _force_inline
void tm_idx_push_profile(tm_idx_gen_t *mii, tm_idx_matcher_t matcher, tm_idx_profile_t *profile)
{
	kv_push(tm_idx_matcher_t, mii->score.matcher, matcher);
	kv_push(tm_idx_profile_t *, mii->score.profile, profile);
	return;
}

static _force_inline
size_t tm_idx_find_profile(tm_idx_gen_t const *mii, char const *name, size_t nlen, char const *comment, size_t clen)
{
	_unused(nlen);
	_unused(clen);

	tm_idx_matcher_t const *matcher = kv_ptr(mii->score.matcher);
	size_t const mcnt = kv_cnt(mii->score.matcher);

	/* regex match; find first */
	for(size_t i = 0; i < mcnt; i++) {
		tm_idx_matcher_t const *m = &matcher[i];

		/* result buffers */
		re_result_t nres = { NULL, NULL }, cres = { NULL, NULL };

		/* do we need to save the results? */
		if(m->name != NULL && !re_match(m->name, name, &nres)) {
			return(i);
		}
		if(m->comment != NULL && !re_match(m->comment, comment, &cres)) {
			return(i);
		}
	}

	/* not found; return the last one */
	return(kv_cnt(mii->score.profile) - 1);
}

static _force_inline
toml_table_t *tm_idx_dump_toml(char const *fn)
{
	FILE *fp = fopen(fn, "r");
	if(fp == NULL) {
		error("failed to open file: `%s'. check file path and permission.", fn);
		return(NULL);
	}

	/* parse table */
	size_t const elen = 256;
	char error[elen];
	toml_table_t *table = toml_parse_file(fp, error, elen);
	if(table == NULL) {
		error("error occurred on parsing `%s': %s", fn, error);
		return(NULL);
	}

	/* done */
	fclose(fp);
	return(table);
}

static _force_inline
tm_idx_matcher_t tm_idx_parse_matcher(char const *pname, toml_table_t const *table)
{
	_unused(pname);

	tm_idx_matcher_t matcher = {
		.name    = NULL,
		.comment = NULL
	};

	/* construct sequence name and comment matcher */
	char const *sname   = toml_raw_in(table, "name");
	char const *comment = toml_raw_in(table, "comment");

	if(sname   != NULL) { matcher.name = re_compile(sname); }
	if(comment != NULL) { matcher.comment = re_compile(comment); }
	return(matcher);
}


/* K-A stat and related */
typedef struct {
	double escore;
	double lambda;
	double h;
	double identity;
} tm_idx_stat_t;

static
double tm_idx_calc_stat_lambda(double pp, double lambda, double score)
{
	return(pp * exp(lambda * score));
}

static
double tm_idx_calc_stat_escore(double pp, double unused, double score)
{
	_unused(unused);
	return(pp * score);
}

static
double tm_idx_calc_stat_id(double pp, double lambda, double score)
{
	return(score > 0.0 ? pp * exp(lambda * score) : 0.0);
}

static
double tm_idx_calc_stat_h(double id, double lambda, double score)
{
	if(score > 0.0) {
		return(id * lambda * score);
	} else {
		return((1.0 - id) * lambda * score / 3.0);
	}
}

static _force_inline
double tm_idx_calc_stat_core(double p0, double p1, int8_t const *score_matrix, double (*callback)(double, double, double))
{
	static uint8_t const ridx[4] = { tA, tC, tG, tT };

	double sum = 0.0;
	for(size_t q = 0; q < 4; q++) {
		int8_t const *qrow = &score_matrix[q * DZ_REF_MAT_SIZE];

		for(size_t r = 0; r < 4; r++) {
			double const score = (double)qrow[ridx[r]];
			sum += callback(p0, p1, score);
		}
	}
	return(sum);
}

static _force_inline
double tm_idx_estimate_lambda(double pp, int8_t const *score_matrix)
{
	double est = 1.0;
	double bounds[2] = { 0.0, 2.0 };

	while(bounds[1] - bounds[0] > 0.00001) {
		double const sum = tm_idx_calc_stat_core(pp, est, score_matrix, tm_idx_calc_stat_lambda);
		uint64_t const sup = sum > 1.0;
		bounds[sup] = est;
		est = (est + bounds[1 - sup]) / 2.0;
	}
	return(est);
}

static _force_inline
tm_idx_stat_t tm_idx_calc_stat(int8_t const *score_matrix)
{
	double const pp = 0.25 * 0.25;

	double const lambda   = tm_idx_estimate_lambda(pp, score_matrix);
	double const identity = tm_idx_calc_stat_core(pp, lambda, score_matrix, tm_idx_calc_stat_id);

	return((tm_idx_stat_t){
		.escore = tm_idx_calc_stat_core(pp, 0.0, score_matrix, tm_idx_calc_stat_escore),
		.lambda = lambda,
		.h = tm_idx_calc_stat_core(identity, lambda, score_matrix, tm_idx_calc_stat_h) / 4.0,
		.identity = identity
	});
}

static _force_inline
void tm_idx_fill_stat(tm_idx_profile_t *profile, int8_t const *score_matrix)
{
	_unused(profile);

	tm_idx_stat_t const stat = tm_idx_calc_stat(score_matrix);
	profile->stat.rlambda = 1.0 / stat.lambda;

	debug("e(%f), lambda(%f), H(%f), id(%f)", stat.escore, stat.lambda, stat.h, stat.identity);
	return;
}


typedef void (*tm_idx_score_foreach_t)(void *opaque, int8_t *p, size_t qidx, size_t ridx);

static _force_inline
void tm_idx_score_foreach(tm_idx_profile_t *profile, void *opaque, tm_idx_score_foreach_t fp)
{
	/* do not touch latter half */
	for(size_t q = 0; q < TM_QUERY_ALPH_SIZE; q++) {
		for(size_t r = 0; r < TM_REF_ALPH_SIZE; r++) {
			fp(opaque, &profile->extend.score_matrix[q * DZ_REF_MAT_SIZE + r], q, r);
		}
	}
	return;
}

static
void tm_idx_fill_score_callback(int64_t *s, int8_t *p, size_t qidx, size_t ridx)
{
	uint64_t const matching = tm_qrch_is_match(tm_2bit_to_qch(qidx), tm_rch_pack(ridx));
	*p = s[1 - matching];
	return;
}

static _force_inline
void tm_idx_fill_score(tm_idx_profile_t *profile, int64_t m, int64_t x)
{
	int64_t const s[2] = { m, x };
	tm_idx_score_foreach(profile,
		(void *)s,
		(tm_idx_score_foreach_t)tm_idx_fill_score_callback
	);
	return;
}

// static _force_inline
tm_idx_score_t tm_idx_get_score(tm_idx_profile_t const *profile)
{
	/* we suppose the score matrix be simple match-mismatch model */
	return((tm_idx_score_t){
		.match    = (uint32_t) profile->extend.score_matrix[(nA>>2) * DZ_REF_MAT_SIZE + tA],
		.mismatch = (uint32_t)-profile->extend.score_matrix[(nA>>2) * DZ_REF_MAT_SIZE + tG]
	});
}

static _force_inline
void tm_idx_fill_bonus(tm_idx_profile_t *profile, uint32_t bonus)
{
	debug("called");

	v16i8_t const bv = _set_v16i8(bonus);

	/* copy former half to latter */
	int8_t const *src = profile->extend.score_matrix;
	int8_t *dst = &profile->extend.score_matrix[TM_QUERY_ALPH_SIZE * TM_REF_ALPH_SIZE];

	for(size_t q = 0; q < TM_QUERY_ALPH_SIZE; q++) {
		v16i8_t const x = _loadu_v16i8(&src[q * DZ_REF_MAT_SIZE]);
		v16i8_t const y = _adds_v16i8(x, bv);

		_storeu_v16i8(&dst[q * DZ_REF_MAT_SIZE], y);

		_print_v16i8(x);
		_print_v16i8(y);
	}
	return;
}

static _force_inline
void tm_idx_set_kadj(tm_idx_profile_t *profile, uint32_t k)
{
	profile->chain.kadj[0] = k + 1;
	profile->chain.kadj[1] = k + 1;
	return;
}

static _force_inline
void tm_idx_fill_default(tm_idx_profile_t *profile)
{
	/* default params */
	size_t const k = 4;
	profile->kbits = 2 * k;

	/* seeding and chaining */
	tm_idx_set_kadj(profile, k);
	profile->chain.window.sep.u = 32;
	profile->chain.window.sep.v = 32;
	profile->chain.min_scnt = 4;

	/* filtering */
	profile->filter.span_thresh = UINT32_MAX;	/* will be overridden in tm_idx_calc_filter_thresh */

	/* extension */
	profile->extend.bonus = 10;
	profile->extend.min_score = 30;
	profile->extend.use_raw = 0;
	profile->extend.giv = 5;
	profile->extend.gev = 1;
	profile->extend.gih = 5;
	profile->extend.geh = 1;
	profile->extend.glim = 16;
	tm_idx_fill_score(profile, 2, -3);
	return;
}

static _force_inline
void tm_idx_override_default(tm_idx_profile_t *profile, tm_idx_conf_t const *conf)
{
	/* load values specified by args */

	if(!tm_idx_is_default(conf->kmer)) {
		size_t const kmer = tm_idx_unwrap(conf->kmer);
		profile->kbits = 2 * kmer;
		tm_idx_set_kadj(profile, kmer);
	}

	/* chaining */
	if(!tm_idx_is_default(conf->window)) {
		profile->chain.window.sep.u = tm_idx_unwrap(conf->window);
		profile->chain.window.sep.v = tm_idx_unwrap(conf->window);
	}
	if(!tm_idx_is_default(conf->min_scnt)) {
		uint32_t const min_scnt = tm_idx_unwrap(conf->min_scnt) - 1;
		profile->chain.min_scnt = min_scnt;
	}

	/* filter */
	if(!tm_idx_is_default(conf->span_thresh)) {
		profile->filter.span_thresh = tm_idx_unwrap(conf->span_thresh);
	}

	/* score matrix */
	if(!tm_idx_is_default(conf->match) || !tm_idx_is_default(conf->mismatch)) {
		tm_idx_score_t const s = tm_idx_get_score(profile);
		int64_t const m = tm_idx_is_default(conf->match)    ? s.match    : tm_idx_unwrap(conf->match);
		int64_t const x = tm_idx_is_default(conf->mismatch) ? s.mismatch : tm_idx_unwrap(conf->mismatch);
		tm_idx_fill_score(profile, m, -x);		/* negate */
	}

	/* gap penalties */
	if(!tm_idx_is_default(conf->gap_open)) {
		profile->extend.giv = tm_idx_unwrap(conf->gap_open);
		profile->extend.gih = tm_idx_unwrap(conf->gap_open);
	}
	if(!tm_idx_is_default(conf->gap_extend)) {
		profile->extend.gev = tm_idx_unwrap(conf->gap_extend);
		profile->extend.geh = tm_idx_unwrap(conf->gap_extend);
	}
	if(!tm_idx_is_default(conf->max_gap_len)) {
		profile->extend.glim = tm_idx_unwrap(conf->max_gap_len);
	}

	/* anchoring bonus; score matrix refilled afterward */
	if(!tm_idx_is_default(conf->full_length_bonus)) {
		profile->extend.bonus = tm_idx_unwrap(conf->full_length_bonus);
	}

	/* postprocess */
	if(!tm_idx_is_default(conf->min_score)) {
		profile->extend.min_score = tm_idx_unwrap(conf->min_score);
	}
	if(!tm_idx_is_default(conf->use_raw_score)) {
		profile->extend.use_raw = tm_idx_unwrap(conf->use_raw_score);
	}
	return;
}

typedef struct {
	int64_t acc, min, max;
	size_t cnt;
} tm_idx_calc_acc_t;

static
void tm_idx_acc_filter_score(tm_idx_calc_acc_t *acc, int8_t *p, size_t qidx, size_t ridx)
{
	uint64_t const matching = tm_qrch_is_match(tm_2bit_to_qch(qidx), tm_rch_pack(ridx));
	tm_idx_calc_acc_t *ptr = &acc[1 - matching];

	int64_t const s = *p;
	ptr->acc += s;
	ptr->min = MIN2(ptr->min, s);
	ptr->max = MAX2(ptr->max, s);
	ptr->cnt++;
	return;
}

static _force_inline
void tm_idx_calc_filter_score(tm_idx_profile_t *profile)
{
	/* count and accumulate */
	tm_idx_calc_acc_t c[2] = {
		{ 0, INT32_MAX, INT32_MIN, 0 },
		{ 0, INT32_MAX, INT32_MIN, 0 }
	};
	tm_idx_score_foreach(profile,
		(void *)c,
		(tm_idx_score_foreach_t)tm_idx_acc_filter_score
	);

	/* take average */
	int64_t const mave = (int64_t)c[0].acc / (int64_t)c[0].cnt;
	int64_t const m = (c[0].min == mave
		? c[0].min
		: MAX2(1, mave * 3 / 4)
	);
	debug("(%ld, %ld, %ld, %zu)", c[0].acc, c[0].max, c[0].min, c[0].cnt);

	int64_t const xave = (int64_t)c[1].acc / (int64_t)c[1].cnt;
	int64_t const x = (c[1].max == xave
		? c[1].max
		: MIN2(-1, xave * 3 / 4)
	);
	debug("(%ld, %ld, %ld, %zu)", c[1].acc, c[1].max, c[1].min, c[1].cnt);

	debug("mave(%ld), m(%ld), xave(%ld), x(%ld)", mave, m, xave, x);

	for(size_t i = 0; i < 16; i++) {
		profile->filter.score_matrix[i] = (i == 0) ? x : m;		/* signed */
	}
	return;
}

static _force_inline
void tm_idx_calc_filter_gap(tm_idx_profile_t *profile)
{
	int32_t const ge = (profile->extend.gev + profile->extend.geh) / 2;
	int32_t const gi = (profile->extend.giv + profile->extend.gih) / 8;
	v16i8_t const gv = _set_v16i8(0 - gi - ge);			/* signed; negated without offset */
	_storeu_v16i8(profile->filter.gap, gv);

	debug("ge(%d)", (int8_t)(0 - gi - ge));
	// _print_v16i8(gv);
	return;
}

static _force_inline
void tm_idx_calc_filter_ivec(tm_idx_profile_t *profile)
{
	int32_t const s = 2 * profile->filter.gap[1] - profile->filter.score_matrix[1];
	for(size_t i = 0; i < 8; i++) {
		profile->filter.init[8 + i] = 128 + s *  i;
		profile->filter.init[7 - i] = 128 + s * (i + 1);
	}
	return;
}

static _force_inline
void tm_idx_calc_filter_thresh(tm_idx_profile_t *profile)
{
	tm_idx_stat_t const stat = tm_idx_calc_stat(profile->extend.score_matrix);
	double const min_identity = stat.identity * 0.85;
	int32_t const x = (int8_t)profile->filter.score_matrix[0];
	int32_t const m = (int8_t)profile->filter.score_matrix[1];

	double const min_per_pair = (double)m * min_identity + (double)x * (1.0 - min_identity);
	int32_t const min_score = (int32_t)(16.0 * min_per_pair);

	profile->filter.min_score = MAX2(0, min_score);
	// fprintf(stderr, "id(%f), m(%d), x(%d), min_per_pair(%f), min_score(%d)\n", min_identity, m, x, min_per_pair, min_score);

	/* overwrite span_thresh if the value is the default one */
	if(profile->filter.span_thresh == UINT32_MAX) {
		uint32_t const span_thresh = profile->chain.window.sep.u * 2;
		profile->filter.span_thresh = span_thresh;
	}
	return;
}

static _force_inline
void tm_idx_calc_filter_params(tm_idx_profile_t *profile)
{
	tm_idx_calc_filter_score(profile);
	tm_idx_calc_filter_gap(profile);
	tm_idx_calc_filter_ivec(profile);
	tm_idx_calc_filter_thresh(profile);
	return;
}

static _force_inline
uint64_t tm_idx_check_profile(tm_idx_profile_t const *profile)
{
	#define tm_idx_assert(_cond, ...) ({ if(!(_cond)) { error("" __VA_ARGS__); ecnt++; } })

	size_t ecnt = 0;

	/* kmer */ {
		uint32_t const k = profile->kbits>>1;
		tm_idx_assert(k >= TM_MIN_KMER && k <= TM_MAX_KMER, "k-mer size (-k) must be >= %u and <= %u.", (uint32_t)TM_MIN_KMER, (uint32_t)TM_MAX_KMER);
	}

	/* window size */ {
		uint32_t const w = MAX2(profile->chain.window.sep.u, profile->chain.window.sep.v);
		tm_idx_assert(w <= 256, "chain window size (-w) must be >= 0 and <= 256.");		/* zero to disable chain */
	}

	/* seed count */ {
		uint32_t const c = profile->chain.min_scnt + 1;
		tm_idx_assert(c >= 1 && c <= 1024, "minimum seed count (-c) must be >= 1 and <= 1024.");
	}

	/* score matrix */ {
		tm_idx_calc_acc_t c[2] = {
			{ 0, INT32_MAX, INT32_MIN, 0 },
			{ 0, INT32_MAX, INT32_MIN, 0 }
		};
		tm_idx_score_foreach((tm_idx_profile_t *)profile,
			(void *)c,
			(tm_idx_score_foreach_t)tm_idx_acc_filter_score
		);
		// fprintf(stderr, "%ld, %ld, %ld, %ld\n", c[0].min, c[0].max, c[1].min, c[1].max);

		tm_idx_assert(c[0].min >= 1   && c[0].max <= 31, "match award (-a) must be >= 1 and <= 31.");
		tm_idx_assert(c[1].min >= -31 && c[1].min <= -1, "mismatch penalty (-b) must be >= 1 and <= 31.");
	}
	/* gaps */ {
		uint8_t const giv = profile->extend.giv, gev = profile->extend.gev;
		uint8_t const gih = profile->extend.gih, geh = profile->extend.geh;
		// fprintf(stderr, "%u, %u, %u, %u\n", giv, gev, gih, geh);

		tm_idx_assert(/* giv >= 0 && */ giv <= 31, "gap open penalty (-p) must be >= 0 and <= 31.");
		tm_idx_assert(/* gih >= 0 && */ gih <= 31, "gap open penalty (-p) must be >= 0 and <= 31.");
		tm_idx_assert(gev >= 1 && gev <= 31, "gap extension penalty (-q) must be >= 1 and <= 31.");
		tm_idx_assert(geh >= 1 && geh <= 31, "gap extension penalty (-q) must be >= 1 and <= 31.");
	}

	/* max gap len */ {
		uint32_t const glim = profile->extend.glim;
		// fprintf(stderr, "%u\n", glim);

		tm_idx_assert(glim >= 1 && glim <= 256, "max gap length (-g) must be >= 1 and <= 256.");
	}

	/* full-length bonus */ {
		tm_idx_assert(profile->extend.bonus < 128, "full-length bonus must be smaller than 127.");
	}

	/* min_score */ {
		uint32_t const min_score = profile->extend.min_score;
		tm_idx_assert(min_score <= INT32_MAX, "minimum score threshold must be smaller than INT32_MAX.");
	}
	return(ecnt);
}

static _force_inline
uint64_t tm_idx_finalize_profile(tm_idx_profile_t *profile)
{
	if(tm_idx_check_profile(profile)) {
		return(1);
	}

	/* calc squash interval */
	uint32_t const u    = profile->chain.window.sep.u;
	uint32_t const intv = 0x40000000>>_lzc_u32(u);		/* divide by 2 */
	profile->chain.squash_intv = MAX2(1, intv);			/* load 1 if u == 1 */
	debug("wu(%u), intv(%u)", u, profile->chain.squash_intv);

	/* instanciate dz */
	dz_allocator_t alloc = {
		.ctx = NULL,
		.fp = tm_idx_malloc
	};
	dz_score_conf_t const conf = {
		.score_matrix = profile->extend.score_matrix,
		.ins_open   = profile->extend.giv,
		.ins_extend = profile->extend.gev,
		.del_open   = profile->extend.gih,
		.del_extend = profile->extend.geh,

		/* fixed */
		// .max_ins_len = profile->extend.vlim,
		// .max_del_len = profile->extend.hlim,
		.max_gap_len = profile->extend.glim,
		.full_length_bonus = profile->extend.bonus
	};
	// fprintf(stderr, "(%u, %u, %u, %u)\n", profile->extend.giv, profile->extend.gev, profile->extend.gih, profile->extend.geh);

	profile->extend.dz = dz_init_profile(&alloc, &conf);
	return(0);
}


static _force_inline
size_t tm_idx_profile_size(size_t plen)
{
	size_t const size = sizeof(tm_idx_profile_t) + _roundup(plen + 1, 16);
	return(size);
}

static _force_inline
char *tm_idx_profile_name_ptr(tm_idx_profile_t *profile)
{
	size_t const ofs = sizeof(tm_idx_profile_t);
	return(_add_offset(profile, ofs));
}

static _force_inline
tm_idx_profile_t *tm_idx_allocate_profile(tm_idx_profile_t const *template, char const *pname)
{
	/* base size + profile name length */
	size_t const plen = mm_strlen(pname);
	size_t const size = tm_idx_profile_size(plen);

	/* malloc and fill zero */
	tm_idx_profile_t *profile = malloc(size);	/* calloc? */

	/* copy base image */
	if(template == NULL) {
		memset(profile, 0, size);
	} else {
		memcpy(profile, template, sizeof(tm_idx_profile_t));
	}

	/* name and tail */
	char *p = tm_idx_profile_name_ptr(profile);
	char *t = _add_offset(profile, size);

	/* override save size and name */
	profile->size = size;
	profile->meta.name = p;
	memcpy(p, pname, plen);

	/* clear remaining to avoid use-of-uninitialized-value error; including tail cap ('\0') */
	memset(&p[plen], 0, t - &p[plen]);
	return(profile);
}

// static _force_inline
tm_idx_profile_t *tm_idx_build_profile(tm_idx_conf_t const *conf)
{
	/* name */
	char const *pname = "default";

	/* allocate profile object */
	tm_idx_profile_t *profile = tm_idx_allocate_profile(NULL, pname);

	/* load default params */
	tm_idx_fill_default(profile);
	tm_idx_override_default(profile, conf);

	/* always refill bonus */
	tm_idx_fill_bonus(profile, profile->extend.bonus);
	tm_idx_fill_stat(profile, profile->extend.score_matrix);

	/* derive filtering parameters from extension parameters */
	tm_idx_calc_filter_params(profile);

	/* instanciate dz */
	if(tm_idx_finalize_profile(profile)) {
		free(profile);
		return(NULL);
	}
	return(profile);
}



/* argument parsers */

typedef uint64_t (*tm_idx_set_t)(tm_idx_profile_t *profile, void *str);
typedef void *(*tm_idx_toml_in_t)(toml_table_t const *tab, char const *key);
typedef struct {
	char const *key;
	size_t prefix;
	tm_idx_set_t set;
	tm_idx_toml_in_t in;
} tm_idx_parse_t;

typedef struct {
	uint8_t idx[TM_REF_ALPH_SIZE];
	uint64_t has_header;
	size_t qsize, rsize;
} tm_idx_dim_t;


static _force_inline
uint64_t tm_idx_is_valid_base(uint8_t const *conv, char const *str)
{
	return(str[0] == '\0' || str[1] != '\0' || conv[(size_t)(uint8_t)str[0]] == 0);
}

static _force_inline
tm_idx_dim_t tm_idx_parse_header(toml_array_t const *arr, tm_idx_dim_t dim)
{
	/* index conversion table */
	static uint8_t const conv[256] = {
		['A'] = A, ['C'] = C, ['G'] = G, ['T'] = T,
		['R'] = R, ['Y'] = Y, ['S'] = S,
		['W'] = W, ['K'] = K, ['M'] = M,
		['B'] = B, ['D'] = D, ['H'] = H, ['V'] = V,
		['N'] = 16		/* to distinguish 'N' with invalid characters */
	};

	tm_idx_dim_t const error = { { 0 }, 0 };		/* qsize == 0 for error */
	toml_array_t const *hdr = toml_array_at(arr, 0);

	/* update dim */
	dim.has_header = 1;
	dim.qsize--;
	dim.rsize = toml_array_nelem(hdr);		/* any number is acceptable */

	/* save index for each base in the header */
	for(size_t ridx = 0; ridx < dim.rsize; ridx++) {
		char const *bstr = toml_raw_at(hdr, ridx);
		if(tm_idx_is_valid_base(conv, bstr)) {
			error("invalid base character in score matrix header: `%s'.", bstr);
			return(error);
		}
		dim.idx[ridx] = conv[(size_t)(uint8_t)bstr[0]];
	}

	/* the header is confirmed valid */
	return(dim);
}

static _force_inline
size_t tm_idx_determine_rsize(toml_array_t const *arr, tm_idx_dim_t dim)
{
	size_t rsize[TM_QUERY_ALPH_SIZE + 1] = { 0 };

	/* copy lengths */
	for(size_t i = 0; i < dim.qsize + dim.has_header; i++) {
		rsize[i] = toml_array_nelem(toml_array_at(arr, i));
	}

	/* check column lengths are matched */
	for(size_t i = 1; i < dim.qsize + dim.has_header; i++) {
		if(rsize[i - 1] != rsize[i]) {
			error("unmatched row lengths of score matrix.");
			return(0);
		}
	}

	/* matched; anything else? */
	return(rsize[0]);
}

static _force_inline
tm_idx_dim_t tm_idx_fill_idx(toml_array_t const *arr, tm_idx_dim_t dim)
{
	tm_idx_dim_t const error = { { 0 }, 0 };		/* qsize == 0 for error */

	/* here header is missing; we expect matrix is 4 x 4 or 16 x 4 */
	dim.rsize = tm_idx_determine_rsize(arr, dim);
	if(dim.rsize != 4 || dim.rsize != 16) { return(error); }

	/* prebuilt conversion table */
	static uint8_t const idx[2][TM_REF_ALPH_SIZE] __attribute__(( aligned(16) )) = {
		{ 1, 2, 4, 8 },		/* A, C, G, T */
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }	/* FIXME? id(x) */
	};
	memcpy(&dim.idx, &idx[dim.rsize == 16], TM_REF_ALPH_SIZE);
	return(dim);
}

static _force_inline
tm_idx_dim_t tm_idx_parse_dim(toml_array_t const *arr)
{
	tm_idx_dim_t const error = { { 0 }, 0 };		/* qsize == 0 for error */
	tm_idx_dim_t dim = {
		.has_header = 0,
		.qsize = toml_array_nelem(arr),
		.rsize = 0
	};

	if(dim.qsize == TM_QUERY_ALPH_SIZE + 1) {
		/* table seems to have header; inspect the header row */
		return(tm_idx_parse_header(arr, dim));

	} else if(dim.qsize == TM_QUERY_ALPH_SIZE) {
		/* header missing; we expect the matrix is 4 x 4 or 16 x 4 */
		return(tm_idx_fill_idx(arr, dim));
	}

	/* no header available and #rows does not match */
	error("wrong #rows for score_matrix. four rows (for each of query-side A, C, G, and T) expected.");
	return(error);
}

static _force_inline
uint64_t tm_idx_parse_row(tm_idx_profile_t *profile, toml_array_t const *row, size_t qidx, tm_idx_dim_t dim)
{
	static uint8_t const rconv[16] __attribute__(( aligned(16) )) = {
		[tA] = A, [tC] = C, [tG] = G, [tT] = T,
		[tR] = R, [tY] = Y, [tS] = S, [tS ^ 0x0f] = S,
		[tK] = K, [tM] = M, [tW] = W, [tW ^ 0x0f] = W,
		[tB] = B, [tD] = D, [tH] = H, [tV] = V
	};
	v16i8_t const mv = _load_v16i8(rconv);

	/* load current row */
	int8_t *p = &profile->extend.score_matrix[qidx * DZ_REF_MAT_SIZE];
	v16i8_t v = _loadu_v16i8(p);

	for(size_t ridx = 0; ridx < dim.rsize; ridx++) {
		/* retrieve string at (qidx, ridx) */
		char const *val = toml_raw_at(row, ridx);

		/* convert to integer, INT64_MIN is returned when unparsable */
		int64_t const n = mm_atoi(val, 0);
		if(n == INT64_MIN) { return(1); }

		/* save (others are left -DZ_SCORE_OFS) */
		uint8_t const ch = dim.idx[ridx];	/* 4bit encoded */
		v16i8_t const cv = _set_v16i8(ch);
		v16i8_t const eq = _eq_v16i8(mv, cv);

		/* fold char into vector */
		v16i8_t const nv = _set_v16i8(n);
		v = _sel_v16i8(eq, nv, v);
	}

	/* write back */
	_storeu_v16i8(p, v);

	/* we need to patch bonus afterward */
	return(0);
}

static
uint64_t tm_idx_parse_matrix(tm_idx_profile_t *profile, toml_array_t const *arr)
{
	/* must be two-dimensional array of integer */
	if(toml_array_kind(arr) == 'a') {
		error("two-dimensional integer matrix expected for score_matrix in parsing `%s'.", profile->meta.name);
		return(1);
	}

	/* calc dimension */
	tm_idx_dim_t const dim = tm_idx_parse_dim(arr);
	if(dim.qsize == 0) {
		goto _tm_idx_parse_matrix_fail;		/* error occurred */
	}

	/* clear matrix */
	memset(profile->extend.score_matrix, -DZ_SCORE_OFS, DZ_REF_MAT_SIZE * DZ_QUERY_MAT_SIZE);

	/* convert to integer array */
	for(size_t qidx = 0; qidx < dim.qsize; qidx++) {
		toml_array_t const *row = toml_array_at(arr, qidx + dim.has_header);

		if(tm_idx_parse_row(profile, row, qidx, dim)) {
			goto _tm_idx_parse_matrix_fail;
		}
	}

	/* done */
	return(0);

_tm_idx_parse_matrix_fail:;
	error("in parsing `%s'.", profile->meta.name);
	return(1);
}

static uint64_t tm_idx_set_kmer(tm_idx_profile_t *profile, char const *str)
{
	int64_t const k = mm_atoi(str, 0);
	profile->kbits = 2 * k;
	tm_idx_set_kadj(profile, k);
	return(0);
}
static uint64_t tm_idx_set_window(tm_idx_profile_t *profile, char const *str) {
	int64_t const w = mm_atoi(str, 0);
	profile->chain.window.sep.u = w;
	profile->chain.window.sep.v = w;
	return(0);
}
static uint64_t tm_idx_set_ccnt(tm_idx_profile_t *profile, char const *str) {
	profile->chain.min_scnt = mm_atoi(str, 0);
	return(0);
}
static uint64_t tm_idx_set_filter(tm_idx_profile_t *profile, char const *str) {
	if(mm_strcmp(str, "false") == 0) {
		profile->filter.span_thresh = UINT32_MAX;
	} else if(mm_strcmp(str, "true") != 0) {
		error("boolean (true or false) expected for enable-filter (%s) in parsing `%s'.", str, profile->meta.name);
		return(1);
	}
	return(0);
}
static uint64_t tm_idx_set_gap_max(tm_idx_profile_t *profile, char const *str) {
	profile->extend.glim = mm_atoi(str, 0);
	return(0);
}
static uint64_t tm_idx_set_ins_open(tm_idx_profile_t *profile, char const *str) {
	profile->extend.giv = mm_atoi(str, 0);
	return(0);
}
static uint64_t tm_idx_set_ins_extend(tm_idx_profile_t *profile, char const *str) {
	profile->extend.gev = mm_atoi(str, 0);
	return(0);
}
static uint64_t tm_idx_set_del_open(tm_idx_profile_t *profile, char const *str) {
	profile->extend.gih = mm_atoi(str, 0);
	return(0);
}
static uint64_t tm_idx_set_del_extend(tm_idx_profile_t *profile, char const *str) {
	profile->extend.geh = mm_atoi(str, 0);
	return(0);
}
static uint64_t tm_idx_set_min_score(tm_idx_profile_t *profile, char const *str) {
	profile->extend.min_score = mm_atoi(str, 0);
	return(0);
}



static _force_inline
uint64_t tm_idx_match_key(char const *key, char const *expected, size_t min_len)
{
	size_t i = 0;
	for(; key[i] != '\0' && expected[i] != '\0'; i++) {
		if(key[i] != expected[i]) { return(0); }
	}
	return(i >= min_len && key[i] == '\0');		/* regard as match if key is shorter than full-length key */
}

static _force_inline
uint64_t tm_idx_parse_key(tm_idx_profile_t *profile, toml_table_t const *table, char const *key)
{
	/* use the same parser for command-line arguments for consistency */
	static tm_idx_parse_t const keys[] = {
		/* short option */
		{ "kmer_size",     4, (tm_idx_set_t)tm_idx_set_kmer,       (tm_idx_toml_in_t)toml_raw_in },
		{ "window_size",   6, (tm_idx_set_t)tm_idx_set_window,     (tm_idx_toml_in_t)toml_raw_in },
		{ "min_seed_cnt",  8, (tm_idx_set_t)tm_idx_set_ccnt,       (tm_idx_toml_in_t)toml_raw_in },
		{ "enable_filter", 6, (tm_idx_set_t)tm_idx_set_filter,     (tm_idx_toml_in_t)toml_raw_in },		/* disable filter if arg == "false" */
		{ "gap_max_len",   7, (tm_idx_set_t)tm_idx_set_gap_max,    (tm_idx_toml_in_t)toml_raw_in },
		{ "ins_open",      6, (tm_idx_set_t)tm_idx_set_ins_open,   (tm_idx_toml_in_t)toml_raw_in },
		{ "ins_extend",    6, (tm_idx_set_t)tm_idx_set_ins_extend, (tm_idx_toml_in_t)toml_raw_in },
		{ "del_open",      6, (tm_idx_set_t)tm_idx_set_del_open,   (tm_idx_toml_in_t)toml_raw_in },
		{ "del_extend",    6, (tm_idx_set_t)tm_idx_set_del_extend, (tm_idx_toml_in_t)toml_raw_in },
		{ "min_score",     6, (tm_idx_set_t)tm_idx_set_min_score,  (tm_idx_toml_in_t)toml_raw_in },
		{ "score_matrix",  5, (tm_idx_set_t)tm_idx_parse_matrix,   (tm_idx_toml_in_t)toml_table_in },

		{ NULL, 0, NULL, NULL }		/* sentinel */
	};

	/* parse int / boolean */
	for(tm_idx_parse_t const *p = keys; p->key != NULL; p++) {
		if(!tm_idx_match_key(key, p->key, p->prefix)) { continue; }

		return(p->set(profile, (void *)p->in(table, key)));
	}

	error("unrecognized key `%s' in parsing `%s'.", key, profile->meta.name);
	return(1);		/* not found */
}

static _force_inline
tm_idx_profile_t *tm_idx_parse_table(tm_idx_profile_t const *template, char const *pname, toml_table_t const *table)
{
	/* copy params */
	tm_idx_profile_t *profile = tm_idx_allocate_profile(template, pname);

	/* parse table */
	char const *key = NULL;
	for(size_t i = 0; (key = toml_key_in(table, i)) != NULL; i++) {
		if(tm_idx_parse_key(profile, table, key)) {
			error("error in parsing table: `%s'.", profile->meta.name);
			return(NULL);
		}
	}

	/* always recalc filtering scores and refill bonus */
	tm_idx_fill_bonus(profile, profile->extend.bonus);
	tm_idx_fill_stat(profile, profile->extend.score_matrix);
	tm_idx_calc_filter_params(profile);

	/* instanciate dz */
	if(tm_idx_finalize_profile(profile)) {
		free(profile);
		return(NULL);
	}
	return(profile);
}

static _force_inline
uint64_t tm_idx_load_profile_core(tm_idx_gen_t *mii, tm_idx_profile_t const *template, toml_table_t const *root)
{
	char const *key;
	for(size_t i = 0; (key = toml_key_in(root, i)) != NULL; i++) {
		toml_table_t const *table = toml_table_in(root, key);

		/* check if the key has a corresponding table */
		if(table == NULL) {
			error("missing tabular content for key `%s'. ignoring.", key);
			continue;
		}

		/* parse name / comment matcher */
		tm_idx_matcher_t matcher = tm_idx_parse_matcher(key, table);
		if(matcher.name == NULL && matcher.comment == NULL) { continue; }	/* error when both NULL */

		/* parse score matrix */
		tm_idx_profile_t *profile = tm_idx_parse_table(template, key, table);	/* NULL if error */
		if(profile == NULL) { return(1); }	/* unignorable error */

		/* save */
		tm_idx_push_profile(mii, matcher, profile);
	}
	return(0);			/* no error */
}

static _force_inline
uint64_t tm_idx_load_profile(tm_idx_gen_t *mii, tm_idx_profile_t const *template, char const *fn)
{
	/* open file as toml */
	toml_table_t *root = tm_idx_dump_toml(fn);	/* error message printed inside */
	if(root == NULL) { return(1); }

	/* parse */
	uint64_t const state = tm_idx_load_profile_core(mii, template, root);
	toml_free(root);
	return(state);
}



/* profile builder API and multithreading */

static _force_inline
uint64_t tm_idx_gen_parr(tm_idx_gen_t *mii, tm_idx_conf_t const *conf, char const *fn, FILE *lfp)
{
	/* compose default (fallback) profile as template */
	tm_idx_profile_t *fallback = tm_idx_build_profile(conf);
	if(fallback == NULL) {
		error("failed load default score profile.");
		return(1);
	}

	if(fn != NULL) {
		message(lfp, "reading score matrices...");
		if(tm_idx_load_profile(mii, fallback, fn)) {
			error("failed construct score profile.");
			return(1);
		}
	} else {
		message(lfp, "loading default score matrix...");
	}

	/* append default score matrix as fallback */
	tm_idx_matcher_t m = { NULL, NULL };
	tm_idx_push_profile(mii, m, fallback);
	return(0);
}


static
tm_idx_batch_t *tm_idx_fetch(uint32_t tid, tm_idx_gen_t *self)
{
	_unused(tid);		/* always run on the main thread */

	if(bseq_is_eof(self->col.fp)) { return(NULL); }	/* reached tail */

	bseq_batch_t *bin = bseq_read(self->col.fp);
	if(bin == NULL) { return(NULL); }		/* reached tail (or something is wrong?) */

	/* assign block id */
	tm_idx_batch_t *batch = _sub_offset(bin, offsetof(tm_idx_batch_t, bin));

	/* update sequence count */
	size_t const scnt = kv_cnt(self->col.bin);
	batch->base_rid = scnt;

	/* expand array */
	size_t const nscnt = scnt + bseq_meta_cnt(&batch->bin);
	kv_reserve(tm_idx_sketch_t *, self->col.bin, nscnt);
	kv_cnt(self->col.bin) = nscnt;

	return(batch);
}

static _force_inline
size_t tm_idx_roundup(size_t len)
{
	/* keep separation larger than 12 bases */
	return(_roundup(len + 12, 16));
}

static _force_inline
tm_ref_sketch_t *tm_idx_build_sketch(tm_idx_gen_t *self, tm_ref_tbuf_t *ref, size_t pid, char const *name, size_t nlen, uint8_t const *seq, size_t slen)
{
	_unused(name);

	/* calc margin size for saving name and seq */
	size_t const margin = tm_idx_roundup(nlen) + tm_idx_roundup(slen);

	/* retrieve profile */
	tm_idx_profile_t const *profile = kv_ptr(self->score.profile)[pid];

	/* if fail (due to too many k-mers), try again with reduced ambiguity */
	for(size_t i = 0; i < TM_MAX_SKETCH_TRIAL; i++) {
		/* pack args */
		tm_ref_conf_t const conf = {
			.kbits       = profile->kbits,
			.max_amb_col = TM_MAX_AMB_COLUMNS>>i,
			.margin = {
				.head = sizeof(tm_idx_sketch_hdr_t),
				.tail = margin + 32			/* add 32 for vectorized access */
			}
		};
		tm_ref_sketch_t *sk = tm_ref_sketch(ref, &conf, seq, slen);
		if(sk != NULL) { return(sk); }

		warn("failed to pack link for `%.*s' with max ambiguity = %u. %s.",
			(int)nlen, name, (uint32_t)(TM_MAX_AMB_COLUMNS>>i),
			i == (TM_MAX_SKETCH_TRIAL - 1) ? "give up" : "try again"
		);
	}

	/* failure */
	return(NULL);
}

static _force_inline
tm_idx_sketch_t *tm_idx_save_seq(tm_ref_sketch_t *sr, uint32_t pid, char const *name, size_t nlen, uint8_t const *seq, size_t slen)
{
	v32i8_t const z = _zero_v32i8();

	/* copy name */
	uint8_t *n = _add_offset(sr, sr->size);
	for(size_t i = 0; i < nlen; i++) {
		v16i8_t const v = _loadu_v16i8(&name[i]);
		_storeu_v16i8(&n[i], v);
	}
	_storeu_v32i8(&n[nlen], z);

	/* copy sequence */
	uint8_t *s = _add_offset(n, tm_idx_roundup(nlen));
	for(size_t i = 0; i < slen; i += 16) {
		v16i8_t const v = _loadu_v16i8(&seq[i]);
		_storeu_v16i8(&s[i], v);
	}
	_storeu_v32i8(&s[slen], z);

	/* build header */
	tm_idx_sketch_t *si = _sub_offset(sr, sizeof(tm_idx_sketch_hdr_t));
	size_t const size = (
		  sizeof(tm_idx_sketch_hdr_t)	/* header */
		+ sr->size						/* body */
		+ tm_idx_roundup(nlen)			/* name (margined) */
		+ tm_idx_roundup(slen)			/* sequence (margined) */
	);
	si->h = (tm_idx_sketch_hdr_t){
		.size = size,
		.rid  = 0,
		.pid  = pid,
		.sofs = size - tm_idx_roundup(slen)
	};
	return(si);
}

static
tm_idx_batch_t *tm_idx_collect(uint32_t tid, tm_idx_gen_t *self, tm_idx_batch_t *batch)
{
	/* load pointers */
	tm_ref_tbuf_t *ref = self->ref[tid];
	bseq_meta_t *meta  = bseq_meta_ptr(&batch->bin);

	for(size_t i = 0; i < bseq_meta_cnt(&batch->bin); i++) {
		bseq_meta_t *p = &meta[i];

		p->u.ptr = NULL;		/* clear first */
		if(bseq_seq_len(p) > INT16_MAX) {
			debugblock({ fprintf(stderr, "i(%zu), len(%zu), seq(%s)\n", batch->bin.base_id + i, bseq_seq_len(p), bseq_name(p)); });
			continue;
		}

		/*
		debug("name(%s)", bseq_name(p));
		for(size_t j = 0; j < bseq_seq_len(p); j++) {
			fprintf(stderr, "%c", "NACMGRSVTWYHKDBN"[bseq_seq(p)[j]]);
		}
		fprintf(stderr, "\n");
		*/

		/* get score matrix for this sequence from its name with regex matching */
		size_t const pid = tm_idx_find_profile(self,
			bseq_name(p),    bseq_name_len(p),
			bseq_comment(p), bseq_comment_len(p)
		);
		// debug("pid(%zu)", pid);

		/* build hash table */
		tm_ref_sketch_t *sr = tm_idx_build_sketch(self,
			ref, pid,
			bseq_name(p), bseq_name_len(p),
			bseq_seq(p),  bseq_seq_len(p)
		);
		debug("done, sr(%p)", sr);
		if(sr == NULL) {
			debugblock({ fprintf(stderr, "i(%zu), len(%zu), seq(%s)\n", batch->bin.base_id + i, bseq_seq_len(p), bseq_name(p)); });
			continue;
		}

		/* copy sequence */
		tm_idx_sketch_t *si = tm_idx_save_seq(
			sr, pid,	/* save profile id */
			bseq_name(p), bseq_name_len(p),
			bseq_seq(p),  bseq_seq_len(p)			
		);

		/* save ptr */
		p->u.ptr = si;
	}
	return(batch);
}

static
void tm_idx_record(uint32_t tid, tm_idx_gen_t *self, tm_idx_batch_t *batch)
{
	_unused(tid);

	/* we expect tm_idx_fetch and tm_idx_record are not run concurrently */
	tm_idx_sketch_t **arr = kv_ptr(self->col.bin);
	bseq_meta_t const *meta = bseq_meta_ptr(&batch->bin);

	/* no need for sorting batches */
	size_t const base_rid = batch->base_rid;
	for(size_t i = 0; i < bseq_meta_cnt(&batch->bin); i++) {
		bseq_meta_t const *p = &meta[i];
		if(_unlikely(p->u.ptr == NULL && bseq_seq_len(p))) {
			error("sequence too long (%.*s: length = %zu). removed.", (int)bseq_name_len(p), bseq_name(p), bseq_seq_len(p));
		}
		arr[base_rid + i] = p->u.ptr;		/* just copy */
	}

	/* done */
	bseq_free(&batch->bin);
	return;
}

static _force_inline
uint64_t tm_idx_check_sanity(tm_idx_gen_t *mii, tm_idx_conf_t const *conf, char const *fn)
{
	tm_idx_sketch_t const **arr = (tm_idx_sketch_t const **)kv_ptr(mii->col.bin);

	/* make sure at least one sequence is valid */
	size_t const total = kv_cnt(mii->col.bin);
	size_t valid = 0;
	for(size_t i = 0; i < total; i++) {
		valid += arr[i] != NULL;
	}

	/* if strict mode; make sure #valid == #seq */
	if(conf->ignore_overflow == 0 && valid != total) {
		error("pass -L to ignore removed sequences.");
		return(1);
	}

	/* otherwise at least one valid sequence is enough */
	if(valid == 0) {
		error("no valid sequence found in `%s'. sequence(s) might be too long.", fn);
		return(1);
	}
	return(0);
}

static _force_inline
uint64_t tm_idx_gen_core_mt(tm_idx_gen_t *mii, pt_t *pt)
{
	/* init dst array */
	kv_init(mii->col.bin);
	for(size_t i = 0; i < pt_nth(pt); i++) {
		mii->ref[i] = tm_ref_tbuf_init();
	}

	/* stream */
	pt_stream(pt, mii,
		(pt_source_t)tm_idx_fetch,
		(pt_worker_t)tm_idx_collect,
		(pt_drain_t)tm_idx_record		/* no need for sorting */
	);

	/* done */
	for(size_t i = 0; i < pt_nth(pt); i++) {
		tm_ref_tbuf_destroy(mii->ref[i]);
	}
	return(0);
}

static _force_inline
uint64_t tm_idx_gen_core(tm_idx_gen_t *mii, tm_idx_conf_t const *conf, char const *fn, pt_t *pt)
{
	_unused(conf);

	/* do not keep comment nor qual, do not use head margin */
	bseq_conf_t const bseq_conf = {
		.batch_size  = TM_BATCH_SIZE,
		.head_margin = sizeof(tm_idx_batch_t),

		/* irregular 4bit encoding for reference */
		#define _c(x)		( tm_rch_pack(x) )		/* make different from '\0' */
		.conv_table  = {
			[A] = _c(tA), [C] = _c(tC), [G] = _c(tG), [T] = _c(tT),

			[R] = _c(tR), [Y] = _c(tY), [S] = _c(tS),
			[K] = _c(tK), [M] = _c(tM), [W] = _c(tW),

			[B] = _c(tB), [D] = _c(tD), [H] = _c(tH), [V] = _c(tV),

			[N] = _c(tR)		/* FIXME */
		}
		#undef _c
	};
	mii->col.fp = bseq_open(&bseq_conf, fn);
	if(mii->col.fp == NULL) {
		error("failed to open file `%s', please check path and permission.", fn);
		return(1);
	}

	/* multithreading */
	if(tm_idx_gen_core_mt(mii, pt)) {
		return(2);
	}
	bseq_close_t const c = bseq_close(mii->col.fp);

	/* sanity check */
	uint64_t e = tm_idx_check_sanity(mii, conf, fn);
	if(c.status != 0) {
		warn("broken file format detected for `%s'.", fn);
	}
	if(c.cnt != kv_cnt(mii->col.bin)) {
		error("#sequences do not match (may be a bug).")
		e = 1;
	}

	// debug("read(%zu, %zu)", c.cnt, kv_cnt(mii->col.bin));
	return(e);
}

static _force_inline
void tm_idx_destroy_matcher(tm_idx_gen_t *mii)
{
	tm_idx_matcher_t const *matcher = kv_ptr(mii->score.matcher);
	size_t const mcnt = kv_cnt(mii->score.matcher);

	for(size_t i = 0; i < mcnt; i++) {
		tm_idx_matcher_t const *m = &matcher[i];

		if(m->name != NULL) { re_free(m->name); }
		if(m->comment != NULL) { re_free(m->comment); }
	}
	kv_destroy(mii->score.matcher);
	return;
}

size_t tm_idx_squash_invalid(tm_idx_sketch_t **sk, size_t scnt)
{
	size_t cnt = 0;
	for(size_t i = 0; i < scnt; i++) {
		if(sk[i] == NULL) { continue; }

		/* copy and assign rid */
		sk[cnt] = sk[i];
		sk[cnt]->h.rid = cnt;
		cnt++;
	}
	return(cnt);
}


static _force_inline
tm_idx_t *tm_idx_gen_finalize(tm_idx_gen_t *mii, char const *ref, char const *profile)
{
	/* cleanup regex matcher */
	tm_idx_destroy_matcher(mii);

	/* save filename */
	tm_idx_t *mi = &mii->mi;
	mi->filename.ref     = mm_strdup(ref);
	mi->filename.profile = mm_strdup(profile);

	/* copy pointers */
	mi->profile.arr = kv_ptr(mii->score.profile);
	mi->profile.cnt = kv_cnt(mii->score.profile);

	mi->sketch.arr = kv_ptr(mii->col.bin);
	mi->sketch.cnt = tm_idx_squash_invalid(		/* remove NULL elements */
		kv_ptr(mii->col.bin),
		kv_cnt(mii->col.bin)
	);

	/* clear working buffers */
	memset(_add_offset(mi, sizeof(tm_idx_t)), 0, sizeof(tm_idx_gen_t) - sizeof(tm_idx_t));
	return(mi);
}

// static _force_inline
tm_idx_t *tm_idx_gen(tm_idx_conf_t const *conf, pt_t *pt, char const *ref, char const *profile, FILE *lfp)
{
	/* allocate thread-local buffers */
	size_t size = sizeof(tm_idx_gen_t) + pt_nth(pt) * sizeof(tm_ref_tbuf_t *);
	tm_idx_gen_t *mii = aligned_malloc(size);
	memset(mii, 0, size);

	/* load score matrices and build regex matcher */
	if(tm_idx_gen_parr(mii, conf, profile, lfp)) {
		goto _tm_idx_gen_error;
	}

	/* read sequence and collect k-mers */
	message(lfp, "reading sequences and collecting k-mers...");
	if(tm_idx_gen_core(mii, conf, ref, pt)) {
		goto _tm_idx_gen_error;
	}

	message(lfp, "done.");
	return(tm_idx_gen_finalize(mii, ref, profile));

_tm_idx_gen_error:;
	kv_destroy(mii->col.bin);
	tm_idx_destroy(&mii->mi);
	return(NULL);
}


/* index I/O */
#define TM_IDX_MAGIC				( 0x35494d54 )

typedef struct {
	uint64_t magic;
	size_t size;
} tm_idx_hdr_t;

typedef struct {
	void *fp;
	write_t wfp;
	size_t ofs;
} tm_idx_dump_t;

/* dump utilities; offset accumulator */
static _force_inline
void *tm_idx_acc_ofs(size_t *p, size_t size)
{
	size_t ofs = *p;
	*p += size;
	return((void *)ofs);
}

static _force_inline
size_t tm_idx_padding_size(size_t ofs)
{
	return((64ULL - ofs) & (64ULL - 1));
}

static _force_inline
size_t tm_idx_string_size(char const *s)
{
	return(mm_strlen(s) + 1);
}

/* dump block and return offset from dumped object as void * */
static _force_inline
void *tm_idx_dump_block(tm_idx_dump_t *w, void const *p, size_t size)
{
	w->wfp(w->fp, p, size);
	return(tm_idx_acc_ofs(&w->ofs, size));
}

static _force_inline
void *tm_idx_dump_pad(tm_idx_dump_t *w)
{
	size_t size = tm_idx_padding_size(w->ofs);
	if(size == 0) { return((void *)((uintptr_t)w->ofs)); }

	static uint8_t const pad[64] = { 0 };
	return(tm_idx_dump_block(w, pad, size));
}

static _force_inline
void *tm_idx_dump_string(tm_idx_dump_t *w, char const *s)
{
	char const c = '\0';
	char const *p = s == NULL ? &c : s;
	size_t const l = tm_idx_string_size(s);
	return(tm_idx_dump_block(w, p, l));
}


static _force_inline
size_t tm_idx_scan_profile(tm_idx_t const *mi)
{
	size_t size = 0;

	for(size_t i = 0; i < mi->profile.cnt; i++) {
		size += mi->profile.arr[i]->size;		/* no padding between profiles */
	}
	size += tm_idx_padding_size(size);

	/* pointer table */
	size += sizeof(tm_idx_profile_t *) * mi->profile.cnt;
	size += tm_idx_padding_size(size);
	return(size);
}

static _force_inline
void *tm_idx_dump_profile_core(tm_idx_dump_t *w, tm_idx_profile_t *profile)
{
	uint8_t buf[profile->size + 1];
	tm_idx_profile_t *p = (tm_idx_profile_t *)buf;
	memcpy(p, profile, profile->size);

	/* ignore dz_profile_t object */
	p->meta.name = NULL;
	p->extend.dz = NULL;
	return(tm_idx_dump_block(w, p, p->size));
}

static _force_inline
void *tm_idx_dump_profile(tm_idx_dump_t *w, tm_idx_t const *mi)
{
	/* copy pointer array */
	tm_idx_profile_t **buf = (tm_idx_profile_t **)malloc(sizeof(tm_idx_profile_t *) * mi->profile.cnt);
	memcpy(buf, mi->profile.arr, sizeof(tm_idx_profile_t *) * mi->profile.cnt);

	/* dump body */
	debug("pcnt(%zu)", mi->profile.cnt);
	for(size_t i = 0; i < mi->profile.cnt; i++) {
		buf[i] = tm_idx_dump_profile_core(w, buf[i]);
		debug("i(%zu), ofs(%zu)", i, (size_t)buf[i]);
		/* we don't add padding between profiles */
	}
	tm_idx_dump_pad(w);

	/* dump offsets (pointers) */
	void *ofs = tm_idx_dump_block(w, buf, sizeof(tm_idx_profile_t *) * mi->profile.cnt);
	tm_idx_dump_pad(w);

	/* done */
	free(buf);
	return(ofs);
}

static _force_inline
size_t tm_idx_scan_sketch(tm_idx_t const *mi)
{
	size_t size = 0;

	for(size_t i = 0; i < mi->sketch.cnt; i++) {
		size += mi->sketch.arr[i]->h.size;
		size += tm_idx_padding_size(size);
	}

	/* pointer table */
	size += sizeof(tm_idx_sketch_t *) * mi->sketch.cnt;
	size += tm_idx_padding_size(size);
	return(size);
}

static _force_inline
void *tm_idx_dump_sketch(tm_idx_dump_t *w, tm_idx_t const *mi)
{
	/* copy pointer array */
	tm_idx_sketch_t **buf = (tm_idx_sketch_t **)malloc(sizeof(tm_idx_sketch_t *) * mi->sketch.cnt);
	memcpy(buf, mi->sketch.arr, sizeof(tm_idx_sketch_t *) * mi->sketch.cnt);

	/* dump body */
	for(size_t i = 0; i < mi->sketch.cnt; i++) {
		buf[i] = tm_idx_dump_block(w, buf[i], buf[i]->h.size);
		tm_idx_dump_pad(w);		/* make aligned */
	}

	/* dump offsets */
	void *ofs = tm_idx_dump_block(w, buf, sizeof(tm_idx_sketch_t *) * mi->sketch.cnt);
	tm_idx_dump_pad(w);

	/* done */
	free(buf);
	return(ofs);
}

static _force_inline
size_t tm_idx_calc_size(tm_idx_t const *mi)
{
	size_t size = 0;

	size += tm_idx_scan_profile(mi);
	size += tm_idx_scan_sketch(mi);

	size += tm_idx_string_size(mi->filename.ref);
	size += tm_idx_string_size(mi->filename.profile);
	size += tm_idx_padding_size(size);		/* make aligned */

	size += sizeof(tm_idx_t);
	return(size);
}

// static _force_inline
size_t tm_idx_dump(tm_idx_t const *mi, void *fp, write_t wfp)
{
	/* copy on stack */
	tm_idx_t cmi = *mi;

	/* first compute entire size */
	tm_idx_hdr_t hdr = {
		.magic = TM_IDX_MAGIC,
		.size  = tm_idx_calc_size(&cmi)
	};
	wfp(fp, &hdr, sizeof(tm_idx_hdr_t));

	/* init dump context */
	tm_idx_dump_t w = {
		.fp = fp,
		.wfp = wfp,
		.ofs = 0
	};
	cmi.profile.arr = tm_idx_dump_profile(&w, &cmi);
	cmi.sketch.arr  = tm_idx_dump_sketch(&w, &cmi);

	/* input seqence names */
	cmi.filename.ref     = tm_idx_dump_string(&w, cmi.filename.ref);
	cmi.filename.profile = tm_idx_dump_string(&w, cmi.filename.profile);

	/* done */
	tm_idx_dump_pad(&w);
	tm_idx_dump_block(&w, &cmi, sizeof(tm_idx_t));

	debug("size(%zu, %zu)", w.ofs, hdr.size);
	return(w.ofs);
}


/* laod */
static _force_inline
tm_idx_t *tm_idx_slice_root(uint8_t *base, size_t size)
{
	debug("size(%zu)", size);
	tm_idx_t *mi = _sub_offset(&base[size], sizeof(tm_idx_t));

	/* adjust bkt and seq pointers */
	mi->profile.arr = _add_offset(mi->profile.arr, base);
	mi->sketch.arr  = _add_offset(mi->sketch.arr, base);

	/* save base pointer (to be freed correctly) */
	mi->base = (void *)base;
	return(mi);
}

static _force_inline
void tm_idx_patch_ptr(tm_idx_t *mi, uint8_t *base)
{
	/* profile */
	tm_idx_profile_t **p = mi->profile.arr;
	size_t const pcnt = mi->profile.cnt;
	for(size_t i = 0; i < pcnt; i++) {
		p[i] = _add_offset(p[i], base);
		p[i]->meta.name = tm_idx_profile_name_ptr(p[i]);
		_unused(tm_idx_finalize_profile(p[i]));		/* never fail; instanciate dz */
	}

	/* sketch */
	tm_idx_sketch_t **s = mi->sketch.arr;
	size_t const scnt = mi->sketch.cnt;
	for(size_t i = 0; i < scnt; i++) {
		s[i] = _add_offset(s[i], base);
	}
	return;
}

// static _force_inline
tm_idx_t *tm_idx_load(void *fp, read_t rfp)
{
	tm_idx_hdr_t hdr;
	size_t const hdr_size = rfp(fp, &hdr, sizeof(tm_idx_hdr_t));
	if(hdr_size != sizeof(tm_idx_hdr_t)) {
		if(hdr_size > 0) { error("failed to read header"); }
		return(NULL);
	}

	if(hdr.magic != TM_IDX_MAGIC) {
		uint64_t base = ((hdr.magic ^ TM_IDX_MAGIC) & 0xffffff) == 0;
		error("invalid header%s. please rebuild the index.", base ? " (possibly index format is old)" : "");
		return(NULL);
	}

	/* bulk load */
	size_t idx_size = _roundup(hdr.size, ARCH_HUGEPAGE_SIZE);
	uint8_t *base = aligned_malloc(idx_size);
	size_t const bytes = rfp(fp, base, hdr.size);

	/* sanity check */
	if(bytes != hdr.size) {
		error("truncated index file. please rebuild the index (expected: %zu, got: %zu).", hdr.size, bytes);
		free(base);
		return(NULL);
	}

	/* restore pointers */
	tm_idx_t *mi = tm_idx_slice_root(base, hdr.size);
	mi->filename.ref     = _add_offset(mi->filename.ref, base);
	mi->filename.profile = _add_offset(mi->filename.profile, base);
	tm_idx_patch_ptr(mi, base);
	return(mi);
}

/**
 * end of index.c
 */

