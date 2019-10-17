
/**
 * @file index.h
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#ifndef _INDEX_H_INCLUDED
#define _INDEX_H_INCLUDED


#ifndef _DZ_H_INCLUDED
#  error "#include \"dozeu.h\" must be before #include \"index.h\""
#endif
#ifndef _TM_H_INCLUDED
#  error "#include \"tinymasker.h\" must be before #include \"index.h\""
#endif
#ifndef _DBG_H_INCLUDED
#  error "#include \"dbg.h\" must be before #include \"index.h\""
#endif


/* configurations: make user-provided value different from default one */
#define tm_idx_wrap(x)					( (x) + 1 )
#define tm_idx_unwrap(x)				( (x) - 1 )
#define tm_idx_is_default(x)			( (x) == 0 )

typedef struct {
	uint64_t ignore_overflow;

	/* fallback parameters */
	size_t kmer, window;		/* k-mer length and chain window size */
	size_t min_scnt;			/* minimum seed count for chain */
	size_t span_thresh;			/* 0 to skip filter */

	/* extension params */
	uint64_t match, mismatch;
	uint64_t gap_open, gap_extend;
	uint64_t max_gap_len;
	uint64_t full_length_bonus;			/* anchoring bonus */
	int64_t min_score;
	uint64_t use_raw_score;
} tm_idx_conf_t;


/* matcher object; k-mer length, score matrix, and score thresholds */
typedef struct {
	uint32_t size;

	/* k-mer size */
	uint32_t kbits;

	/* chaining */
	struct {
		union {
			struct { uint32_t u, v; } sep;
			uint64_t all;
		} window;

		/* the following four must be in this order */
		uint32_t kadj[2];
		uint32_t min_scnt;
		uint32_t squash_intv;
	} chain;

	/* filtering */
	struct {
		uint8_t score_matrix[16];	/* match-mismatch score matrix */
		uint8_t gap[16];			/* gap */
		uint8_t init[16];			/* p = 0 */
		uint32_t span_thresh;		/* -1 to disable */
		int32_t min_score;			/* discard if sum of extension scores is smaller than min_score */
	} filter;

	/* metadata */
	struct {
		char *name;
	} meta;

	/* stats */
	struct {
		double rlambda;				/* 1 / lambda */
	} stat;

	/* extension */
	struct {
		dz_profile_t *dz;

		/* Smith-Waterman params */
		int32_t min_score;
		uint16_t use_raw;
		uint16_t bonus;				/* reference-side bonus; handled manually */
		uint16_t glim, unused;		/* max_ins_len and max_del_len */
		uint8_t giv, gev, gih, geh;
		int8_t score_matrix[DZ_QUERY_MAT_SIZE * DZ_REF_MAT_SIZE];
	} extend;
} tm_idx_profile_t;
_static_assert((sizeof(tm_idx_profile_t) % 16) == 0);


/* profile construction */
tm_idx_profile_t *tm_idx_build_profile(tm_idx_conf_t const *conf);
void tm_idx_destroy_profile(tm_idx_profile_t *profile);



/* simple match-mismatch score parameter, converted from tm_idx_profile_t */
typedef struct {
	uint32_t match, mismatch;
} tm_idx_score_t;

tm_idx_score_t tm_idx_get_score(tm_idx_profile_t const *profile);



/* k-mer index object (wrapper of ref_sketch_t) */
typedef struct {
	uint32_t size;					/* entire object size */
	uint32_t rid, pid;				/* reference id and profile id */
	uint32_t sofs;					/* sequence base position (sequence length is saved inside tm_ref_sketch_t) */
} tm_idx_sketch_hdr_t;

typedef struct {
	tm_idx_sketch_hdr_t h;
	tm_ref_sketch_t ref;			/* wraps ref_sketch_t */
} tm_idx_sketch_t;

#define tm_idx_sketch(x)			( (tm_ref_sketch_t const *)&(x)->ref )
#define tm_idx_ref_rid(x)			( (x)->h.rid )
#define tm_idx_ref_pid(x)			( (x)->h.pid )
#define tm_idx_ref_seq_ptr(x)		( (uint8_t const *)_add_offset((x), (x)->h.sofs) )
#define tm_idx_ref_seq_len(x)		( (x)->ref.slen )
#define tm_idx_ref_name_ptr(x)		( (uint8_t const *)_add_offset((x), sizeof(tm_idx_sketch_hdr_t) + (x)->ref.size) )
#define tm_idx_ref_name_len(x)		( strlen((char const *)tm_idx_ref_name_ptr(x)) )


typedef struct {
	/* save global configuration at the head */
	void *base;
	struct {
		char const *ref;
		char const *profile;
	} filename;

	/* score matrices */
	struct {
		tm_idx_profile_t **arr;
		size_t cnt;
	} profile;

	/* sketch array */
	struct {
		tm_idx_sketch_t **arr;
		size_t cnt;
	} sketch;
} tm_idx_t;


/* index construction */
tm_idx_t *tm_idx_gen(tm_idx_conf_t const *conf, pt_t *pt, char const *ref, char const *profile, FILE *lfp);
void tm_idx_destroy(tm_idx_t *mi);


/* I/O */
size_t tm_idx_dump(tm_idx_t const *mi, void *fp, write_t wfp);
tm_idx_t *tm_idx_load(void *fp, read_t rfp);



#endif	/* _INDEX_H_INCLUDED */
/**
 * end of index.h
 */
