
/**
 * @file align.h
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#ifndef _ALIGN_H_INCLUDED
#define _ALIGN_H_INCLUDED

#ifndef _DZ_H_INCLUDED
#  error "#include \"dozeu.h\" must be before #include \"align.h\""
#endif
#ifndef _TM_H_INCLUDED
#  error "#include \"tinymasker.h\" must be before #include \"align.h\""
#endif
#ifndef _INDEX_H_INCLUDED
#  error "#include \"index.h\" must be before #include \"align.h\""
#endif



/* general position */
typedef struct {
	uint32_t r, q;
} tm_pair_t;

static _force_inline
tm_pair_t tm_add_pair(tm_pair_t a, tm_pair_t b)
{
	return((tm_pair_t){
		.r = a.r + b.r,
		.q = a.q + b.q
	});
}



/* position with direction */
typedef struct {
	uint16_t r;
	uint8_t dir, unused;		/* dir == 1 for reverse, 0 for forward */
	uint32_t q;
} tm_pos_t;

static _force_inline
char *tm_pos_to_str(tm_pos_t const *p)
{
	return(xbprintf("dir(%u), pos(%u, %u), unused(%u)", p->dir, p->r, p->q, p->unused));
}

static _force_inline
tm_pair_t tm_pos_to_pair(tm_pos_t pos)
{
	/* drop dir and unused */
	return((tm_pair_t){
		.r = pos.r,
		.q = pos.q
	});
}

static _force_inline
tm_pos_t tm_add_pos(tm_pos_t a, tm_pair_t b)
{
	return((tm_pos_t){
		.r   = a.r + (a.dir ? -b.r : b.r),
		.q   = a.q + b.q,
		.dir = a.dir,
		.unused = a.unused
	});
}

static _force_inline
tm_pos_t tm_sub_pos(tm_pos_t a, tm_pair_t b)
{
	return((tm_pos_t){
		.r   = a.r + (a.dir ? b.r : -b.r),
		.q   = a.q - b.q,
		.dir = a.dir,
		.unused = a.unused
	});
}

static _force_inline
tm_pos_t tm_canon_pos(tm_pos_t pos, tm_pair_t span)
{
	return((tm_pos_t){
		.r   = pos.r - (pos.dir ? span.r : 0),	/* convert to head position */
		.q   = pos.q,
		.dir = pos.dir,
		.unused = pos.unused
	});
}



/* alignment scores */
typedef struct {
	int32_t raw;			/* raw SW score + end bonus */
	int32_t patched;		/* raw SW score + end bonus + adjustment */
} tm_score_t;


/* alignment result */
typedef struct {
	/* rbt */
	rbt_header_t h;
	uint32_t qmax;

	/* positions */
	struct {
		uint32_t rid : 24;
		uint32_t max_weight : 8;
	} attr;

	tm_pos_t pos;
	tm_pair_t span;

	/* scores */
	tm_score_t score;

	/* alignment path */
	struct {
		uint8_t const *ptr;
		size_t len;
	} path;

	/* everything else in dz_alignment_t */
	dz_alignment_t const *aln;
} tm_aln_t;
_static_assert(sizeof(tm_aln_t) == 64);




/* working buffer */
typedef struct tm_scan_s tm_scan_t;

tm_scan_t *tm_scan_init(void);
void tm_scan_destroy(tm_scan_t *self);
void tm_scan_start_batch(tm_scan_t *self, dz_arena_t *mem);


/* alignment calc */
tm_alnv_t *tm_scan_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *seq, size_t slen);


/* utils */
typedef struct {
	double identity;
} tm_aln_stat_t;

tm_aln_stat_t tm_aln_calc_stat(tm_aln_t const *aln, uint64_t flip);


#endif	/* _ALIGN_H_INCLUDED */
/**
 * end of align.h
 */
