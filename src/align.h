
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
#ifndef __H_INCLUDED
#  error "#include \"baln.h\" must be before #include \"align.h\""
#endif
#ifndef _INDEX_H_INCLUDED
#  error "#include \"index.h\" must be before #include \"align.h\""
#endif


/* working buffer */
typedef struct tm_scan_s tm_scan_t;

tm_scan_t *tm_scan_init(void);
void tm_scan_destroy(tm_scan_t *self);
void tm_scan_start_batch(tm_scan_t *self, dz_arena_t *mem);


/* alignment calc */
baln_alnv_t *tm_scan_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *seq, size_t slen);


#if 0
/* utils */
typedef struct {
	double identity;
} tm_aln_stat_t;

tm_aln_stat_t tm_aln_calc_stat(tm_aln_t const *aln, uint64_t flip);
#endif

#endif	/* _ALIGN_H_INCLUDED */
/**
 * end of align.h
 */
