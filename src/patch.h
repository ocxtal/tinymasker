
/**
 * @file patch.h
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#ifndef _PATCH_H_INCLUDED
#define _PATCH_H_INCLUDED


#include "bseq.h"
#include "baln.h"


/* configurations: make user-provided value different from default one */
#define tm_patch_wrap(x)				( (x) + 1 )
#define tm_patch_unwrap(x)				( (x) - 1 )
#define tm_patch_is_default(x)			( (x) == 0 )


/* configurations */
typedef struct {
	size_t max_gap_len;
	size_t margin_len;
} tm_patch_conf_t;

/* working buffer */
typedef struct {
	struct {
		size_t gap;
		size_t margin;
	} len;

	struct {
		bseq_dump_t dump;
		bseq_conf_t conf;
	} bseq;

	rbread_t *fp;
} tm_patch_tbuf_t;

uint64_t tm_patch_tbuf_init_static(tm_patch_tbuf_t *w, tm_patch_conf_t const *conf, char const *aln_file);
uint64_t tm_patch_tbuf_destroy_static(tm_patch_tbuf_t *w);


/* patcher */
uint64_t tm_patch_seq(bseq_meta_t *seq, baln_aln_t const *aln, size_t cnt);
uint64_t tm_patch_file(tm_patch_tbuf_t *w, char const *seq_file, pt_t *pt);


#endif	/* _PATCH_H_INCLUDED */
/**
 * end of patch.h
 */
