
/**
 * @file patch.c
 * @brief convert masked FASTA/Q bases to lowercases
 *
 * @author Hajime Suzuki
 * @license MIT
 */


#ifndef UNITTEST
#  define UNITTEST				( 1 )
#endif
#define UNITTEST_UNIQUE_ID		6


#include "utils/utils.h"		/* include all */
#include "tinymasker.h"
unittest_config( .name = "patch" );

#include "bseq.h"
#include "baln.h"


/* forward decls */
#include "patch.h"



typedef struct {
	uint64_t unused;
} tm_mask_t;


static _force_inline
uint64_t tm_patch_sort_aln(baln_aln_t const *aln, size_t cnt)
{
	return(0);
}

static _force_inline
uint64_t tm_patch_merge_aln(baln_aln_t const *aln, size_t cnt)
{
	return(0);
}

/* we need working buffer for this function?? */
// static _force_inline
uint64_t tm_patch_seq(bseq_seq_t *seq, baln_aln_t const *aln, size_t cnt)
{
	/* we expect aln be sorted by qpos */

	return(0);
}




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

// static _force_inline
uint64_t tm_patch_tbuf_init_static(tm_patch_tbuf_t *w, tm_patch_conf_t const *conf, char const *aln_file)
{
	memset(w, 0, sizeof(tm_patch_tbuf_t));

	/* constants */
	w->len.gap    = conf->max_gap_len;
	w->len.margin = conf->margin_len;

	/* sequence I/O */
	w->bseq.conf = {
		.batch_size  = TM_BATCH_SIZE,
		.head_margin = 0,
		.conv_table  = {
			[A] = A,
			[C] = C,
			[G] = G,
			[T] = T
		}
	};
	bseq_dump_init_static(&w->bseq.dump, &w->bseq.conf, NULL);

	/* init stream */
	baln_conf_t const baln_conf = {
		.batch_size = TM_BATCH_SIZE
	};
	w->fp = baln_open(&baln_conf, aln_file);
	return(w->fp != NULL);
}

// static _force_inline
uint64_t tm_patch_tbuf_destroy_static(tm_patch_tbuf_t *w)
{
	bseq_dump_destroy_static(&w->bseq.dump);

	baln_close_t const r = baln_close(w->fp);
	_unused(r);			/* FIXME */
	return(0);
}




/* sequence fetcher */
typedef struct {
	bseq_file_t *fp;

	/* current bin */
	size_t idx;
	bseq_batch_t *batch;
} tm_patch_fetch_t;

static _force_inline
uint64_t tm_patch_fetch_init(tm_patch_fetch_t *self, bseq_conf_t const *conf, char const *filename)
{
	/* open query sequence file */
	self->fp = bseq_open(conf, filename);
	if(self->fp == NULL) { return(1); }

	/* fetch first bin */
	self->batch = bseq_read(self->fp);
	if(self->batch == NULL) {
		bseq_close(self->fp);
		return(2);
	}
	return(0);
}

static _force_inline
uint64_t tm_patch_fetch_destroy(tm_patch_fetch_t *self)
{
	bseq_close_t const r = bseq_close(self->fp);
	_unused(r);
	return(0);
}




/* batch patching */
static _force_inline
bseq_meta_t *tm_patch_fetch_seq(tm_patch_fetch_t *self, char const *qname)
{
	do {
		bseq_meta_t *seq = bseq_meta_ptr(self->batch);
		size_t const scnt = bseq_meta_cnt(self->batch);

		while(self->idx < scnt) {
			bseq_meta_t *s = &seq[self->idx];

			/* forward and save index before comparison */
			self->idx++;

			/* if names match, return it */
			if(mm_strcmp(bseq_name(s), qname) == 0) {
				return(s);
			}
		}

		/* not found in this batch, try next */
		bseq_free(self->batch);
	} while((self->batch = bseq_read(self->fp)) != NULL);

	/* not found */
	return(NULL);
}

// static _force_inline
uint64_t tm_patch_file(tm_patch_tbuf_t *w, char const *seq_file, pt_t *pt)
{
	_unused(pt);

	tm_patch_fetch_t fetch;
	if(tm_patch_fetch_init(&fetch, w->bseq.conf, seq_file)) {
		return(1);
	}

	uint64_t error = 2;
	baln_alnv_t *v = NULL;
	while((v = baln_read(w->fp)) != NULL) {
		baln_aln_t const *aln = baln_alnv_ptr(v);
		size_t const acnt = baln_alnv_cnt(v);

		/* get corresponding sequence for the alignment batch */
		bseq_meta_t *seq = tm_patch_fetch_seq(&batch, fp, aln[0].name.q);
		if(seq == NULL) { goto _tm_patch_file_error; }

		/* patch and dump */
		tm_patch_seq(seq, aln, acnt);
		bseq_dump_seq(&w->bseq.dump, seq);

		/* done */
		baln_alnv_destroy(v);
		v = NULL;
	}
	error = 0;
_tm_patch_file_error:;

	/* done */
	baln_alnv_destroy(v);
	tm_patch_fetch_destroy(&fetch);
	return(error);
}


/**
 * end of patch.c
 */
