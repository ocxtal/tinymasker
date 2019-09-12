
/**
 * @file aligner.c
 * @brief (fasta, fasta) -> paf API and multithreading
 *
 * @author Hajime Suzuki
 * @license MIT
 */


#ifndef UNITTEST
#  define UNITTEST				( 1 )
#endif
#define UNITTEST_UNIQUE_ID		2


#include "utils/utils.h"		/* include all */
#include "tinymasker.h"
unittest_config( .name = "aligner" );


#include "dozeu.h"
#include "dbg.h"
#include "index.h"		/* index wrapper */
#include "align.h"
#include "bseq.h"		/* FASTA/Q parser */


/* forward declarations */
#include "aligner.h"


/* printer */
typedef struct {
	struct {
		char const *ptr;
		size_t len;				/* converted to int when passed to printf */
	} name;
	struct {
		size_t len, spos, epos;
	} seq;
} tm_print_seq_t;


// static _force_inline
void tm_print_destory_static(tm_print_t *self)
{
	_unused(self);

	return;
}

// static _force_inline
void tm_print_init_static(tm_print_t *self, tm_print_conf_t const *conf, char const *args)
{
	_unused(args);

	/* copy flip flag */
	self->flip = conf->flip != 0;		/* make sure self->flip is either 1 or 0 */
	return;
}

#if 0
static _force_inline
void tm_print_cigar_forward(tm_print_t *self, uint8_t const *path, size_t len)
{
	_unused(self);

	debug("path(%p), len(%zu)", path, len);

	uint8_t const *p = path, *t = &path[len];
	v16i8_t v = _set_v16i8(*p);
	_print_v16i8(v);

	while(p < t) {
		uint8_t const *q = p + 1;
		while(q < t) {
			v16i8_t const w = _loadu_v16i8(q);
			v16i8_t const eq = _eq_v16i8(v, w);
			_print_v16i8(w);
			_print_v16i8(eq);

			uint64_t const mask = ((v16_masku_t){ .mask = _mask_v16i8(eq) }).all;
			ZCNT_RESULT size_t raw = _tzc_u64(~mask);
			size_t const cnt = MIN2(raw, (size_t)(t - q));
			debug("cnt(%zu)", cnt);

			q += cnt;
			if(cnt < 16) { break; }
		}
		printf("%zu%c", (size_t)(q - p), *p);

		p = q;
		v = _set_v16i8(*p);
	}
	return;
}
#endif

static _force_inline
void tm_print_cigar_reverse(tm_print_t *self, uint8_t const *path, size_t len)
{
	_unused(self);

	debug("path(%p), len(%zu), %s", path, len, path);

	uint8_t const *p = &path[len], *t = path;
	while(p > t) {
		uint8_t const *q = p;
		uint8_t const ch = p[-1];
		v16i8_t const v = _set_v16i8(ch);
		// _print_v16i8(v);

		while(q > t) {
			v16i8_t const w = _loadu_v16i8(q - 16);
			v16i8_t const eq = _eq_v16i8(v, w);
			// _print_v16i8(w);
			// _print_v16i8(eq);

			uint64_t const mask = ((v16_masku_t){ .mask = _mask_v16i8(eq) }).all;
			uint64_t const shifted = mask<<48;
			ZCNT_RESULT size_t raw = _lzc_u64(~shifted);
			size_t const cnt = MIN2(raw, (size_t)(q - t));
			// debug("cnt(%zu)", cnt);

			q -= cnt;
			if(cnt < 16) { break; }
		}

		printf("%zu%c", (size_t)(p - q), ch);
		p = q;
	}
	return;
}

static _force_inline
tm_print_seq_t tm_print_compose_query(bseq_meta_t const *query, tm_aln_t const *aln)
{
	tm_pos_t const pos = tm_canon_pos(aln->pos, aln->span);
	tm_pair_t const span = aln->span;

	tm_print_seq_t const q = {
		.name = {
			.len = bseq_name_len(query),
			.ptr = (char const *)bseq_name(query)
		},
		.seq = {
			.len = bseq_seq_len(query),
			.spos = pos.q,
			.epos = pos.q + span.q - 1
		}
	};
	return(q);
}

static _force_inline
tm_print_seq_t tm_print_compose_ref(tm_idx_sketch_t const *ref, tm_aln_t const *aln)
{
	tm_pos_t const pos = tm_canon_pos(aln->pos, aln->span);
	tm_pair_t const span = aln->span;

	tm_print_seq_t const r = {
		.name = {
			.len = tm_idx_ref_name_len(ref),
			.ptr = (char const *)tm_idx_ref_name_ptr(ref)
		},
		.seq = {
			.len = tm_idx_ref_seq_len(ref),
			.spos = pos.r,
			.epos = pos.r + span.r - 1,
		}
	};
	return(r);
}

static _force_inline
void tm_print_aln(tm_print_t *self, tm_idx_sketch_t const **si, bseq_meta_t const *query, tm_aln_t const *aln)
{
	_unused(self);

	/* load reference sequence object */
	tm_idx_sketch_t const *ref = si[aln->attr.rid];

	/* compose alignment coordinate object */
	tm_print_seq_t const seq[2] = {
		tm_print_compose_query(query, aln),			/* query side won't change except for alignment span */
		tm_print_compose_ref(ref, aln)				/* reference side needs reloaded every time */
	};

	/* leftside and rightside; query comes first if flip == 0 */
	size_t const flip = self->flip & 0x01;
	tm_print_seq_t const *l = &seq[flip];
	tm_print_seq_t const *r = &seq[flip ^ 0x01];
	tm_aln_stat_t const stat = tm_aln_calc_stat(aln, flip);

	/*
	 * in PAF
	 * qname, qlen, qstart (0-based), qend (0-based, inclusive), strand,
	 * rname, rlen, rstart (0-based), rend (0-based, inclusive), #matches, block len, mapq
	 *
	 * we use size_t for (seq len, pos, span) for future compatibility.
	 */
	printf(
		/* query */ "%.*s\t%zu\t%zu\t%zu\t"
		/* dir   */ "%c\t"
		/* ref   */ "%.*s\t%zu\t%zu\t%zu\t"
		/* stats */ "*\t%u\t255\tAS:i:%u\tXS:i:%u\tXI:f:%0.4f\tCG:Z:",	/* and cigar */
		(int)l->name.len, l->name.ptr, l->seq.len, l->seq.spos, l->seq.epos,
		aln->pos.dir ? '-' : '+',
		(int)r->name.len, r->name.ptr, r->seq.len, r->seq.spos, r->seq.epos,
		aln->span.r,
		(uint32_t)aln->score.patched,		/* patched score */
		(uint32_t)aln->score.raw,			/* raw score (with bonus) */
		(double)stat.identity
	);
	tm_print_cigar_reverse(self, aln->path.ptr, aln->path.len);
	printf("\n");
	return;
}


/* sanity checker for sequences read from file */

static _force_inline
uint64_t tm_validate_query(uint8_t const *query, size_t qlen)
{
	/*
	 * we only allow { 0x0, 0x4, 0x8, 0xc } for query-side input.
	 * non-zero when valid, zero otherwise.
	 */
	v32i8_t const mv = _set_v32i8(0xf3);
	v32i8_t cv = _zero_v32i8();			/* accumulator */

	/* we don't have margin at the tail */
	size_t qpos = 0;
	while((qpos += 32) < qlen) {
		v32i8_t const v = _loadu_v32i8(&query[qpos - 32]);
		v32i8_t const t = _and_v32i8(mv, v);
		cv = _or_v32i8(cv, t);			/* accumulate non-zero */
	}

	/* accumulate tail */ {
		v32i8_t const v = _loadu_v32i8(&query[qpos - 32]);
		v32i8_t const t = _and_v32i8(mv, v);

		/* create mask and clear out the tail */
		static uint8_t const inc[32] __attribute__(( aligned(32) )) = {
			 -1,  -2,  -3,  -4,  -5,  -6,  -7,  -8,  -9, -10, -11, -12, -13, -14, -15, -16,
			-17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30, -31, -32
		};
		v32i8_t const iv = _load_v32i8(inc);
		v32i8_t const rv = _set_v32i8(qlen & 0x1f);		/* remainder length */
		v32i8_t const sv = _add_v32i8(iv, rv);			/* first <rem> elements are non-negative */

		/* apply mask */
		v32i8_t const z = _zero_v32i8();
		v32i8_t const s = _sel_v32i8(sv, z, t);
		cv = _or_v32i8(cv, s);
	}

	/* 1 if all zero (valid), 0 otherwise */
	v32i8_t const eq = _eq_v32i8(cv, _zero_v32i8());
	uint64_t const mask = ((v32_masku_t){ .mask = _mask_v32i8(eq) }).all;
	return(mask == 0xffffffffULL);		/* all the columns are zero */
}



/* multithreaded scan-and-mask */

typedef struct {
	size_t id;
	dz_arena_t *mem;
	bseq_batch_t seq_bin;
} tm_mtscan_batch_t;

typedef struct {
	size_t id;
	tm_mtscan_batch_t *batch;
} mm_mtscan_drain_t;
#define mtscan_drain_cmp(a, b)		( (int64_t)(a).id - (int64_t)(b).id )

struct tm_mtscan_s {
	tm_idx_t const *mi;
	bseq_file_t *fp;
	tm_print_t *printer;
	pt_t *pt;

	/* counters */
	size_t icnt, ocnt;
	kvechq_t(mm_mtscan_drain_t) hq;

	tm_scan_t *scan[];
};

// static _force_inline
tm_mtscan_t *tm_mtscan_init(tm_idx_t const *mi, tm_print_t *printer, pt_t *pt)
{
	size_t const size = sizeof(tm_mtscan_t) + pt_nth(pt) * sizeof(tm_scan_t *);
	tm_mtscan_t *self = malloc(size);
	*self = (tm_mtscan_t){
		.mi = mi,
		.printer = printer,
		.pt = pt
	};

	for(size_t i = 0; i < pt_nth(pt); i++) {
		self->scan[i] = tm_scan_init();
	}
	return(self);
}

// static _force_inline
void tm_mtscan_destroy(tm_mtscan_t *self)
{
	if(self == NULL) { return; }

	for(size_t i = 0; i < pt_nth(self->pt); i++) {
		tm_scan_destroy(self->scan[i]);
	}

	kv_hq_destroy(self->hq);		/* sorter */
	free(self);
	return;
}

static
void *tm_mtscan_source(uint32_t tid, tm_mtscan_t *self)
{
	_unused(tid);

	/* NULL if reached tail */
	if(bseq_is_eof(self->fp)) { return(NULL); }

	/* fetch */
	bseq_batch_t *seq_bin = bseq_read(self->fp);
	if(seq_bin == NULL) { return(NULL); }

	/* init working buffer */
	tm_mtscan_batch_t *batch = _sub_offset(seq_bin, offsetof(tm_mtscan_batch_t, seq_bin));
	batch->id  = self->icnt++;
	batch->mem = dz_arena_init(DZ_MEM_INIT_SIZE);
	return(batch);
}

static
void *tm_mtscan_worker(uint32_t tid, tm_mtscan_t *self, tm_mtscan_batch_t *batch)
{
	/* working buffers */
	tm_scan_t *scan   = self->scan[tid];
	bseq_meta_t *meta = bseq_meta_ptr(&batch->seq_bin);

	/* setup result buffer and shrink memory */
	tm_scan_start_batch(scan, batch->mem);

	/* for each query sequence */
	for(size_t i = 0; i < bseq_meta_cnt(&batch->seq_bin); i++) {
		bseq_meta_t *seq = &meta[i];

		/* print info */
		debugblock({
			fprintf(stderr, "begin, tid(%u), i(%zu), len(%zu), seq(%s)\n", tid, i, bseq_seq_len(seq), bseq_name(seq));
		});

		if(tm_validate_query(bseq_seq(seq), bseq_seq_len(seq)) == 0) {
			error("invalid character found in `%.*s'. skipped.", (int)bseq_name_len(seq), bseq_name(seq));
			continue;
		}
		seq->u.ptr = tm_scan_all(scan, self->mi, bseq_seq(seq), bseq_seq_len(seq));

		/* print info */
		debugblock({
			fprintf(stderr, "done, tid(%u), i(%zu), cnt(%zu)\n", tid, i, seq->u.ptr != NULL ? ((tm_alnv_t const *)seq->u.ptr)->cnt : 0);
			tm_alnv_t const *v = seq->u.ptr;
			debug("v(%p), cnt(%zu)", v, v != NULL ? v->cnt : 0);
		});
	}
	return(batch);
}

static _force_inline
void tm_mtscan_drain_core(tm_mtscan_t *self, tm_mtscan_batch_t *batch)
{
	bseq_meta_t const *meta = bseq_meta_ptr(&batch->seq_bin);					/* query */
	tm_idx_sketch_t const **si = (tm_idx_sketch_t const **)self->mi->sketch.arr;		/* ref */

	for(size_t i = 0; i < bseq_meta_cnt(&batch->seq_bin); i++) {
		bseq_meta_t const *seq = &meta[i];
		tm_alnv_t const *v = seq->u.ptr;
		debug("v(%p), cnt(%zu)", v, v != NULL ? v->cnt : 0);

		if(v == NULL) { continue; }

		for(size_t j = 0; j < v->cnt; j++) {
			debug("idx(%zu)", j);
			tm_print_aln(self->printer, si, seq, &v->arr[j]);
		}
	}

	dz_arena_destroy(batch->mem);
	bseq_free(&batch->seq_bin);
	return;
}

static
void tm_mtscan_drain(uint32_t tid, tm_mtscan_t *self, tm_mtscan_batch_t *batch)
{
	_unused(tid);

	/* push batch to heapqueue for sorting */
	kv_hq_push(mm_mtscan_drain_t, mtscan_drain_cmp, self->hq,
		((mm_mtscan_drain_t){
			.id = batch->id,
			.batch = batch
		})
	);

	/* pop */
	while(kv_hq_cnt(self->hq) > 0 && kv_hq_head(self->hq).id == self->ocnt) {
		self->ocnt++;
		batch = kv_hq_pop(mm_mtscan_drain_t, mtscan_drain_cmp, self->hq).batch;
		tm_mtscan_drain_core(self, batch);
	}
	return;
}


// static _force_inline
int tm_mtscan_file(tm_mtscan_t *self, char const *fn)
{
	bseq_conf_t const bseq_conf = {
		.batch_size  = TM_BATCH_SIZE,
		.head_margin = offsetof(tm_mtscan_batch_t, seq_bin),

		/* use 2bit encoding for query */
		.conv_table  = {
			[A] = nA,
			[C] = nC,
			[G] = nG,
			[T] = nT
		}
	};
	self->fp = bseq_open(&bseq_conf, fn);
	if(self->fp == NULL) { return(1); }

	/* multithreaded scan-and-mask */
	kv_hq_clear(self->hq);
	pt_stream(self->pt, self,
		(pt_source_t)tm_mtscan_source,
		(pt_worker_t)tm_mtscan_worker,
		(pt_drain_t)tm_mtscan_drain		/* sort */
	);

	/* done */
	bseq_close(self->fp);
	self->fp = NULL;
	return(0);
}


/**
 * end of aligner.c
 */
