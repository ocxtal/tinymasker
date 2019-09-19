
/**
 * @file baln.h
 * @brief general alignment object, and PAF parser-formatter
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#include "utils/utils.h"		/* include all */


/* general alignment and alignment array */
typedef struct {
	size_t r, q;			/* we might need more bits because sequences can be longer */
} baln_pair_t;

typedef struct {
	int32_t raw, patched;	/* need 64bits? */
} baln_score_t;

typedef struct {
	/* sequence by name */
	struct {
		char const *r;
		char const *q;
	} name;

	/* lengths; for paf compatibility */
	baln_pair_t len;

	/* positions */
	baln_pair_t pos;		/* 0-based inclusive */
	baln_pair_t span;		/* always positive for canonical alignment representation */

	/* orientation */
	uint32_t dir;			/* 0 for forward, 1 for reverse */
	uint32_t unused;

	/* scores */
	baln_score_t score;

	/* alignment path */
	struct {
		uint8_t const *ptr;	/* MMMMMIMMDMMMMM...; no discrimination between match and mismatch */
		size_t len;
	} path;
} baln_aln_t;
_static_assert(sizeof(baln_aln_t) == 96);


/* array of alignments (and builder) */
typedef void (*baln_free_t)(void *ptr);

static
void baln_free(void *ctx, void *ptr)
{
	_unused(ctx);
	free(ptr);
	return;
}

typedef struct {
	size_t cnt;
	void *bin;				/* string bin */

	struct {
		baln_free_t body;
		baln_free_t bin;
	} free;
	baln_aln_t arr[];
} baln_alnv_t;

static _force_inline
void baln_alnv_destroy(baln_alnv_t *alnv)
{
	if(alnv->free.bin != NULL) {
		alnv->free.bin(alnv->bin);
	}

	if(alnv->free == NULL) { return; }
	alnv->free.body(alnv);
	return;
}


typedef struct {
	kvecm_t(baln_aln_t) buf;
} baln_builder_t;

static _force_inline
void baln_builder_init(baln_builder_t *self)
{
	kvm_init(self->buf, sizeof(baln_alnv_t), 0);
	return;
}

static _force_inline
void baln_builder_push(baln_builder_t *self, baln_aln_t const *aln)
{
	kvm_push(self->buf, *aln);
	return;
}

static _force_inline
baln_alnv_t *baln_builder_finalize(baln_builder_t *self)
{
	baln_alnv_t *alnv = kvm_base_ptr(self->buf);

	alnv->cnt  = kvm_cnt(self->buf);
	alnv->free = baln_free;
	return(alnv);
}


/* paf integer field index -> baln_aln_t field offset mapping */
#define tm_paf_idx_to_ofs(_idx)		( 4 * (_idx) - 20 + ((_idx) < 4 ? 22 : 0) )

_static_assert(offsetof(baln_aln_t, len.q)  == sizeof(uint32_t) * tm_paf_idx_to_ofs(1));	/* 6 */
_static_assert(offsetof(baln_aln_t, spos.q) == sizeof(uint32_t) * tm_paf_idx_to_ofs(2));	/* 10 */
_static_assert(offsetof(baln_aln_t, epos.q) == sizeof(uint32_t) * tm_paf_idx_to_ofs(3));	/* 14 */
_static_assert(offsetof(baln_aln_t, len.r)  == sizeof(uint32_t) * tm_paf_idx_to_ofs(6));	/* 4 */
_static_assert(offsetof(baln_aln_t, spos.r) == sizeof(uint32_t) * tm_paf_idx_to_ofs(7));	/* 8 */
_static_assert(offsetof(baln_aln_t, epos.r) == sizeof(uint32_t) * tm_paf_idx_to_ofs(8));	/* 12 */


/* context */
typedef struct {
	char const *ptr;
	size_t len;
} tm_paf_parse_t;


/* returns nonzero on error */
typedef uint64_t (*tm_paf_callback_t)(baln_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin);

static
uint64_t tm_paf_parse_nop(baln_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
{
	/* do nothing */
	_unused(aln);
	_unused(field);
	_unused(str);
	_unused(slen);
	_unused(bin);
	return(0);
}

static
uint64_t tm_paf_parse_name(baln_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
{
	uint64_t const x = field == 0;
	uint64_t const y = field == 5;
	if((x + y) != 1 || slen == 0) { return(1); }

	char const **dst = &aln->name.r;
	dst[field == 0] = mm_bin_strdup(bin, str, slen);
	return(0);
}

static
uint64_t tm_paf_parse_int(baln_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
{
	int64_t const n  = mm_atoi(str, slen);
	size_t const ofs = tm_paf_idx_to_ofs(field);

	/* returns nonzero on parsing error */
	if(n == INT64_MIN) { return(1); }

	/* from paf field index to baln_aln_t field */
	size_t *dst = &((size_t *)aln)[ofs];
	_storeu_u64(dst, n);
	return(0);
}

static
uint64_t tm_paf_parse_dir(baln_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
{
	/* either x or y is 1 (exclusive) */
	uint64_t const x = str[0] == '-';
	uint64_t const y = str[0] == '+';

	/* assert(x + y == 1) and assert(slen == 1) (FIXME) */
	if((x + y) != 1 || slen != 1) { return(1); }

	/* save dir */
	aln->dir = x;
	return(0);
}

static _force_inline
tm_paf_parse_t tm_paf_parse_body(baln_aln_t *aln, char const *line, size_t llen, mm_bin_t *bin)
{
	/* parser for each field */
	static tm_paf_callback_t const parser[16] = {
		/* qname, qlen, qstart, qend */
		tm_paf_parse_name, tm_paf_parse_int, tm_paf_parse_int, tm_paf_parse_int,
		/* direction */
		tm_paf_parse_dir,
		/* rname, rlen, rstart, rend */
		tm_paf_parse_name, tm_paf_parse_int, tm_paf_parse_int, tm_paf_parse_int,
		/* score and miscellaneous */
		tm_paf_parse_nop, tm_paf_parse_int, tm_paf_parse_nop
	};

	/* delimiters */
	static uint8_t const delim[16] __attribute__(( aligned(16) )) = {
		'\t', '\v', '\r', '\n', '\0'		/* space ' ' is skipped */
	};

	mm_split_foreach(line, llen, delim, {	/* parse-and-fill loop */
		if(parser[i] == NULL) {
			mm_split_break(line, llen);		/* update pointer and length for optional fields */
			break;
		}
		if(parser[i](aln, i, p, l, bin)) {
			error("unparsable PAF element `%.*s' (index: %zu).", (int)l, p, i);
			return((tm_paf_parse_t){
				.ptr = NULL,
				.len = 0
			});
		}
	});

	/* done without error */
	return((tm_paf_parse_t){
		.ptr = line,
		.len = llen
	});
}

static _force_inline
tm_paf_parse_t tm_paf_parse_opt(baln_aln_t *aln, char const *line, size_t llen, mm_bin_t *bin)
{
	/* needs modified if collide */
	#define HASH_SIZE			( 16 )
	#define _key(_str)			( (uint16_t)((_str)[0]) | (((uint16_t)((_str)[1]))<<8) )
	#define _hash(_k)			( (3 * (_k) + ((_k)>>8)) & (HASH_SIZE - 1) )
	#define _elem(_str, _fp)	[_hash(_key(_str))] = { .key = _key(_str) .fp = (_fp) }

	static struct { uint16_t key; tm_paf_callback_t fp; } const parser[] = {
		_elem("AS", tm_paf_parse_score),
		_elem("XS", tm_paf_parse_score),
		_elem("XI", NULL),
		_elem("CG", tm_paf_parse_cigar)
	};

	/* delimiters */
	static uint8_t const delim[16] __attribute__(( aligned(16) )) = {
		'\t', '\v', '\r', '\n', '\0'		/* space ' ' is skipped */
	};

	mm_split_foreach(line, llen, delim, {
		uint16_t const key = _key(p);
		size_t const hash  = _hash(key);

		/* skip unknown record */
		if(parser[hash].fp == NULL || parser[hash].key != key) { continue; }

		/* 0 if successful */
		if(parser[hash].fp(aln, key, p, l, bin)) {
			error("unparsable PAF optional record `%.*s'.", (int)l, p);
			return((tm_paf_parse_t){
				.ptr = NULL,
				.len = 0
			});
		}
	});

	/* llen == 0 if everything is done */
	return((tm_paf_parse_t){
		.ptr = line,
		.len = llen
	});

	#undef HASH_SIZE
	#undef _hash
	#undef _elem
}

static _force_inline
uint64_t tm_paf_to_aln(baln_aln_t *aln, char const *line, size_t llen, mm_bin_t *bin)
{
	/* working buffer and context */
	tm_paf_parse_t w = {
		.ptr = line,
		.len = llen
	};
	memset(aln, 0, sizeof(baln_aln_t));

	/* body */
	w = tm_paf_parse_body(aln, w.ptr, w.len, bin);
	if(w.ptr == NULL) { return(1); }

	/* optional field */
	w = tm_paf_parse_opt(aln, w.ptr, w.len, bin);
	if(w.ptr == NULL) { return(1); }
	if(w.len != 0) {
		error("unknown error on parsing PAF record `%.*s'.", line, llen);
		return(1);
	}

	/* done without error */
	return(0);
}

static _force_inline
tm_alnv_t *tm_paf_to_alnv(char const *str, size_t slen)
{
	/* delimiters */
	static uint8_t const delim[16] __attribute__(( aligned(16) )) = {
		'\v', '\r', '\n', '\0'		/* space and tab skipped */
	};

	/* init margined array */
	baln_builder_t alnv;			/* alignment bin */
	baln_builder_init(&alnv);

	mm_bin_t *bin = mm_bin_init();	/* string bin */

	mm_split_foreach(str, slen, delim, {
		baln_aln_t aln;				/* working buffer; cleared inside */

		if(tm_paf_to_aln(&aln, p, l, bin)) {
			error("at line %zu.", i);
			break;
		}
		baln_builder_push(&alnv, &aln);
	});

	/* create header */
	tm_alnv_t *alnv = kvm_base_ptr(arr);
	alnv->size = kvm_cnt(arr);
	return(alnv);
}

static _force_inline
size_t baln_aln_to_paf(baln_aln_t const *aln)
{
	return;
}



static _force_inline
void tm_print_aln(tm_print_t *self, tm_idx_sketch_t const **si, bseq_meta_t const *query, baln_aln_t const *aln)
{
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



// static _force_inline
uint64_t tm_paf_seq(bseq_seq_t *seq, baln_aln_t const **aln, size_t cnt)
{
	/* we expect aln be sorted by qpos */

	return(0);
}



/**
 * end of baln.h
 */
