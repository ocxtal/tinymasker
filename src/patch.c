
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


#include "align.h"


/* forward decls */
#include "patch.h"


/* paf integer field index -> tm_aln_t field offset mapping */
#define tm_patch_idx_to_ofs(_idx)		( 2 * (_idx) - 8 + ((_idx) < 4 ? 11 : 0) )

_static_assert(offsetof(tm_aln_t, len.q)  == sizeof(uint32_t) * tm_patch_idx_to_ofs(1));
_static_assert(offsetof(tm_aln_t, spos.q) == sizeof(uint32_t) * tm_patch_idx_to_ofs(2));
_static_assert(offsetof(tm_aln_t, epos.q) == sizeof(uint32_t) * tm_patch_idx_to_ofs(3));
_static_assert(offsetof(tm_aln_t, len.r)  == sizeof(uint32_t) * tm_patch_idx_to_ofs(6));
_static_assert(offsetof(tm_aln_t, spos.r) == sizeof(uint32_t) * tm_patch_idx_to_ofs(7));
_static_assert(offsetof(tm_aln_t, epos.r) == sizeof(uint32_t) * tm_patch_idx_to_ofs(8));



typedef struct {
	;
} tm_mask_t;

typedef struct {
	kvec_t(tm_mask_t) mask;
} tm_patch_tbuf_t;


static _force_inline
uint64_t tm_patch_sort_aln(tm_aln_t **aln, size_t cnt)
{
	return(0);
}

static _force_inline
uint64_t tm_patch_merge_aln(tm_aln_t **aln, size_t cnt)
{
	return(0);
}


/* context */
typedef struct {
	char const *ptr;
	size_t len;
} tm_paf_parse_t;


/* returns nonzero on error */
typedef uint64_t (*tm_paf_callback_t)(tm_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin);

static
uint64_t tm_paf_parse_nop(tm_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
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
uint64_t tm_paf_parse_name(tm_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
{
	uint64_t const x = field == 0;
	uint64_t const y = field == 5;
	if((x + y) != 1 || slen == 0) { return(1); }

	char const **dst = &aln->name.r;
	dst[field == 0] = mm_bin_strdup(bin, str, slen);
	return(0);
}

static
uint64_t tm_paf_parse_int(tm_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
{
	int64_t const n  = mm_atoi(str, slen);
	size_t const ofs = tm_patch_idx_to_ofs(field);

	/* returns nonzero on parsing error */
	if(n == INT64_MIN) { return(1); }

	/* from paf field index to tm_aln_t field */
	uint32_t *dst = &((uint32_t *)aln)[ofs];

	/* read-modify-write to conserve dst flag in rpos field */
	uint32_t const x = _loadu_u32(dst);
	_storeu_u32(dst, x + n);		/* 32bit */
	return(0);
}

static
uint64_t tm_paf_parse_dir(tm_aln_t *aln, size_t field, char const *str, size_t slen, mm_bin_t *bin)
{
	/* either x or y is 1 (exclusive) */
	uint64_t const x = str[0] == '-';
	uint64_t const y = str[0] == '+';

	/* assert(x + y == 1) and assert(slen == 1) (FIXME) */
	if((x + y) != 1 || slen != 1) { return(1); }

	/* save dir */
	aln->spos.dir = x;
	return(0);
}

static _force_inline
tm_paf_parse_t tm_paf_parse_body(tm_aln_t *aln, char const *line, size_t llen, mm_bin_t *bin)
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
tm_paf_parse_t tm_paf_parse_opt(tm_aln_t *aln, char const *line, size_t llen, mm_bin_t *bin)
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
uint64_t tm_paf_to_aln(tm_aln_t *aln, char const *line, size_t llen, mm_bin_t *bin)
{
	/* working buffer and context */
	tm_paf_parse_t w = {
		.ptr = line,
		.len = llen
	};
	memset(aln, 0, sizeof(tm_aln_t));

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
	kvecm_t(tm_aln_t) arr;
	kvm_init(arr, sizeof(tm_alnv_t), 0);

	mm_split_foreach(str, slen, delim, {
		tm_aln_t aln;				/* working buffer; cleared inside */

		if(tm_paf_to_aln(&aln, p, l, bin)) {
			error("at line %zu.", i);
			break;
		}
		kvm_push(arr, aln);
	});

	/* create header */
	tm_alnv_t *alnv = kvm_base_ptr(arr);
	alnv->size = kvm_cnt(arr);
	return(alnv);
}

static _force_inline
size_t tm_aln_to_paf(tm_aln_t const *aln)
{
	return;
}



static _force_inline
void tm_print_aln(tm_print_t *self, tm_idx_sketch_t const **si, bseq_meta_t const *query, tm_aln_t const *aln)
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
uint64_t tm_patch_seq(bseq_seq_t *seq, tm_aln_t const **aln, size_t cnt)
{
	/* we expect aln be sorted by qpos */

	return(0);
}



/**
 * end of patch.c
 */
