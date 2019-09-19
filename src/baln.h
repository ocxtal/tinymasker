
/**
 * @file baln.h
 * @brief general alignment object, and PAF parser-formatter
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#ifndef _BALN_H_INCLUDED
#define _BALN_H_INCLUDED

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

	/* scores */
	float identity;
	baln_score_t score;

	/* alignment path */
	struct {
		uint8_t const *ptr;	/* MMMMMIMMDMMMMM...; no discrimination between match and mismatch */
		size_t len;
	} path;
} baln_aln_t;
_static_assert(sizeof(baln_aln_t) == 96);


/*
 * PAF record -> baln_aln_t parser
 */

/* strdup; string bin */
typedef char const *(*baln_bin_strdup_t)(void *bin, char const *str, size_t len);

typedef struct {
	void *bin;					/* mm_bin_t* in mmstring.h */
	baln_bin_strdup_t strdup;	/* mm_bin_strdup */
} baln_bin_t;


/* paf integer field index -> baln_aln_t field offset mapping */
#define baln_paf_ofs(_idx)		( 4 * (_idx) - 20 + ((_idx) < 4 ? 22 : 0) )

_static_assert(offsetof(baln_aln_t, len.q)  == sizeof(uint32_t) * baln_paf_ofs(1));	/* 6 */
_static_assert(offsetof(baln_aln_t, spos.q) == sizeof(uint32_t) * baln_paf_ofs(2));	/* 10 */
_static_assert(offsetof(baln_aln_t, epos.q) == sizeof(uint32_t) * baln_paf_ofs(3));	/* 14 */
_static_assert(offsetof(baln_aln_t, len.r)  == sizeof(uint32_t) * baln_paf_ofs(6));	/* 4 */
_static_assert(offsetof(baln_aln_t, spos.r) == sizeof(uint32_t) * baln_paf_ofs(7));	/* 8 */
_static_assert(offsetof(baln_aln_t, epos.r) == sizeof(uint32_t) * baln_paf_ofs(8));	/* 12 */


/* internal context */
typedef struct {
	char const *ptr;
	size_t len;
} baln_paf_parse_t;


/* returns nonzero on error */
typedef uint64_t (*baln_paf_callback_t)(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin);

static
uint64_t baln_paf_parse_nop(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
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
uint64_t baln_paf_parse_name(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
{
	uint64_t const x = field == 0;
	uint64_t const y = field == 5;
	if((x + y) != 1 || slen == 0) { return(1); }

	char const **dst = &aln->name.r;
	dst[field == 0] = mm_bin_strdup(bin, str, slen);
	return(0);
}

static
uint64_t baln_paf_parse_int(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
{
	int64_t const n  = mm_atoi(str, slen);
	size_t const ofs = baln_paf_ofs(field);

	/* returns nonzero on parsing error */
	if(n == INT64_MIN) { return(1); }

	/* from paf field index to baln_aln_t field */
	size_t *dst = &((size_t *)aln)[ofs];
	_storeu_u64(dst, n);
	return(0);
}

static
uint64_t baln_paf_parse_dir(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
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
baln_paf_parse_t baln_paf_parse_body(baln_aln_t *aln, char const *line, size_t llen, baln_bin_t *bin)
{
	/* parser for each field */
	static baln_paf_callback_t const parser[16] = {
		/* qname, qlen, qstart, qend */
		baln_paf_parse_name, baln_paf_parse_int, baln_paf_parse_int, baln_paf_parse_int,
		/* direction */
		baln_paf_parse_dir,
		/* rname, rlen, rstart, rend */
		baln_paf_parse_name, baln_paf_parse_int, baln_paf_parse_int, baln_paf_parse_int,
		/* score and miscellaneous */
		baln_paf_parse_nop, baln_paf_parse_int, baln_paf_parse_nop
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
			return((baln_paf_parse_t){
				.ptr = NULL,
				.len = 0
			});
		}
	});

	/* done without error */
	return((baln_paf_parse_t){
		.ptr = line,
		.len = llen
	});
}

static
uint64_t baln_paf_parse_score(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
{
	return(0);
}

static
uint64_t baln_paf_parse_cigar(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
{
	return(0);
}

static _force_inline
baln_paf_parse_t baln_paf_parse_opt(baln_aln_t *aln, char const *line, size_t llen, baln_bin_t *bin)
{
	/* needs modified if collide */
	#define HASH_SIZE			( 16 )
	#define _key(_str)			( (uint16_t)((_str)[0]) | (((uint16_t)((_str)[1]))<<8) )
	#define _hash(_k)			( (3 * (_k) + ((_k)>>8)) & (HASH_SIZE - 1) )
	#define _elem(_str, _fp)	[_hash(_key(_str))] = { .key = _key(_str) .fp = (_fp) }

	static struct { uint16_t key; baln_paf_callback_t fp; } const parser[] = {
		_elem("AS", baln_paf_parse_score),
		_elem("XS", baln_paf_parse_score),
		_elem("XI", NULL),
		_elem("CG", baln_paf_parse_cigar)
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
			return((baln_paf_parse_t){
				.ptr = NULL,
				.len = 0
			});
		}
	});

	/* llen == 0 if everything is done */
	return((baln_paf_parse_t){
		.ptr = line,
		.len = llen
	});

	#undef HASH_SIZE
	#undef _hash
	#undef _elem
}

static _force_inline
uint64_t baln_paf_parse(baln_aln_t *aln, char const *line, size_t llen, baln_bin_t *bin)
{
	/* working buffer and context */
	baln_paf_parse_t w = {
		.ptr = line,
		.len = llen
	};
	memset(aln, 0, sizeof(baln_aln_t));

	/* body */
	w = baln_paf_parse_body(aln, w.ptr, w.len, bin);
	if(w.ptr == NULL) { return(1); }

	/* optional field */
	w = baln_paf_parse_opt(aln, w.ptr, w.len, bin);
	if(w.ptr == NULL) { return(1); }
	if(w.len != 0) {
		error("unknown error on parsing PAF record `%.*s'.", line, llen);
		return(1);
	}

	/* done without error */
	return(0);
}



/*
 * alignment bin: array of alignments (and builder)
 */
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
	void *bin;
	baln_bin_strdup_t strdup;

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
baln_alnv_t *baln_builder_finalize(baln_builder_t *self, void *bin, baln_free_t bin_free)
{
	baln_alnv_t *alnv = kvm_base_ptr(self->buf);

	alnv->cnt = kvm_cnt(self->buf);
	alnv->bin = bin;
	alnv->free.body = baln_free;
	alnv->free.bin  = bin_free;
	return(alnv);
}


/* cached string bin for duplicated names */
typedef struct {
	mm_bin_t *bin;

	struct {
		char const *ptr;
		size_t len;
	} cache;
} baln_cbin_t;


static _force_inline
void baln_cbin_init(baln_cbin_t *self)
{
	self->bin = mm_bin_init();
	self->ptr = NULL;
	self->len = 0;
	return;
}

static _force_inline
mm_bin_t *baln_cbin_finalize(baln_cbin_t *self)
{
	return(self->bin);
}

static _force_inline
void baln_cbin_update(baln_cbin_t *self, char const *str, size_t len)
{
	self->cache.ptr = str;
	self->cache.len = len == 0 ? mm_strlen(str) : len;
	return;
}

static _force_inline
uint64_t baln_cbin_is_dup(baln_cbin_t *self, char const *str, size_t len)
{
	if(self->cache.len != len) {
		return(0);
	}
	if(mm_strncmp(self->cache.ptr, str, MIN2(self->cache.len, len)) != 0) {
		return(0);
	}
	return(1);
}

static
char const *baln_cbin_strdup(baln_cbin_t *self, char const *str, size_t len)
{
	if(baln_cbin_is_dup(self, str, len)) {
		return(self->cache.ptr);
	}
	return(mm_bin_strdup(self->bin, str, len));
}


/* streamed parser */

typedef struct {
	size_t batch_size;
} baln_conf_t;

typedef struct {
	size_t cnt;			/* paf record count */
	uint64_t status;	/* nonzero if broken */
	rbread_t rb;
} baln_file_t;

typedef struct {
	size_t cnt;
	uint64_t status;
} baln_close_t;

static _force_inline
baln_file_t *baln_open(baln_conf_t const *conf, char const *fn)
{
	/* create instance */
	baln_file_t *fp = (baln_file_t *)calloc(1, sizeof(baln_file_t));

	/* open stream */
	size_t const batch_size = conf->batch_size == 0 ? BALN_BATCH_SIZE : conf->batch_size;
	if(rbopen_bulk_static(&fp->rb, fn, batch_size) != 0) {
		error("failed to open file `%s'.", fn);
		goto _baln_open_fail;
	}
	return(fp);

_baln_open_fail:;
	free(fp);
	return(NULL);
}

static _force_inline
baln_close_t baln_close(baln_file_t *fp)
{
	baln_close_t c = { 0 };
	if(fp == NULL) { return(c); }

	c.cnt    = fp->cnt;
	c.status = fp->status;

	rbclose_static(&fp->rb);
	free(fp);
	return(c);
}

static _force_inline
uint64_t baln_read_is_split(baln_cbin_t *bin, char const *line, size_t llen)
{
	/* delimiters */
	static uint8_t const delim[16] __attribute__(( aligned(16) )) = {
		'\t', '\v', '\r', '\n', '\0'		/* space ' ' is skipped */
	};

	/* compare first name, split if the names are different */
	mm_split_foreach(line, llen, delim, {
		/* regarded as complete match for the first element */
		if(bin->cache.ptr != NULL) {
			return(baln_cbin_is_dup(bin, p, l));
		}

		char const *qname = mm_bin_strdup(bin->bin, p, l);
		baln_cbin_update(bin, qname, l);
		return(0);
	});
	return(0);
}

static _force_inline
baln_alnv_t *baln_read(baln_file_t *fp)
{
	/* init string bin */
	baln_cbin_t bin;
	baln_cbin_init(&bin);
	baln_bin_t b = {
		.bin    = &bin,
		.strdup = baln_cbin_strdup
	};

	/* init alignment bin */
	baln_builder_t alnv;
	baln_builder_init(&alnv);

	rb_readline_t r = fp->prev;
	while(r.ptr != NULL) {
		if(baln_read_is_split(&bin, r.ptr, r.len)) { break; }

		fp->cnt++;

		/* parse; using cached bin */
		baln_aln_t w = { 0 };
		if(baln_paf_parse(&w, r.ptr, r.len, &b)) {
			error("at line %zu.", fp->cnt);
			break;
		}

		/* push */
		baln_builder_push(&alnv, &w);

		/* fetch next */
		r = rbreadline(&fp->rb);
	}

	return(baln_builder_finalize(&alnv,
		(void *)baln_cbin_finalize(&bin),
		(baln_free_t)mm_bin_free
	));
}


static _force_inline
tm_alnv_t *baln_paf_to_alnv(char const *str, size_t slen)
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

		if(baln_paf_to_aln(&aln, p, l, bin)) {
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
uint64_t baln_paf_seq(bseq_seq_t *seq, baln_aln_t const **aln, size_t cnt)
{
	/* we expect aln be sorted by qpos */

	return(0);
}

#endif		/* _BALN_H_INCLUDED */

/**
 * end of baln.h
 */
