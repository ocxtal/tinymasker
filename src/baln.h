
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


/* memory management */
#define BALN_BATCH_SIZE				( 2ULL * 1024 * 1024 )


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
	baln_bin_strdup_t dup;		/* mm_bin_strdup */
} baln_bin_t;


/* paf integer field index -> baln_aln_t field offset mapping */
#define baln_paf_ofs(_idx)		( 2 * (_idx) - 10 + ((_idx) < 4 ? 11 : 0) )

_static_assert(offsetof(baln_aln_t, len.q)  == sizeof(size_t) * baln_paf_ofs(1));	/* 3 */
_static_assert(offsetof(baln_aln_t, pos.q)  == sizeof(size_t) * baln_paf_ofs(2));	/* 5 */
_static_assert(offsetof(baln_aln_t, span.q) == sizeof(size_t) * baln_paf_ofs(3));	/* 7 */
_static_assert(offsetof(baln_aln_t, len.r)  == sizeof(size_t) * baln_paf_ofs(6));	/* 2 */
_static_assert(offsetof(baln_aln_t, pos.r)  == sizeof(size_t) * baln_paf_ofs(7));	/* 4 */
_static_assert(offsetof(baln_aln_t, span.r) == sizeof(size_t) * baln_paf_ofs(8));	/* 6 */


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
	debug("save str(%.*s)", (int)slen, str);

	uint64_t const x = field == 0;
	uint64_t const y = field == 5;
	if((x + y) != 1 || slen == 0) { return(1); }

	char const **dst = &aln->name.r;
	dst[field == 0] = bin->dup(bin->bin, str, slen);
	debug("name(%s)", dst[field == 0]);
	return(0);
}

static
uint64_t baln_paf_parse_int(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
{
	_unused(bin);

	debug("parse str(%.*s) to int", (int)slen, str);

	int64_t const n  = mm_atoi(str, slen);
	size_t const ofs = baln_paf_ofs(field);

	/* returns nonzero on parsing error */
	if(n == INT64_MIN) { return(1); }

	/* from paf field index to baln_aln_t field */
	size_t *dst = &((size_t *)aln)[ofs];
	_storeu_u64(dst, n);
	debug("int(%lu)", n);
	return(0);
}

static
uint64_t baln_paf_parse_dir(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
{
	debug("parse str(%.*s) to dir", (int)slen, str);

	_unused(field);
	_unused(bin);

	/* either x or y is 1 (exclusive) */
	uint64_t const x = str[0] == '-';
	uint64_t const y = str[0] == '+';

	/* assert(x + y == 1) and assert(slen == 1) (FIXME) */
	if((x + y) != 1 || slen != 1) { return(1); }

	/* save dir */
	aln->dir = x;
	debug("dir(%lu)", x);
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
		baln_paf_parse_nop, baln_paf_parse_nop, baln_paf_parse_nop
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
	debug("parse str(%.*s) to score", (int)slen, str);

	_unused(aln);
	_unused(field);
	_unused(str);
	_unused(slen);
	_unused(bin);

	return(0);
}

static
uint64_t baln_paf_parse_cigar(baln_aln_t *aln, size_t field, char const *str, size_t slen, baln_bin_t *bin)
{
	debug("parse str(%.*s) to cigar", (int)slen, str);

	_unused(aln);
	_unused(field);
	_unused(str);
	_unused(slen);
	_unused(bin);

	return(0);
}

static _force_inline
baln_paf_parse_t baln_paf_parse_opt(baln_aln_t *aln, char const *line, size_t llen, baln_bin_t *bin)
{
	/* needs modified if collide */
	#define HASH_SIZE			( 16 )
	#define _key(_a, _b)		( (uint16_t)(_a) | ((uint16_t)(_b)<<8) )
	#define _hash(_a, _b)		( (3 * (_a) + (_b)) & (HASH_SIZE - 1) )
	#define _elem(_a, _b, _fp)	[_hash(_a, _b)] = { .key = _key(_a, _b), .fp = (_fp) }

	static struct { uint16_t key; baln_paf_callback_t fp; } const parser[] = {
		_elem('A', 'S', baln_paf_parse_score),
		_elem('X', 'S', baln_paf_parse_score),
		_elem('X', 'I', NULL),
		_elem('C', 'G', baln_paf_parse_cigar)
	};

	/* delimiters */
	static uint8_t const delim[16] __attribute__(( aligned(16) )) = {
		'\t', '\v', '\r', '\n', '\0'		/* space ' ' is skipped */
	};

	mm_split_foreach(line, llen, delim, {
		uint16_t const key = _key(p[0], p[1]);
		size_t const hash  = _hash(p[0], p[1]);

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
		.ptr = &line[llen],
		.len = 0
	});

	#undef HASH_SIZE
	#undef _hash
	#undef _elem
}

static _force_inline
void baln_paf_fixup_pos(baln_aln_t *aln)
{
	debug("fixup, spos(%lu, %lu), epos(%lu, %lu)", aln->pos.q, aln->pos.r, aln->span.q, aln->span.r);

	v2i64_t const spos = _loadu_v2i64(&aln->pos);
	v2i64_t const epos = _loadu_v2i64(&aln->span);
	v2i64_t const ones = _set_v2i64(1);

	/* inclusive -> exclusive end pos */
	v2i64_t const x = _add_v2i64(epos, ones);

	/* epos -> span */
	v2i64_t const y = _sub_v2i64(x, spos);
	_storeu_v2i64(&aln->span, y);

	debug("fixup, spos(%lu, %lu), span(%lu, %lu)", aln->pos.q, aln->pos.r, aln->span.q, aln->span.r);
	return;
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
		error("unknown error on parsing PAF record `%.*s'.", (int)llen, line);
		return(1);
	}

	/* done without error */
	baln_paf_fixup_pos(aln);
	return(0);
}



/*
 * alignment bin: array of alignments (and builder)
 */
typedef void (*baln_free_t)(void *ptr);

static
void baln_free(void *ptr)
{
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

#define baln_alnv_ptr(_x)		( (_x)->arr )
#define baln_alnv_cnt(_x)		( (_x)->cnt )

static _force_inline
void baln_alnv_destroy(baln_alnv_t *alnv)
{
	if(alnv == NULL) { return; }

	if(alnv->free.bin != NULL) {
		alnv->free.bin(alnv->bin);
	}

	if(alnv->free.body == NULL) { return; }
	alnv->free.body(alnv);
	return;
}


typedef struct {
	void *bin;
	baln_bin_strdup_t dup;

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
	kvm_push(baln_aln_t, self->buf, *aln);
	// debug("cnt(%zu)", kvm_cnt(self->buf));
	return;
}

static _force_inline
baln_alnv_t *baln_builder_finalize(baln_builder_t *self, void *bin, baln_free_t bin_free)
{
	size_t const cnt = kvm_cnt(self->buf);
	// debug("cnt(%zu)", cnt);
	if(cnt == 0) {
		bin_free(bin);
		return(NULL);
	}

	baln_alnv_t *alnv = kvm_base_ptr(self->buf);
	alnv->cnt = cnt;
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
	self->bin       = mm_bin_init();
	self->cache.ptr = NULL;
	self->cache.len = 0;
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
	debug("str(%.*s), cache(%.*s)", (int)len, str, (int)self->cache.len, self->cache.ptr);

	if(self->cache.len != len) {
		return(0);
	}
	if(mm_strncmp(self->cache.ptr, str, MIN2(self->cache.len, len)) != 0) {
		return(0);
	}
	debug("matched");
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
	size_t parsed;	/* paf record count */
	uint64_t status;

	/* stream */
	rb_readline_t line;
	rbread_t rb;
} baln_file_t;

typedef struct {
	size_t cnt;
	uint64_t status;
} baln_close_t;

static _force_inline
uint64_t baln_fetch_line(baln_file_t *fp)
{
	do {
		fp->line = rbreadline_intl(&fp->rb);
		if(fp->line.ptr == NULL) { return(0); }
	} while(fp->line.len == 0);
	debug("line fetched, ptr(%p), str(%s)", fp->line.ptr, fp->line.ptr);
	return(1);
}

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

	if(!baln_fetch_line(fp)) {
		error("no valid record found in `%s'.", fn);
		goto _baln_open_fail;
	}
	return(fp);

_baln_open_fail:;
	rbclose_static(&fp->rb);
	free(fp);
	return(NULL);
}

static _force_inline
baln_close_t baln_close(baln_file_t *fp)
{
	baln_close_t c = { 0 };
	if(fp == NULL) { return(c); }

	c.cnt    = fp->parsed;
	c.status = fp->status;

	rbclose_static(&fp->rb);
	free(fp);
	return(c);
}

static _force_inline
uint64_t baln_read_is_split(baln_cbin_t *cbin, char const *line, size_t llen)
{
	/* delimiters */
	static uint8_t const delim[16] __attribute__(( aligned(16) )) = {
		'\t', '\v', '\r', '\n', '\0'		/* space ' ' is skipped */
	};

	/* compare first name, split if the names are different */
	mm_split_foreach(line, llen, delim, {
		/* regarded as complete match for the first element */
		if(cbin->cache.ptr != NULL) {
			debug("test dup");
			return(!baln_cbin_is_dup(cbin, p, l));
		}

		debug("non dup");
		char const *qname = mm_bin_strdup(cbin->bin, p, l);
		baln_cbin_update(cbin, qname, l);
		return(0);
	});
	return(0);
}

static _force_inline
baln_alnv_t *baln_read(baln_file_t *fp)
{
	/* init string bin */
	baln_cbin_t cbin;
	baln_cbin_init(&cbin);
	baln_bin_t b = {
		.bin = (void *)&cbin,
		.dup = (baln_bin_strdup_t)baln_cbin_strdup
	};

	/* init alignment bin */
	baln_builder_t alnv;
	baln_builder_init(&alnv);

	while(fp->line.ptr != NULL) {
		rb_readline_t const r = fp->line;
		if(baln_read_is_split(&cbin, (char const *)r.ptr, r.len)) {
			debug("split detected at line(%zu)", fp->parsed);
			break;
		}

		/* parse; using cached bin */
		baln_aln_t w = { { 0 } };
		if(baln_paf_parse(&w, (char const *)r.ptr, r.len, &b)) {
			fp->status = 1;
			error("at line %zu.", fp->parsed);
			break;
		}

		/* push */
		fp->parsed++;
		baln_builder_push(&alnv, &w);
		debug("parsed, cnt(%zu)", fp->parsed);

		/* fetch next */
		if(!baln_fetch_line(fp)) { break; }
	}

	debug("done");
	return(baln_builder_finalize(&alnv,
		(void *)baln_cbin_finalize(&cbin),
		(baln_free_t)mm_bin_destroy
	));
}


#if 0
static _force_inline
size_t baln_dump_cigar_forward(uint8_t const *path, size_t len)
{
	debug("path(%p), len(%zu)", path, len);

	uint8_t const *p = path, *t = &path[len];
	v16i8_t v = _set_v16i8(*p);
	_print_v16i8(v);

	size_t acc = 0;
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
		acc += (size_t)printf("%zu%c", (size_t)(q - p), *p);

		p = q;
		v = _set_v16i8(*p);
	}
	return(acc);
}
#endif

static _force_inline
size_t baln_dump_cigar_reverse(uint8_t const *path, size_t len)
{
	debug("path(%p), len(%zu), %s", path, len, path);

	uint8_t const *p = &path[len], *t = path;

	size_t acc = 0;
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

		acc += (size_t)printf("%zu%c", (size_t)(p - q), ch);
		p = q;
	}
	return(acc);
}

static _force_inline
size_t baln_dump_aln(baln_aln_t const *a)
{
	size_t acc = 0;

	/*
	 * PAF
	 * qname, qlen, qstart (0-based), qend (0-based, inclusive), strand,
	 * rname, rlen, rstart (0-based), rend (0-based, inclusive), #matches, block len, mapq
	 *
	 * we use size_t for (seq len, pos, span) for future compatibility.
	 */
	acc += (size_t)printf(
		/* query */ "%s\t%zu\t%zu\t%zu\t"
		/* dir   */ "%c\t"
		/* ref   */ "%s\t%zu\t%zu\t%zu\t"
		/* stats */ "*\t%zu\t255\tAS:i:%u\tXS:i:%u\tXI:f:%0.4f\tCG:Z:",	/* and cigar */
		a->name.q, a->len.q, a->pos.q, a->pos.q + a->span.q - 1,
		a->dir ? '-' : '+',
		a->name.r, a->len.r, a->pos.r, a->pos.r + a->span.r - 1,
		a->span.r,
		(uint32_t)a->score.patched,		/* patched score */
		(uint32_t)a->score.raw,			/* raw score (with bonus) */
		(double)a->identity
	);
	acc += baln_dump_cigar_reverse(a->path.ptr, a->path.len);
	acc += (size_t)printf("\n");
	return(acc);
}

static _force_inline
size_t baln_dump_alnv(baln_alnv_t const *alnv)
{
	baln_aln_t const *p = baln_alnv_ptr(alnv);
	size_t const cnt    = baln_alnv_cnt(alnv);

	size_t acc = 0;
	for(size_t i = 0; i < cnt; i++) {
		acc += baln_dump_aln(&p[i]);
	}
	return(acc);
}

#endif		/* _BALN_H_INCLUDED */

/**
 * end of baln.h
 */
