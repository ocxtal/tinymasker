
/**
 * @file tinymasker.c
 * @brief fast repeat masking tool
 *
 * @author Hajime Suzuki
 * @date 2019/2/25
 * @license MIT
 *
 * Note: some lines were derived from minialign, which is originally fork of minimap.
 */

#ifndef BATCH_SIZE
#  define BATCH_SIZE			( 512ULL * 1024 )
#endif

#ifndef MAX_THREADS
#  define MAX_THREADS			( 256 )
#endif


#ifndef UNITTEST
#  define UNITTEST				( 0 )
#endif

#define UNITTEST_UNIQUE_ID		1
#include "utils/utils.h"		/* include all */

#define DZ_NUCL_4BIT
#define DZ_MAT_SIZE				( 16 )
#include "dozeu.h"


/* toml parser and regex */
#include "toml.h"
#include "re.h"


/* misc */
#ifndef MK_VERSION
#  define MK_VERSION			"tinymasker-0.0.1"
#endif
#define MK_ARCH_NAME			ARCH_NAME			/* SSE4.1 or AVX2 */

/* version string */
static _force_inline
char const *tm_version(void)
{
	char const *prefix = "tinymasker-";

	/* remove prefix */
	uint64_t spos = mm_startswith(MK_VERSION, prefix) ? strlen(prefix) : 0;
	return(&MK_VERSION[spos]);
}

unittest_config( .name = "tinymasker" );


/* alphabets */
enum alphabet_iupac {
	A = 0x01,
	C = 0x02,
	G = 0x04,
	T = 0x08,

	/* IUPAC ambiguous bases */
	R = A | G,
	Y = C | T,
	S = G | C,
	W = A | T,
	K = G | T,
	M = A | C,

	B = C | G | T,
	D = A | G | T,
	H = A | C | T,
	V = A | C | G,

	/* invalid */
	N = 0
};

enum alphabet_2bit {
	nA = 0x00,
	nC = 0x01,
	nG = 0x02,
	nT = 0x03
};

/* FASTA/Q parser */
#include "bseq.h"


/* index data structure; we expect each repeat sequence is shorter than 32kbp */
typedef struct {
	uint32_t kbits;
	struct {
		size_t head, tail;
	} margin;
} tm_ref_conf_t;

/* reference index object */
typedef struct {
	uint32_t size, head_margin;		/* object size */
	uint32_t kbits, slen;			/* kmer size and sequence length */
} tm_ref_sketch_t;

static _force_inline
void tm_ref_destroy(tm_ref_sketch_t *sk)
{
	free(_sub_offset(sk, sk->head_margin));
	return;
}


/* construct reference object */
typedef struct {
	uint32_t kmer, pos;
} tm_ref_kpos_t;
_static_assert(sizeof(tm_ref_kpos_t) == 8);

typedef struct { tm_ref_kpos_t *a; size_t n, m; } tm_ref_kpos_v;

#define tm_ref_kmer(x)			( (x).kmer )
KRADIX_SORT_INIT(kmer, tm_ref_kpos_t, tm_ref_kmer, 4);

typedef struct {
	uint32_t f, r;
} tm_ref_kmer_pair_t;

typedef struct {
	/* constants */
	uint32_t amask, shift;

	/* dst array */
	tm_ref_kpos_t *q;

	/* kmer stack */
	uint32_t size;
	uint32_t branches;
	tm_ref_kmer_pair_t kmer[256];	/* branch stack; max ambiguity = 4^4 */
} tm_ref_work_t;

static _force_inline
uint64_t tm_ref_base_to_2bit(uint8_t base)
{
	/* convert iupac nucleotide to pair of 2-bit encoded bases */
	uint64_t const magic = 0x498e4dc399824100;
	return((magic>>(4 * base)) & 0x0f);
}

static _force_inline
uint64_t tm_ref_is_ambiguous(uint8_t base)
{
	uint64_t const magic = 0x1111111011101000;
	return((magic>>(4 * base)) & 0x01);
}

static _force_inline
size_t tm_ref_dup_stack(tm_ref_work_t *w)
{
	size_t prev_size = w->size;

	/* update size and branch pos */
	w->size = w->size<<1;
	w->branches = (w->branches>>2) | (0x01<<w->shift);

	for(size_t i = 0; i < prev_size; i++) {
		w->kmer[w->size + i] = w->kmer[i];
	}
	return(prev_size);
}

static _force_inline
void tm_ref_shrink_stack(tm_ref_work_t *w)
{
	w->size = w->size>>1;
	for(size_t i = 0; i < w->size; i++) {
		w->kmer[i] = w->kmer[2 * i];
	}
	return;
}

static _force_inline
void tm_ref_update_kmer(tm_ref_kmer_pair_t *kmer, size_t size, uint32_t amask, uint32_t shift, uint8_t base)
{
	for(size_t i = 0; i < size; i++) {
		kmer[i].f = ((kmer[i].f<<2) | base) & amask;
		kmer[i].r =  (kmer[i].r>>2) | ((0x03 ^ base)<<shift);
	}
	return;
}

static _force_inline
void tm_ref_push_pos(tm_ref_work_t *w, uint32_t pos)
{
	for(size_t i = 0; i < w->size; i++) {
		*w->q++ = (tm_ref_kpos_t){
			.kmer = w->kmer[i].f,
			.pos = pos
		};
		*w->q++ = (tm_ref_kpos_t){
			.kmer = w->kmer[i].r,
			.pos = pos			/* FIXME: direction flag needed? */
		};
	}
	return;
}

static _force_inline
void tm_ref_push_base(tm_ref_work_t *w, uint8_t base)
{
	size_t const size = w->size;		/* save current size before duplicating stack */
	uint8_t const b2 = tm_ref_base_to_2bit(base);

	if(_unlikely(tm_ref_is_ambiguous(base))) {
		tm_ref_dup_stack(w);
		tm_ref_update_kmer(&w->kmer[size], size, w->amask, w->shift, b2>>2);
	}

	/* push base */
	tm_ref_update_kmer(&w->kmer[size], size, w->amask, w->shift, b2 & 0x03);
	return;
}

static _force_inline
tm_ref_kpos_t *tm_ref_reserve(tm_ref_kpos_v *buf, tm_ref_kpos_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	tm_ref_kpos_t *p = kv_reserve(tm_ref_kpos_t, *buf, 1024);
	return(&p[kv_cnt(*buf)]);
}

static _force_inline
size_t tm_ref_enumerate_kmer(size_t kbits, uint8_t const *seq, size_t slen, tm_ref_kpos_v *buf)
{
	tm_ref_work_t w = {
		.amask = (0x01<<kbits) - 0x01,
		.shift = kbits - 2,

		/* make room in dst array */
		.q = tm_ref_reserve(buf, kv_ptr(*buf)),

		/* clear stack; keeps all the combination of ambiguous bases */
		.size = 1,
		.branches = 0,
		.kmer = { { 0 } }
	};

	/* push bases at the head */
	size_t const k = kbits>>1;
	for(size_t i = 0; i < k - 1; i++) {
		tm_ref_push_base(&w, seq[i]);
	}

	/* body */
	for(size_t i = k - 1; i < slen; i++) {
		tm_ref_push_base(&w, seq[i]);			/* stack will be expanded if needed */
		tm_ref_push_pos(&w, i);

		if(_unlikely(w.branches & 0x01)) {
			tm_ref_shrink_stack(&w);			/* shrink if needed */
		}

		/*
		 * check need for expansion every 32times
		 * (Skylake's branch predictor successfully learns once-every-32times event, but fails for more)
		 */
		if((i & 0x1f) != 0) { continue; }
		w.q = tm_ref_reserve(buf, w.q);
	}

	/* done; record count */
	kv_cnt(*buf) = w.q - kv_ptr(*buf);
	return(kv_cnt(*buf));
}


/* rotate */
static _force_inline
size_t tm_ref_rotate_kmer(size_t kbits, tm_ref_kpos_t *kpos, size_t kcnt)
{
	uint32_t const mask = (0x01<<kbits) - 0x01;

	/* automatically vectorized? */
	for(size_t i = 0; i < kcnt; i++) {
		uint32_t const kmer = kpos[i].kmer;
		kpos[i].kmer = ((kmer & mask)<<2) | (kmer>>kbits);
	}
	return(kcnt);
}

static _force_inline
size_t tm_ref_collect_kmer(tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen, tm_ref_kpos_v *buf)
{
	/* slice k-mers (with additional succeeding base) */
	if(tm_ref_enumerate_kmer(conf->kbits + 2, seq, slen, buf) == 0) {
		return(0);				/* abort if no kmer found */
	}

	/* rotate by a base before sort */
	tm_ref_rotate_kmer(conf->kbits, kv_ptr(*buf), kv_cnt(*buf));

	/* sort kmers by position */
	radix_sort_kmer(kv_ptr(*buf), kv_cnt(*buf));
	return(kv_cnt(*buf));
}

/* pack kmer-position array */
typedef struct {
	uint32_t tail : 12;			/* pos count for this bin */
	uint32_t plen : 4;			/* unmatching prefix length */
	int32_t next  : 16;			/* diff to next bin */
} tm_ref_link_t;
_static_assert(sizeof(tm_ref_link_t) == sizeof(uint32_t));

typedef struct {
	tm_ref_link_t link[4];		/* de-Bruijn-graph-like link for each nucleotide */
	uint16_t pos[];
} tm_ref_bin_t;
_static_assert(sizeof(tm_ref_bin_t) == 16);		/* smallest bin size == 16 */

#define MK_REF_ALIGN_SIZE		( sizeof(uint64_t) )
#define MK_REF_BASE_OFS			( sizeof(tm_ref_sketch_t) )


typedef struct {
	uint32_t cnt, ofs;
} tm_ref_cnt_t;
typedef struct { tm_ref_cnt_t *a; size_t n, m; } tm_ref_cnt_v;

static _force_inline
size_t tm_ref_count_kmer(tm_ref_cnt_t *p, size_t ksize, tm_ref_kpos_t const *kpos, size_t kcnt)
{
	_unused(ksize);

	for(tm_ref_kpos_t const *k = kpos, *t = &kpos[kcnt]; k < t; k++) {
		p[k->kmer>>2].cnt++;
	}
	return(kcnt);
}

static _force_inline
size_t tm_ref_calc_size(tm_ref_cnt_t const *p, size_t ksize)
{
	size_t const hdr = sizeof(tm_ref_bin_t);

	size_t acc = 0;
	for(size_t i = 0; i < ksize; i++) {
		if(p[i].cnt == 0) { continue; }
		acc += _roundup(hdr + sizeof(uint16_t) * p[i].cnt, MK_REF_ALIGN_SIZE);
	}
	return(acc);
}

static _force_inline
size_t tm_ref_pack_kpos(tm_ref_cnt_t *p, size_t ksize, tm_ref_kpos_t const *kpos, size_t kcnt, tm_ref_sketch_t *sk)
{
	_unused(ksize);

	/* first bin */
	size_t const hdr = sizeof(tm_ref_bin_t);
	size_t ofs = MK_REF_BASE_OFS;

	tm_ref_kpos_t const *k = kpos, *t = &kpos[kcnt];
	while(k < t) {
		/* save current offset */
		uint32_t const kbase = k->kmer>>2;
		p[kbase].ofs = ofs;

		/* slice bin from current offset */
		tm_ref_bin_t *bin = _add_offset(sk, ofs);
		_storeu_v4i32(bin, _zero_v4i32());		/* clear link */

		/* pack pos array */
		size_t i = 0, j = 0;
		do {
			uint32_t const kall = k->kmer;
			do {
				bin->pos[i++] = k++->pos;
			} while(k < t && k->kmer == kall);

			do {
				bin->link[j++].tail = i;
			} while(j < (kall & 0x03));
		} while(k < t && (k->kmer>>2) == kbase);

		/* update offset */
		ofs += _roundup(hdr + sizeof(uint16_t) * p[kbase].cnt, MK_REF_ALIGN_SIZE);
	}
	return(ofs);
}

typedef struct {
	uint32_t len, ofs;
} tm_ref_prefix_t;

static _force_inline
tm_ref_prefix_t tm_ref_find_link(tm_ref_cnt_t const *p, uint32_t kmer, size_t max_depth)
{
	/* breadth-first search on kmer branching tree over profile table */
	size_t depth = 0;
	do {
		size_t const slen = 2 * depth;		/* substitution length (zero for the first iteration) */
		uint32_t const krem = kmer & (0 - (0x01<<slen));

		/* for each substituted kmer */
		for(size_t i = 0; i < 0x01ULL<<slen; i++) {
			uint32_t const ksub = krem + i;
			if(p[ksub].cnt == 0) { continue; }

			/* hit */
			return((tm_ref_prefix_t){
				.len = depth,
				.ofs = p[ksub].ofs
			});
		}
	} while(++depth < max_depth);

	/* never reach here?? */
	return((tm_ref_prefix_t){
		.len = max_depth,
		.ofs = MK_REF_BASE_OFS
	});
}

static _force_inline
void tm_ref_build_link(tm_ref_cnt_t const *p, size_t ksize, tm_ref_sketch_t *sk)
{
	size_t const shift = _tzcnt_u64(ksize);
	for(size_t i = 0; i < ksize; i++) {
		if(p[i].cnt == 0) { continue; }

		tm_ref_bin_t *bin = _add_offset(sk, p[i].ofs);

		/* find next link for each base */
		for(size_t j = 0; j < 4; j++) {
			tm_ref_prefix_t n = tm_ref_find_link(p, (i + (j<<shift))>>2, shift>>1);

			/* save link; preserve tail */
			bin->link[j].plen = n.len;
			bin->link[j].next = n.ofs - p[i].ofs;	/* overflows; keep sign bit */
		}
	}
	return;
}

static _force_inline
tm_ref_sketch_t *tm_ref_build_index(tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen, tm_ref_kpos_t const *kpos, size_t kcnt, tm_ref_cnt_v *buf)
{
	_unused(seq);
	kv_clear(*buf);

	/* reserve profile buffer */
	size_t const ksize = 0x01ULL<<conf->kbits;
	tm_ref_cnt_t *p = kv_reserve(tm_ref_cnt_t, *buf, ksize);

	/* count kmers */
	tm_ref_count_kmer(p, ksize, kpos, kcnt);

	/* accumulate block size */
	size_t size = tm_ref_calc_size(p, ksize);
	if(size > UINT32_MAX) { trap(); }

	/* malloc; use entire page */
	void *base = malloc(_roundup(size + conf->margin.head + conf->margin.tail, 4096));
	if(base == NULL) { return(NULL); }

	/* save metadata */
	tm_ref_sketch_t *sk = _add_offset(base, conf->margin.head);
	*sk = (tm_ref_sketch_t){
		.size = size,
		.head_margin = conf->margin.head,
		.kbits = conf->kbits,
		.slen = slen
	};

	/* pack pos array */
	tm_ref_pack_kpos(p, ksize, kpos, kcnt, sk);
	tm_ref_build_link(p, ksize, sk);

	/* anything else? */
	return(sk);
}


/* working buffer */
typedef struct {
	tm_ref_kpos_v kpos;
	tm_ref_cnt_v prof;
} tm_ref_tbuf_t;

static _force_inline
tm_ref_sketch_t *tm_ref_sketch(tm_ref_tbuf_t *self, tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen)
{
	/* slice kmers from sequence (ambiguous bases are expanded here) */
	if(tm_ref_collect_kmer(conf, seq, slen, &self->kpos) == 0) {		/* with additional base */
		return(NULL);
	}

	/* build index */
	tm_ref_kpos_t *kpos = kv_ptr(self->kpos);
	size_t const kcnt = kv_cnt(self->kpos);
	return(tm_ref_build_index(conf, seq, slen, kpos, kcnt, &self->prof));
}


/* matcher */
typedef struct {
	size_t src, dst, cnt;
} tm_ref_squash_t;

typedef struct {
	size_t prefix, unmatching;
	tm_ref_bin_t const *bin;
} tm_ref_state_t;

typedef struct {
	tm_ref_state_t state;
	tm_ref_squash_t squash;
} tm_ref_next_t;

typedef struct {
	uint16_t const *ptr;
	size_t cnt;
} tm_ref_match_t;


static _force_inline
tm_ref_state_t tm_ref_match_init(tm_ref_sketch_t const *ref)
{
	return((tm_ref_state_t){
		.prefix     = ref->kbits>>1,
		.unmatching = 1,
		.bin        = _add_offset(ref, MK_REF_BASE_OFS)	/* bin for AAA... */
	});
}

static _force_inline
tm_ref_next_t tm_ref_match_next(tm_ref_state_t s, uint8_t next)
{
	int64_t const pitch = MK_REF_ALIGN_SIZE;

	/* get link to the next bin */
	tm_ref_link_t const link = s.bin->link[next];
	tm_ref_bin_t const *bin = _add_offset(s.bin, pitch * link.next);

	/* update matching state */
	size_t const prefix = s.prefix - s.unmatching + link.plen;
	size_t const sidx = bin->link[next - 1ULL].tail;
	return((tm_ref_next_t){
		.state = {
			.prefix     = prefix,
			.unmatching = prefix != 0,
			.bin        = bin
		},

		/* extract squashable subbin */
		.squash = {
			.src = next == 0 ? 0 : sidx,
			.dst = bin->link[next].tail,
			.cnt = bin->link[3].tail
		}
	});
}

static _force_inline
tm_ref_match_t tm_ref_get_arr(tm_ref_state_t s)
{
	return((tm_ref_match_t){
		.ptr = s.bin->pos,
		.cnt = s.bin->link[3].tail
	});
}


/* index data structure */
typedef struct {
	/* fallback parameters */
	size_t kmer, window;		/* k-mer length and chain window size */

	/* extension params */
	uint64_t match, mismatch, gap_open, gap_extend;
	int64_t min_score;
} tm_idx_conf_t;

/* matcher object; k-mer length, score matrix, and score thresholds */
typedef struct {
	uint32_t size;

	/* k-mer size */
	uint32_t kbits;

	/* metadata */
	struct {
		char *name;
	} meta;

	/* chaining */
	struct {
		union {
			struct { uint32_t u, v; } sep;
			uint64_t all;
		} window;
		uint32_t min_scnt;
		uint32_t reserved;
	} chain;

	/* filtering */
	struct {
		uint8_t score_matrix[16];		/* match-mismatch score matrix */
		struct {
			uint8_t cv[16], pv[16], gv[16];	/* p = 0, p = 1, gap */
		} init;

		/* extension test if < test_cnt, do inward extension if longer than uspan_thresh */
		uint32_t test_cnt;
		int32_t min_score;		/* discard if sum of extension scores is smaller than min_score */
	} filter;

	/* extension */
	struct {
		dz_profile_t *dz;
		int32_t min_score;
		uint8_t giv, gev, gih, geh;
		int8_t score_matrix[DZ_MAT_SIZE * DZ_MAT_SIZE];
	} extend;
} tm_idx_profile_t;
// _static_assert(sizeof(tm_idx_profile_t) == 16);


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
#define tm_idx_ref_seq_ptr(x)		( (uint8_t const *)_add_offset((x), sizeof(tm_idx_sketch_hdr_t) + (x)->h.sofs) )
#define tm_idx_ref_seq_len(x)		( (x)->ref.slen )
#define tm_idx_ref_name_ptr(x)		( (uint8_t const *)_add_offset((x), sizeof(tm_idx_sketch_hdr_t) + (x)->ref.size) )
#define tm_idx_ref_name_len(x)		( strlen(tm_idx_ref_name_ptr(x)) )


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


static _force_inline
void tm_idx_destroy(tm_idx_t *mi)
{
	if(mi == NULL) { return; }
	if(mi->base != NULL) {
		free(mi->base);
		return;
	}

	for(size_t i = 0; i < mi->profile.cnt; i++) {
		free(mi->profile.arr[i]);
	}
	free(mi->profile.arr);

	for(size_t i = 0; i < mi->sketch.cnt; i++) {
		tm_ref_destroy(&mi->sketch.arr[i]->ref);
	}
	free(mi->sketch.arr);

	free(mi);
	return;
}


/* index builder */
typedef struct {
	re_pattern_t *name;
	re_pattern_t *comment;
} tm_idx_matcher_t;

typedef struct {
	tm_idx_t mi;

	struct {
		kvec_t(tm_idx_profile_t *) profile;
		kvec_t(tm_idx_matcher_t) matcher;
	} score;

	struct {
		bseq_file_t *fp;
		kvec_t(tm_idx_sketch_t *) bin;
	} col;

	/* working buffers */
	tm_ref_tbuf_t ref[];
} tm_idx_gen_t;

typedef struct {
	size_t base_rid, _pad[3];
	bseq_batch_t bin;
} tm_idx_batch_t;


static _force_inline
void tm_idx_push_profile(tm_idx_gen_t *mii, tm_idx_matcher_t matcher, tm_idx_profile_t *profile)
{
	kv_push(tm_idx_matcher_t, mii->score.matcher, matcher);
	kv_push(tm_idx_profile_t *, mii->score.profile, profile);
	return;
}

static _force_inline
size_t tm_idx_find_profile(tm_idx_gen_t const *mii, char const *name, size_t nlen, char const *comment, size_t clen)
{
	_unused(nlen);
	_unused(clen);

	tm_idx_matcher_t const *matcher = kv_ptr(mii->score.matcher);
	size_t const mcnt = kv_cnt(mii->score.matcher);

	/* regex match; find first */
	for(size_t i = 0; i < mcnt; i++) {
		tm_idx_matcher_t const *m = &matcher[i];

		/* result buffers */
		re_result_t nres = { NULL, NULL }, cres = { NULL, NULL };

		/* do we need to save the results? */
		if(m->name != NULL && !re_match(m->name, name, &nres)) {
			return(i);
		}
		if(m->comment != NULL && !re_match(m->comment, comment, &cres)) {
			return(i);
		}
	}

	/* not found; return the last one */
	return(mii->mi.profile.cnt - 1);
}

static _force_inline
toml_table_t *tm_idx_dump_toml(char const *fn)
{
	FILE *fp = fopen(fn, "r");
	if(fp == NULL) {
		error("failed to open file: `%s'. check file path and permission.", fn);
		return(NULL);
	}

	/* parse table */
	size_t const elen = 256;
	char error[elen];
	toml_table_t *table = toml_parse_file(fp, error, elen);
	if(table == NULL) {
		error("error occurred on parsing `%s': %s", fn, error);
		return(NULL);
	}

	/* done */
	fclose(fp);
	return(table);
}

static _force_inline
tm_idx_matcher_t tm_idx_parse_matcher(char const *pname, toml_table_t const *table)
{
	_unused(pname);

	tm_idx_matcher_t matcher = {
		.name    = NULL,
		.comment = NULL
	};

	/* construct sequence name and comment matcher */
	char const *sname   = toml_raw_in(table, "name");
	char const *comment = toml_raw_in(table, "comment");

	if(sname   != NULL) { matcher.name = re_compile(sname); }
	if(comment != NULL) { matcher.comment = re_compile(comment); }
	return(matcher);
}

static _force_inline
size_t tm_idx_profile_size(size_t plen)
{
	size_t const size = sizeof(tm_idx_profile_t) + _roundup(plen + 1, 16);
	return(size);
}

static _force_inline
char *tm_idx_profile_name_ptr(tm_idx_profile_t *profile)
{
	size_t const ofs = sizeof(tm_idx_profile_t);
	return(_add_offset(profile, ofs));
}

typedef void (*tm_idx_score_foreach_t)(void *opaque, int8_t *p, size_t i, size_t j);

static _force_inline
void tm_idx_score_foreach(tm_idx_profile_t *profile, void *opaque, tm_idx_score_foreach_t fp)
{
	for(size_t i = 0; i < DZ_MAT_SIZE; i++) {
		for(size_t j = 0; j < DZ_MAT_SIZE; j++) {
			 fp(opaque, &profile->extend.score_matrix[i * DZ_MAT_SIZE + j], i, j);
		}
	}
	return;
}

static
void tm_idx_fill_score_callback(int64_t *s, int8_t *p, size_t i, size_t j)
{
	*p = s[(i & j) == 0];
	return;
}

static _force_inline
void tm_idx_fill_score(tm_idx_profile_t *profile, int64_t m, int64_t x)
{
	int64_t s[2] = { m, x };
	tm_idx_score_foreach(profile,
		(void *)s,
		(tm_idx_score_foreach_t)tm_idx_fill_score_callback
	);
	return;
}

static _force_inline
void tm_idx_fill_default(tm_idx_profile_t *profile)
{
	/* default params */
	profile->kbits = 2 * 5;

	profile->chain.window.sep.u = 32;
	profile->chain.window.sep.v = 32;
	profile->chain.min_scnt = 2;

	profile->extend.min_score = 16;
	profile->extend.giv = 5;
	profile->extend.gev = 1;
	profile->extend.gih = 5;
	profile->extend.geh = 1;
	tm_idx_fill_score(profile, 2, -3);
	return;
}

static _force_inline
void tm_idx_override_default(tm_idx_profile_t *profile, tm_idx_conf_t const *conf)
{
	if(conf->kmer > 0) {
		profile->kbits = 2 * conf->kmer;
	}

	/* chaining */
	if(conf->window > 0) {
		profile->chain.window.sep.u = conf->window;
		profile->chain.window.sep.v = conf->window;
	}

	/* score matrix */
	if(conf->match > 0 && conf->mismatch > 0) {
		int64_t const m = conf->match;
		int64_t const x = -((int64_t)conf->mismatch);

		tm_idx_fill_score(profile, m, x);
	}

	/* gap penalties */
	if(conf->gap_open > 0) {
		profile->extend.giv = conf->gap_open;
		profile->extend.gih = conf->gap_open;
	}
	if(conf->gap_extend > 0) {
		profile->extend.gev = conf->gap_extend;
		profile->extend.geh = conf->gap_extend;
	}

	/* postprocess */
	if(conf->min_score > 0) {
		profile->extend.min_score = conf->min_score;
	}
	return;
}

typedef struct {
	int64_t acc, min, max;
	size_t cnt;
} tm_idx_calc_acc_t;

static
void tm_idx_acc_filter_score(tm_idx_calc_acc_t *acc, int8_t *p, size_t i, size_t j)
{
	int64_t const s = *p;
	tm_idx_calc_acc_t *q = &acc[(i & j) == 0];
	q->acc += s;
	q->min = MIN2(q->min, s);
	q->max = MAX2(q->max, s);
	q->cnt++;
	return;
}

static _force_inline
void tm_idx_calc_filter_score(tm_idx_profile_t *profile)
{
	/* count and accumulate */
	tm_idx_calc_acc_t c[2] = {
		{ 0, INT32_MAX, INT32_MIN, 0 },
		{ 0, INT32_MAX, INT32_MIN, 0 }
	};
	tm_idx_score_foreach(profile,
		(void *)c,
		(tm_idx_score_foreach_t)tm_idx_acc_filter_score
	);

	/* take average */
	int64_t const mave = c[0].acc / c[0].cnt;
	int64_t const m = (c[0].min == mave
		? c[0].min
		: MAX2(1, mave * 3 / 4)
	);

	int64_t const xave = c[1].acc / c[1].cnt;
	int64_t const x = (c[1].min == mave
		? c[1].min
		: MIN2(-1, xave * 3 / 4)
	);

	for(size_t i = 0; i < 16; i++) {
		profile->filter.score_matrix[i] = (i == 0) ? x : m;
	}
	return;
}

static _force_inline
void tm_idx_calc_filter_ivec(tm_idx_profile_t *profile)
{
	/* FIXME; initial vectors and gap vector */
	_unused(profile);

	return;
}

static _force_inline
void tm_idx_calc_filter_thresh(tm_idx_profile_t *profile)
{
	int64_t min_score = profile->extend.min_score / 4;
	profile->filter.min_score = MAX2(0, min_score);
	return;
}

static _force_inline
void tm_idx_calc_filter_params(tm_idx_profile_t *profile)
{
	tm_idx_calc_filter_score(profile);
	tm_idx_calc_filter_ivec(profile);
	tm_idx_calc_filter_thresh(profile);
	return;
}

static _force_inline
void *tm_idx_malloc(void *unused, size_t size)
{
	_unused(unused);
	return(malloc(size));
}

static _force_inline
void tm_idx_finalize_profile(tm_idx_profile_t *profile)
{
	dz_allocator_t alloc = {
		.ctx = NULL,
		.fp = tm_idx_malloc
	};
	dz_score_conf_t conf = {
		.score_matrix = profile->extend.score_matrix,
		.ins_open   = profile->extend.giv,
		.ins_extend = profile->extend.gev,
		.del_open   = profile->extend.gih,
		.del_extend = profile->extend.geh,

		/* fixed */
		.max_gap_len = 16,
		.full_length_bonus = 0
	};
	profile->extend.dz = dz_init_profile(&alloc, &conf);
	return;
}

static _force_inline
tm_idx_profile_t *tm_idx_default_profile(tm_idx_conf_t const *conf)
{
	/* name */
	char const *pname = "default";
	size_t const plen = strlen(pname);

	/* allocate profile object */
	size_t const size = tm_idx_profile_size(plen);
	tm_idx_profile_t *profile = malloc(size);
	memset(profile, 0, size);

	/* save size and name */
	profile->size = size;
	profile->meta.name = tm_idx_profile_name_ptr(profile);
	memcpy(profile->meta.name, pname, plen + 1);

	/* load default params */
	tm_idx_fill_default(profile);
	tm_idx_override_default(profile, conf);

	/* derive filtering parameters from extension parameters */
	tm_idx_calc_filter_params(profile);

	/* instanciate dz */
	tm_idx_finalize_profile(profile);
	return(profile);
}

static _force_inline
tm_idx_profile_t *tm_idx_parse_profile(tm_idx_profile_t const *template, char const *pname, toml_table_t const *table)
{
	size_t const plen = strlen(pname);	

	/* FIXME: save profile name */
	_unused(template);
	_unused(table);

	/* instanciate dz */
	// tm_idx_finalize_profile(profile);
	return(NULL);
}

static _force_inline
uint64_t tm_idx_load_score_core(tm_idx_gen_t *mii, tm_idx_profile_t const *template, toml_table_t const *root)
{
	char const *key;
	for(size_t i = 0; (key = toml_key_in(root, i)) != NULL; i++) {
		toml_table_t const *table = toml_table_in(root, key);

		/* check if the key has a corresponding table */
		if(table == NULL) {
			error("missing tabular content for key `%s'. ignoring.", key);
			continue;
		}

		/* parse name / comment matcher */
		tm_idx_matcher_t matcher = tm_idx_parse_matcher(key, table);
		if(matcher.name == NULL && matcher.comment == NULL) { continue; }	/* error when both NULL */

		/* parse score matrix */
		tm_idx_profile_t *profile = tm_idx_parse_profile(template, key, table);	/* NULL if error */
		if(profile == NULL) { return(1); }	/* unignorable error */

		/* save */
		tm_idx_push_profile(mii, matcher, profile);
	}
	return(0);			/* no error */
}

static _force_inline
uint64_t tm_idx_load_score(tm_idx_gen_t *mii, tm_idx_profile_t const *template, char const *fn)
{
	/* open file as toml */
	toml_table_t *root = tm_idx_dump_toml(fn);

	/* parse */
	uint64_t state = tm_idx_load_score_core(mii, template, root);
	toml_free(root);
	return(state);
}

static _force_inline
uint64_t tm_idx_gen_profile(tm_idx_gen_t *mii, tm_idx_conf_t const *conf, char const *fn, FILE *log)
{
	/* compose default (fallback) profile as template */
	tm_idx_profile_t *fallback = tm_idx_default_profile(conf);

	if(fn != NULL) {
		message(log, "reading score matrices...");
		if(tm_idx_load_score(mii, fallback, fn)) {
			return(1);
		}
	} else {
		message(log, "loading default score matrix...");
	}

	/* append default score matrix as fallback */
	tm_idx_matcher_t m = { NULL, NULL };
	tm_idx_push_profile(mii, m, fallback);
	return(0);
}


static
tm_idx_batch_t *tm_idx_fetch(uint32_t tid, tm_idx_gen_t *self)
{
	_unused(tid);		/* always run on the main thread */

	if(bseq_is_eof(self->col.fp)) { return(NULL); }	/* reached tail */

	bseq_batch_t *bin = bseq_read(self->col.fp);
	if(bin == NULL) { return(NULL); }		/* reached tail (or something is wrong?) */

	/* assign block id */
	tm_idx_batch_t *batch = _sub_offset(bin, offsetof(tm_idx_batch_t, bin));

	/* update sequence count */
	size_t const scnt = kv_cnt(self->col.bin);
	batch->base_rid = scnt;

	/* expand array */
	size_t const nscnt = scnt + bseq_meta_cnt(&batch->bin);
	kv_reserve(tm_idx_sketch_t *, self->col.bin, nscnt);
	kv_cnt(self->col.bin) = nscnt;

	return(batch);
}

static _force_inline
size_t tm_idx_roundup(size_t len)
{
	/* keep separation larger than 12 bases */
	return(_roundup(len + 12, 16));
}

static _force_inline
tm_ref_sketch_t *tm_idx_build_sketch(tm_idx_gen_t *self, tm_ref_tbuf_t *ref, size_t pid, char const *name, size_t nlen, uint8_t const *seq, size_t slen)
{
	_unused(name);

	/* calc margin size for saving name and seq */
	size_t const margin = tm_idx_roundup(nlen) + tm_idx_roundup(slen);

	/* retrieve profile */
	tm_idx_profile_t const *profile = self->mi.profile.arr[pid];

	/* pack args */
	tm_ref_conf_t conf = {
		.kbits = profile->kbits,
		.margin = {
			.head = sizeof(tm_idx_sketch_hdr_t),
			.tail = margin + 32			/* add 32 for vectorized access */
		}
	};
	return(tm_ref_sketch(ref, &conf, seq, slen));
}

static _force_inline
tm_idx_sketch_t *tm_idx_save_seq(tm_ref_sketch_t *sr, uint32_t rid, uint32_t pid, char const *name, size_t nlen, uint8_t const *seq, size_t slen)
{
	v32i8_t const z = _zero_v32i8();

	/* copy name */
	uint8_t *n = _add_offset(sr, sr->size);
	for(size_t i = 0; i < nlen; i++) {
		v16i8_t const v = _loadu_v16i8(&name[i]);
		_storeu_v16i8(&n[i], v);
	}
	_storeu_v32i8(&n[nlen], z);

	/* copy sequence */
	uint8_t *s = _add_offset(n, tm_idx_roundup(nlen));
	for(size_t i = 0; i < slen; i++) {
		v16i8_t const v = _loadu_v16i8(&seq[i]);
		_storeu_v16i8(&s[i], v);
	}
	_storeu_v32i8(&s[slen], z);

	/* build header */
	tm_idx_sketch_t *si = _sub_offset(sr, sizeof(tm_idx_sketch_hdr_t));
	size_t const size = (
		  sizeof(tm_idx_sketch_hdr_t)	/* header */
		+ sr->size						/* body */
		+ tm_idx_roundup(nlen)			/* name (margined) */
		+ tm_idx_roundup(slen)			/* sequence (margined) */
	);
	si->h = (tm_idx_sketch_hdr_t){
		.size = size,
		.rid = rid,
		.pid = pid,
		.sofs = size - tm_idx_roundup(slen)
	};
	return(si);
}

static
tm_idx_batch_t *tm_idx_collect(uint32_t tid, tm_idx_gen_t *self, tm_idx_batch_t *batch)
{
	/* load pointers */
	tm_ref_tbuf_t *ref = &self->ref[tid];
	bseq_meta_t *meta = bseq_meta_ptr(&batch->bin);
	uint32_t const bid = batch->bin.base_id;

	for(size_t i = 0; i < bseq_meta_cnt(&batch->bin); i++) {
		bseq_meta_t *p = &meta[i];

		/* get score matrix for this sequence from its name with regex matching */
		size_t pid = tm_idx_find_profile(self,
			bseq_name(p),    bseq_name_len(p),
			bseq_comment(p), bseq_comment_len(p)
		);

		/* build hash table */
		tm_ref_sketch_t *sr = tm_idx_build_sketch(self,
			ref, pid,
			bseq_name(p), bseq_name_len(p),
			bseq_seq(p),  bseq_seq_len(p)
		);

		/* copy sequence */
		tm_idx_sketch_t *si = tm_idx_save_seq(sr,
			bid + i, pid,	/* save reference id and profile id */
			bseq_name(p), bseq_name_len(p),
			bseq_seq(p),  bseq_seq_len(p)			
		);

		/* save ptr */
		p->u.ptr = si;
	}
	return(batch);
}

static
void tm_idx_record(uint32_t tid, tm_idx_gen_t *self, tm_idx_batch_t *batch)
{
	_unused(tid);

	/* we expect tm_idx_fetch and tm_idx_record are not run concurrently */
	tm_idx_sketch_t **arr = kv_ptr(self->col.bin);
	bseq_meta_t const *meta = bseq_meta_ptr(&batch->bin);

	/* no need for sorting batches */
	size_t const base_rid = batch->base_rid;
	for(size_t i = 0; i < bseq_meta_cnt(&batch->bin); i++) {
		arr[base_rid + i] = meta[i].u.ptr;			/* just copy */
	}

	/* done */
	bseq_free(&batch->bin);
	return;
}

static _force_inline
uint64_t tm_idx_gen_core(tm_idx_gen_t *mii, tm_idx_conf_t const *conf, char const *fn, pt_t *pt)
{
	_unused(conf);

	/* do not keep comment nor qual, do not use head margin */
	bseq_conf_t const bseq_conf = {
		.batch_size  = BATCH_SIZE,
		.head_margin = sizeof(tm_idx_batch_t)
	};
	mii->col.fp = bseq_open(&bseq_conf, fn);
	if(mii->col.fp == NULL) {
		error("failed to open file `%s', please check path and permission.", fn);
		return(1);
	}

	/* init dst array */
	kv_init(mii->col.bin);

	/* stream */
	pt_stream(pt, mii,
		(pt_source_t)tm_idx_fetch,
		(pt_worker_t)tm_idx_collect,
		(pt_drain_t)tm_idx_record
	);

	/* check error */
	bseq_close_t c = bseq_close(mii->col.fp);
	if(c.status != 0) {
		warn("broken file format detected for `%s'.", fn);
	}

	/* sanity check */
	size_t read = mii->mi.sketch.cnt;
	return(read == 0 || read != c.cnt || read != kv_cnt(mii->col.bin));
}

static _force_inline
void tm_idx_destroy_matcher(tm_idx_gen_t *mii)
{
	tm_idx_matcher_t const *matcher = kv_ptr(mii->score.matcher);
	size_t const mcnt = kv_cnt(mii->score.matcher);

	for(size_t i = 0; i < mcnt; i++) {
		tm_idx_matcher_t const *m = &matcher[i];

		if(m->name != NULL) { re_free(m->name); }
		if(m->comment != NULL) { re_free(m->comment); }
	}
	kv_destroy(mii->score.matcher);
	return;
}

static _force_inline
tm_idx_t *tm_idx_gen_finalize(tm_idx_gen_t *mii, char const *ref, char const *profile)
{
	/* cleanup regex matcher */
	tm_idx_destroy_matcher(mii);

	/* save filename */
	tm_idx_t *mi = &mii->mi;
	mi->filename.ref     = mm_strdup(ref);
	mi->filename.profile = mm_strdup(profile);

	/* copy pointers */
	mi->profile.arr = kv_ptr(mii->score.profile);
	mi->profile.cnt = kv_cnt(mii->score.profile);

	mi->sketch.arr = kv_ptr(mii->col.bin);
	mi->sketch.cnt = kv_cnt(mii->col.bin);

	/* clear working buffers */
	memset(_add_offset(mi, sizeof(tm_idx_t)), 0, sizeof(tm_idx_gen_t) - sizeof(tm_idx_t));
	return(mi);
}

static _force_inline
tm_idx_t *tm_idx_gen(tm_idx_conf_t const *conf, pt_t *pt, char const *ref, char const *profile, FILE *log)
{
	/* allocate thread-local buffers */
	size_t size = sizeof(tm_idx_gen_t) + pt_nth(pt) * sizeof(tm_ref_tbuf_t);
	tm_idx_gen_t *mii = aligned_malloc(size);
	memset(mii, 0, size);

	/* load score matrices and build regex matcher */
	if(tm_idx_gen_profile(mii, conf, profile, log)) {
		goto _tm_idx_gen_error;
	}

	/* read sequence and collect k-mers */
	message(log, "reading sequences and collecting k-mers...");
	if(tm_idx_gen_core(mii, conf, ref, pt)) {
		goto _tm_idx_gen_error;
	}

	message(log, "done.");
	return(tm_idx_gen_finalize(mii, ref, profile));

_tm_idx_gen_error:;
	kv_destroy(mii->col.bin);
	tm_idx_destroy(&mii->mi);
	return(NULL);
}


/* index I/O */
#define TM_IDX_MAGIC				( 0x30494d54 )

typedef struct {
	uint64_t magic;
	size_t size;
} tm_idx_hdr_t;

typedef struct {
	void *fp;
	write_t wfp;
	size_t ofs;
} tm_idx_dump_t;

/* dump utilities; offset accumulator */
static _force_inline
void *tm_idx_acc_ofs(size_t *p, size_t size)
{
	size_t ofs = *p;
	*p += size;
	return((void *)ofs);
}

static _force_inline
size_t tm_idx_padding_size(size_t ofs)
{
	return((64ULL - ofs) & (64ULL - 1));
}

/* dump block and return offset from dumped object as void * */
static _force_inline
void *tm_idx_dump_block(tm_idx_dump_t *w, void const *p, size_t size)
{
	w->wfp(w->fp, p, size);
	return(tm_idx_acc_ofs(&w->ofs, size));
}

static _force_inline
void *tm_idx_dump_pad(tm_idx_dump_t *w)
{
	size_t size = tm_idx_padding_size(w->ofs);
	if(size == 0) { return((void *)((uintptr_t)w->ofs)); }

	static uint8_t const pad[64] = { 0 };
	return(tm_idx_dump_block(w, pad, size));
}


static _force_inline
size_t tm_idx_scan_profile(tm_idx_t const *mi)
{
	size_t size = 0;

	for(size_t i = 0; i < mi->profile.cnt; i++) {
		size += mi->profile.arr[i]->size;		/* no padding between profiles */
	}
	size += tm_idx_padding_size(size);

	/* pointer table */
	size += sizeof(tm_idx_profile_t *) * mi->profile.cnt;
	size += tm_idx_padding_size(size);
	return(size);
}

static _force_inline
void *tm_idx_dump_profile_core(tm_idx_dump_t *w, tm_idx_profile_t *profile)
{
	uint8_t buf[profile->size + 1];
	tm_idx_profile_t *p = (tm_idx_profile_t *)buf;
	memcpy(p, profile, profile->size);

	/* ignore dz_profile_t object */
	p->meta.name = NULL;
	p->extend.dz = NULL;
	return(tm_idx_dump_block(w, p, p->size));
}

static _force_inline
void *tm_idx_dump_profile(tm_idx_dump_t *w, tm_idx_t const *mi)
{
	/* copy pointer array */
	tm_idx_profile_t **buf = (tm_idx_profile_t **)malloc(sizeof(tm_idx_profile_t *) * mi->profile.cnt);
	memcpy(buf, mi->profile.arr, sizeof(tm_idx_profile_t *) * mi->profile.cnt);

	/* dump body */
	for(size_t i = 0; i < mi->profile.cnt; i++) {
		buf[i] = tm_idx_dump_profile_core(w, buf[i]);
		/* we don't add padding between profiles */
	}
	tm_idx_dump_pad(w);

	/* dump offsets (pointers) */
	tm_idx_dump_block(w, buf, sizeof(tm_idx_profile_t *) * mi->profile.cnt);
	tm_idx_dump_pad(w);

	/* done */
	free(buf);
	return((void *)w->ofs);
}

static _force_inline
size_t tm_idx_scan_sketch(tm_idx_t const *mi)
{
	size_t size = 0;

	for(size_t i = 0; i < mi->sketch.cnt; i++) {
		size += mi->sketch.arr[i]->h.size;
		size += tm_idx_padding_size(size);
	}

	/* pointer table */
	size += sizeof(tm_idx_sketch_t *) * mi->sketch.cnt;
	size += tm_idx_padding_size(size);
	return(size);
}

static _force_inline
void *tm_idx_dump_sketch(tm_idx_dump_t *w, tm_idx_t const *mi)
{
	/* copy pointer array */
	tm_idx_sketch_t **buf = (tm_idx_sketch_t **)malloc(sizeof(tm_idx_sketch_t *) * mi->sketch.cnt);
	memcpy(buf, mi->sketch.arr, sizeof(tm_idx_sketch_t *) * mi->sketch.cnt);

	/* dump body */
	for(size_t i = 0; i < mi->sketch.cnt; i++) {
		buf[i] = tm_idx_dump_block(w, buf[i], buf[i]->h.size);
		tm_idx_dump_pad(w);		/* make aligned */
	}

	/* dump offsets */
	tm_idx_dump_block(w, buf, sizeof(tm_idx_sketch_t *) * mi->sketch.cnt);
	tm_idx_dump_pad(w);

	/* done */
	free(buf);
	return((void *)w->ofs);
}

static _force_inline
size_t tm_idx_calc_size(tm_idx_t const *mi)
{
	size_t size = 0;

	size += sizeof(tm_idx_t);
	size += strlen(mi->filename.ref) + 1;
	size += strlen(mi->filename.profile) + 1;
	size += tm_idx_padding_size(size);		/* make aligned */

	size += tm_idx_scan_profile(mi);
	size += tm_idx_scan_sketch(mi);
	return(size);
}

static _force_inline
size_t tm_idx_dump(tm_idx_t const *mi, void *fp, write_t wfp)
{
	/* copy on stack */
	tm_idx_t cmi = *mi;

	/* first compute entire size */
	tm_idx_hdr_t hdr = {
		.magic = TM_IDX_MAGIC,
		.size  = sizeof(tm_idx_hdr_t) + tm_idx_calc_size(&cmi)
	};

	/* init dump context */
	tm_idx_dump_t w = {
		.fp = fp,
		.wfp = wfp,
		.ofs = 0
	};

	tm_idx_dump_block(&w, &hdr, sizeof(tm_idx_hdr_t));
	cmi.profile.arr = tm_idx_dump_profile(&w, &cmi);
	cmi.sketch.arr = tm_idx_dump_sketch(&w, &cmi);

	/* input seqence names */
	cmi.filename.ref     = tm_idx_dump_block(&w, cmi.filename.ref, strlen(cmi.filename.ref) + 1);
	cmi.filename.profile = tm_idx_dump_block(&w, cmi.filename.profile, strlen(cmi.filename.profile) + 1);

	/* done */
	tm_idx_dump_pad(&w);
	tm_idx_dump_block(&w, &cmi, sizeof(tm_idx_t));
	return(w.ofs);
}


/* laod */
static _force_inline
tm_idx_t *tm_idx_slice_root(uint8_t *base, size_t size)
{
	tm_idx_t *mi = _sub_offset(&base[size], sizeof(tm_idx_t));

	/* adjust bkt and seq pointers */
	mi->profile.arr = _add_offset(mi->profile.arr, base);
	mi->sketch.arr  = _add_offset(mi->sketch.arr, base);

	/* save base pointer (to be freed correctly) */
	mi->base = (void *)base;
	return(mi);
}

static _force_inline
void tm_idx_patch_profile(tm_idx_t *mi, uint8_t *base)
{
	tm_idx_profile_t **p = mi->profile.arr;
	size_t const pcnt = mi->profile.cnt;
	for(size_t i = 0; i < pcnt; i++) {
		p[i] = _add_offset(p[i], base);
		p[i]->meta.name = tm_idx_profile_name_ptr(p[i]);
		tm_idx_finalize_profile(p[i]);		/* instanciate dz */
	}
	return;
}

static _force_inline
void tm_idx_patch_sketch(tm_idx_t *mi, uint8_t *base)
{
	tm_idx_sketch_t **s = mi->sketch.arr;
	size_t const scnt = mi->sketch.cnt;
	for(size_t i = 0; i < scnt; i++) {
		s[i] = _add_offset(s[i], base);
	}
	return;
}

static _force_inline
tm_idx_t *tm_idx_load(void *fp, read_t rfp)
{
	tm_idx_hdr_t hdr;
	size_t hdr_size = rfp(fp, &hdr, sizeof(tm_idx_hdr_t));
	if(hdr_size != sizeof(tm_idx_hdr_t)) {
		if(hdr_size > 0) { error("failed to read header"); }
		return(NULL);
	}

	if(hdr.magic != TM_IDX_MAGIC) {
		uint64_t base = ((hdr.magic ^ TM_IDX_MAGIC) & 0xffffff) == 0;
		error("invalid header%s. please rebuild the index.", base ? " (possibly index format is old)" : "");
		return(NULL);
	}

	/* bulk load */
	size_t idx_size = _roundup(hdr.size, ARCH_HUGEPAGE_SIZE);
	uint8_t *base = aligned_malloc(idx_size);
	size_t read = rfp(fp, base, hdr.size);

	/* sanity check */
	if(read != hdr.size) {
		error("truncated index file. please rebuild the index (%zu, %zu).", read, hdr.size);
		free(base);
		return(NULL);
	}

	/* restore pointers */
	tm_idx_t *mi = tm_idx_slice_root(base, hdr.size);
	mi->filename.ref     = _add_offset(mi->filename.ref, base);
	mi->filename.profile = _add_offset(mi->filename.profile, base);
	tm_idx_patch_profile(mi, base);
	tm_idx_patch_sketch(mi, base);
	return(mi);
}


/* seed */
typedef struct {
	uint32_t u, v;
} tm_seed_t;
_static_assert(sizeof(tm_seed_t) == 8);
typedef struct { tm_seed_t *a; size_t n, m; } tm_seed_v;

#define tm_seed_upos(x)			( (x).u )
KRADIX_SORT_INIT(seed, tm_seed_t, tm_seed_upos, 4);


typedef struct {
	uint16_t src, dst, cnt;
} tm_sqiv_t;
_static_assert(sizeof(tm_sqiv_t) == 6);
typedef struct { tm_sqiv_t *a; size_t n, m; } tm_sqiv_v;

/* chain */
typedef struct {
	/* base seed position; converted to u-v coordinate */
	tm_seed_t pos;

	/* reference and query side spans */
	struct {
		uint16_t u, v;
	} span;

	/* attributes */
	union {
		struct {
			uint32_t rid : 24;		/* reference id */
			uint32_t weight : 8;	/* weight; calculated from seed count and sequence length */
		} stat;
		uint32_t scnt;				/* chained seed count */
	} attr;
} tm_chain_t;
_static_assert(sizeof(tm_chain_t) == 16);
typedef struct { tm_chain_t *a; size_t n, m; } tm_chain_v;

#define tm_chain_attr(x)		( (x).attr.scnt )
KRADIX_SORT_INIT(chain, tm_chain_t, tm_chain_attr, 4);


/* alignment result */
typedef struct {
	uint32_t r, q;
} tm_pair_t;

typedef struct {
	/* rbt */
	rbt_header_t h;
	uint32_t qmax;

	/* stats */
	uint32_t min_weight;
	int64_t score;
	dz_alignment_t *aln;

	/* positions */
	tm_pair_t pos, span;

	/* alignment path */
	struct {
		uint8_t const *ptr;
		size_t len;
	} path;
} tm_aln_t;
_static_assert(sizeof(tm_aln_t) == 64);
typedef struct { tm_aln_t *a; size_t n, m; } tm_aln_v;

#define tm_aln_rbt_header(a)		( &(a)->h )
#define tm_aln_rbt_cmp(a, b)		( (a)->pos.q < (b)->pos.q )
#define tm_aln_ivt_cmp_head(a, b)	( (a)->qmax                > (b)->pos.q )
#define tm_aln_ivt_cmp_tail(a, b)	( (a)->pos.q               < (b)->pos.q + (b)->span.q )
#define tm_aln_ivt_cmp_iter(a, b)	( (a)->pos.q + (a)->span.q > (b)->pos.q + (b)->span.q )

static _force_inline
uint64_t tm_aln_ivt_update(tm_aln_t *parent, tm_aln_t *child)
{
	if(child == NULL) { return(1); }

	uint32_t cepos = child->pos.q + child->span.q;
	uint32_t qmax = parent->qmax;

	if(cepos > qmax) {
		parent->qmax = cepos;
		return(1);
	}
	return(0);
}

RBT_INIT_IVT(aln, tm_aln_t, tm_aln_rbt_header,
	tm_aln_rbt_cmp,
	tm_aln_ivt_update
);
RBT_INIT_ITER(aln, tm_aln_t, tm_aln_rbt_header,
	tm_aln_ivt_cmp_head,
	tm_aln_ivt_cmp_tail,
	tm_aln_ivt_cmp_iter
);


/* working buffers */
typedef struct {
	struct {
		tm_seed_v arr;		/* for sorting seeds */
		tm_sqiv_v sqiv;
	} seed;

	struct {
		tm_chain_v arr;
	} chain;

	struct {
		uint32_t min_weight;	/* working variable */
		dz_arena_t *fill, *trace;

		/* result bin */
		tm_aln_v arr;
	} extend;
} tm_scan_t;


static _force_inline
void tm_scan_init_static(tm_scan_t *self, size_t cnt)
{
	memset(self, 0, sizeof(tm_scan_t) * cnt);

	size_t const size = DZ_MEM_INIT_SIZE - _roundup(sizeof(dz_arena_t), DZ_MEM_ALIGN_SIZE);
	for(size_t i = 0; i < cnt; i++) {
		self->extend.fill = dz_arena_init(size);
		self->extend.trace = NULL;
	}
	return;
}

static _force_inline
void tm_scan_destroy_static(tm_scan_t *self, size_t cnt)
{
	_unused(self);
	_unused(cnt);

	return;
}

static _force_inline
void tm_scan_clear(tm_scan_t *self)
{
	kv_clear(self->seed.arr);
	kv_clear(self->seed.sqiv);
	kv_clear(self->chain.arr);
	kv_clear(self->extend.arr);
	return;
}


static _force_inline
tm_seed_t *tm_seed_reserve(tm_seed_v *buf, tm_seed_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	tm_seed_t *p = kv_reserve(tm_seed_t, *buf, 65536);
	return(&p[kv_cnt(*buf)]);
}

static _force_inline
tm_sqiv_t *tm_sqiv_reserve(tm_sqiv_v *buf, tm_sqiv_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	tm_sqiv_t *p = kv_reserve(tm_sqiv_t, *buf, 64);
	return(&p[kv_cnt(*buf)]);
}

static _force_inline
size_t tm_save_sqiv(tm_ref_state_t s, tm_ref_squash_t sq, tm_sqiv_t *r)
{
	/* just copy */
	*r = (tm_sqiv_t){
		.src = sq.src,
		.dst = sq.dst,
		.cnt = sq.cnt
	};
	return(1 - s.unmatching);
}

static _force_inline
size_t tm_expand_seed(tm_ref_state_t s, v4i32_t uofs, v4i32_t vofs, tm_seed_t *q)
{
	v4i32_t const wmask = _set_v4i32(0x20000000);

	tm_ref_match_t m = tm_ref_get_arr(s);
	for(size_t i = 0; i < m.cnt; i += 4) {
		v4i32_t const z = ({
			__m128i const x = _mm_loadl_epi64((__m128i const *)&m.ptr[i]);	/* movq (mem), xmm */
			__m128i const y = _mm_cvtepi16_epi32(x);	/* pmovsxwd; sign expansion */
			((v4i32_t){ .v1 = y });
		});
		v4i32_t const w = _and_v4i32(z, wmask);			/* disjoin forward and reverse */

		v4i32_t const u = _sub_v4i32(uofs, w);
		v4i32_t const v = _add_v4i32(vofs, _add_v4i32(w, w));

		/* save */
		_storeu_v4i32(&q[i],     _lo_v4i32(u, v));
		_storeu_v4i32(&q[i + 2], _hi_v4i32(u, v));
	}
	return(m.cnt);
}

static _force_inline
size_t tm_collect_seed(tm_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, tm_seed_v *seed, tm_sqiv_v *sqiv)
{
	/* load coordinate constants */
	v4i32_t const uinc = _set_v4i32(2), vinc = _set_v4i32(-1);
	v4i32_t uofs = _set_v4i32(0x20000), vofs = _set_v4i32(0x10000);

	tm_seed_t *q = tm_seed_reserve(seed, kv_ptr(*seed));
	tm_sqiv_t *r = tm_sqiv_reserve(sqiv, kv_ptr(*sqiv));

	/* initial state (of three general-purpose registers) */
	tm_ref_state_t s = tm_ref_match_init(ref);

	/* for each base... */
	for(size_t i = 0; i < qlen; i++) {
		/* update coordinates */
		uofs = _add_v4i32(uofs, uinc);
		vofs = _add_v4i32(vofs, vinc);

		/* update matching status for the next bin (prefetch) */
		tm_ref_next_t n = tm_ref_match_next(s, query[i]);
		r += tm_save_sqiv(n.state, n.squash, r);

		/* skip current bin if not matching */
		if(!s.unmatching) {
			q += tm_expand_seed(s, uofs, vofs, q);
		}

		/* alias new state (we expect the next bin has already arrived) */
		s = n.state;

		/* test array expansion every 32 times */
		if((i & 0x1f) != 0) { continue; }
		q = tm_seed_reserve(seed, q);
		r = tm_sqiv_reserve(sqiv, r);
	}
	return(kv_cnt(*seed));
}

static _force_inline
size_t tm_squash_seed(tm_sqiv_t const *sqiv, size_t icnt, tm_seed_v *seed)
{
	tm_seed_t *q = kv_ptr(*seed);

	/* squash succeeding match */
	size_t idx = 0;
	for(size_t i = 0; i < icnt; i++) {
		tm_sqiv_t const *p = &sqiv[i];

		/* determine source and destination indices */
		size_t const src = idx + p->src;
		size_t const dst = idx + p->dst;
		size_t const tail = idx + p->cnt;

		/* copy with ymm */
		for(size_t i = src, j = dst; j < tail; i += 8, j += 8) {
			v32i8_t const v = _loadu_v32i8(&q[i]);
			_storeu_v32i8(&q[j], v);
		}

		/* update index for next iteration */
		idx += tail - dst + src;
	}

	kv_cnt(*seed) = idx;
	return(idx);
}

static _force_inline
int64_t tm_chain_test_ptr(tm_seed_t const *p, tm_seed_t const *t)
{
	return((int64_t)((ptrdiff_t)(t - p - 1)));
}

static _force_inline
tm_seed_t *tm_chain_find_first(tm_seed_t *p, tm_seed_t *t, uint64_t lb, uint64_t window)
{
	/* constants */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const tmask = 0x8000000080000000;	/* extract two sign bit pair */

	/* load bounds */
	uint64_t const uv = _loadu_u64(p);			/* (ulb, vlb) */
	uint64_t ub = (uv & umask) + window;		/* (uub, vub - vlb) */

	int64_t cont = 0;		/* continuous flag */
	while(1) {
		/* check if reached tail */
		if((cont | tm_chain_test_ptr(++p, t)) < 0) { return(NULL); }

		uint64_t const v = _loadu_u64(p) - lb;	/* 1: (u, v - vlb) for inclusion test */
		if(((ub - v) & tmask) == 0) { break; }	/* 2,3: break if chainable (first chainable seed found) */

		/* unchainable; test if out of uub */
		cont = ub - v;		/* save diff for testing MSb */
	}
	return(p);
}

static _force_inline
tm_seed_t *tm_chain_find_alt(tm_seed_t *p, tm_seed_t *t, uint64_t lb)
{
	/* constants */
	uint64_t const tmask = 0x8000000080000000;	/* extract two sign bit pair */

	/* calculate bounds */
	uint64_t rv = _loadu_u64(p) - lb;
	uint64_t ub = rv + (rv<<32);

	/* keep nearest seed */
	tm_seed_t *n = p;

	int64_t cont = 0;
	while((cont | tm_chain_test_ptr(++p, t)) >= 0) {
		uint64_t const v = _loadu_u64(p) - lb;	/* 1: (u, v - vlb) for inclusion test */
		cont = ub - v;
		if((ub - v) & tmask) { continue; }		/* skip if unchainable */

		/* chainable; test if the seed is nearer than the previous */
		uint64_t const w = v + (v<<32);			/* 2,3: (u + v - vlb, v - vlb) */
		if((ub - w) & tmask) { break; }			/* 4,5: further than previous */

		/* nearer seed found */
		ub = w;			/* update bounds */
		n = p;			/* save pointer */
	}
	return(n);
}

static _force_inline
tm_seed_t *tm_chain_find_nearest(tm_seed_t *p, tm_seed_t *t, uint64_t window)
{
	/* load root positions */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const uv = _loadu_u64(p);			/* (ulb, vlb) */
	uint64_t lb = (uv & ~umask);				/* (0, vlb) */

	/* find first chainable seed */
	tm_seed_t *n = tm_chain_find_first(p, t, lb, window);
	if(n == NULL) { return(NULL); }				/* we expect control paths are kept isolated */

	/* then inspect alternative chainable seeds */
	return(tm_chain_find_alt(n, t, lb));
}

static _force_inline
size_t tm_chain_record(v2i32_t root, v2i32_t tail, v2i32_t scnt, v2i32_t min_scnt, tm_chain_t *q)
{
	/* avoid use of general-purpose register to minimize spills */
	v2i32_t const span = _sub_v2i32(tail, root);

	/* FIXME: raw SSE operations */
	v4i32_t const chain = ({
		__m128i const x = _mm_shufflelo_epi16(span.v1, 0xd4);	/* (uint32_t [4]){ -, -, -, (vspan, uspan) } */
		__m128i const y = _mm_blend_epi16(x, scnt.v1, 0x0c);	/* (uint32_t [4]){ -, -, scnt, (vspan, uspan) } */
		__m128i const z = _mm_unpacklo_epi64(root.v1, y);		/* (uint32_t [4]){ scnt, (vspan, uspan), vpos, upos } */
		((v4i32_t){ .v1 = z });
	});
	_storeu_v4i32(q, chain);

	/* forward if scnt > min_scnt */
	v2i32_t const inc = _set_v2i32(1);		/* kept on register */
	v2i32_t const fwd = _and_v2i32(inc, _gt_v2i32(scnt, min_scnt));
	return(_ext_v2i32(fwd, 0));				/* movq */
}

static _force_inline
size_t tm_chain_seed(tm_idx_profile_t const *profile, tm_seed_t *seed, size_t scnt, tm_chain_v *chain)
{
	/*
	 * greedy extension; relatively smaller window size (~32) is preferred for good performance.
	 * do all of the chainablity test on general-purpose register for smaller latency of core loop.
	 * keep seed and counter on xmm regsters and do everything of the chain composing on xmm registers
	 * to minimize spill of general-purpose registers.
	 */
	uint64_t const chained = 0x40000000;				/* offset */
	uint64_t const window = profile->chain.window.all;	/* window sizes */

	v2i32_t const inc = _set_v2i32(1);
	v2i32_t const min_scnt = _load_v2i32(&profile->chain.min_scnt);

	/* src pointers */
	tm_seed_t *p = seed - 1, *t = &seed[scnt];

	/* reserve mem for chain */
	tm_chain_t *q = kv_reserve(tm_chain_t, *chain, kv_cnt(*chain) + scnt);

	/* for each seed */
	while(++p < t) {
		/* skip chained seed to find next root */
		if(p->u >= chained) { continue; }

		/* root found; keep seed on xmm register */
		v2i32_t const root = _loadu_v2i32(p);	/* we expect this won't be spilled */
		v2i32_t scnt = _zero_v2i32();

		/* iteratively link nearest seed */
		tm_seed_t *s = p;
		while(1) {
			tm_seed_t *n = tm_chain_find_nearest(s, t, window);
			if(n == NULL) { break; }

			/* increment seed count */
			scnt = _add_v2i32(scnt, inc);

			n->u += chained;	/* mark chained */
			s = n;				/* save last chained seed */
		}

		q += tm_chain_record(root, _loadu_v2i32(s), scnt, min_scnt, q);
	}

	/* update chain count */
	size_t const ccnt = q - &kv_ptr(*chain)[kv_cnt(*chain)];
	kv_cnt(*chain) = ccnt;
	return(ccnt);
}

static _force_inline
v16i8_t tm_load_forward(uint8_t const *ptr)
{
	return(_loadu_v16i8(ptr));
}

static _force_inline
v16i8_t tm_load_reverse(uint8_t const *ptr)
{
	static uint8_t const revcomp[16] __attribute__(( aligned(16) )) = {
		0
	};

	v16i8_t const rc = _load_v16i8(revcomp);
	v16i8_t const v = _loadu_v16i8(ptr);
	return(_shuf_v16i8(rc, v));
}

typedef struct {
	uint8_t const *ptr;
	v16i8_t (*load)(uint8_t const *);
	ptrdiff_t fwd;
} tm_filter_t;

static _force_inline
int64_t tm_filter_extend(tm_idx_profile_t const *profile, tm_filter_t r, tm_filter_t q)
{
	/* load sequence vectors */
	v16i8_t rv = r.load(r.ptr), rn = r.load(r.ptr + r.fwd);
	v16i8_t qv = q.load(q.ptr), qn = q.load(q.ptr + q.fwd);

	/* load constants */
	v16i8_t const score_matrix = _load_v16i8(&profile->filter.score_matrix);
	v16i8_t dh = _load_v16i8(&profile->filter.init);
	v16i8_t dv = _load_v16i8(&profile->filter.init);
	v16i8_t sv = _zero_v16i8(), mv = _zero_v16i8();

	#define _calc_a(_rv, _qv, _dh, _dv) ({ \
		/* calc score profile vector with shuffling score matrix vector */ \
		v16i8_t const match = _and_v16i8(_rv, _qv);		/* cmpeq */ \
		v16i8_t const a = _max_v16i8( \
			_shuf_v16i8(score_matrix, match), \
			_max_v16i8(_dv, _dh) \
		); \
		a; \
	})

	/* 32 vector updates */
	for(size_t j = 0; j < 16; j++) {
		/* #0 (even) */
		v16i8_t const a0 = _calc_a(rv, qv, dh, dv);
		v16i8_t const tdv = _sub_v16i8(a0, dh);
		v16i8_t const tdh = _sub_v16i8(a0, dv);
		sv = _add_v16i8(sv, tdv);
		mv = _max_v16i8(mv, sv);

		/* move rightward */
		rv = _bsrd_v16i8(rv, rn, 1); rn = _shr_v16i8(rn, 1);

		/* #1 (odd) */
		v16i8_t const a1 = _calc_a(rv, qv, tdh, tdv);
		dv = _sub_v16i8(a1, tdh);
		dh = _sub_v16i8(a1, tdv);
		sv = _add_v16i8(sv, dh);
		mv = _max_v16i8(mv, sv);

		/* move downward */
		qv = _bsrd_v16i8(qv, qn, 1); qn = _shr_v16i8(qn, 1);
	}

	/* fold scores */
	int64_t const score = 0;
	return(score);

	#undef _calc_a
}

static _force_inline
int64_t tm_filter_extend_inward(tm_idx_profile_t const *profile, uint8_t const *ref, uint8_t const *query, tm_chain_t const *chain)
{
	_unused(profile);
	_unused(ref);
	_unused(query);
	_unused(chain);
	return(0);
}

static _force_inline
int64_t tm_filter_extend_outward(tm_idx_profile_t const *profile, uint8_t const *ref, uint8_t const *query, tm_chain_t const *chain)
{
	_unused(profile);
	_unused(ref);
	_unused(query);
	_unused(chain);
	return(0);
}

static _force_inline
size_t tm_filter_save_chain(uint32_t rid, tm_chain_t *q, tm_chain_t const *p)
{
	/* load */
	v16i8_t v = _loadu_v16i8(p);

	/* calc weight */
	ZCNT_RESULT size_t zc = _lzcnt_u32(p->attr.scnt);
	uint32_t const weight = 32 - zc;

	/* save */
	_storeu_v16i8(q, v);
	q->attr.stat.rid = rid;
	q->attr.stat.weight = weight;
	return(1);
}

static _force_inline
size_t tm_filter_chain(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, tm_chain_t *chain, size_t ccnt)
{
	/* alias pointer */
	tm_chain_t const *src = chain;
	tm_chain_t *dst = chain;

	uint32_t const rid = tm_idx_ref_rid(si);
	uint8_t const *ref = tm_idx_ref_seq_ptr(si);

	/* load constants */
	uint64_t const test_cnt = profile->filter.test_cnt;
	int64_t const min_score = profile->filter.min_score;
	size_t const uspan_thresh =	32;		/* FIXME: profile->filter.uspan_thresh; */

	for(size_t i = 0; i < ccnt; i++) {
		tm_chain_t const *p = &src[i];

		/* try short extension if not heavy enough */
		if(p->attr.scnt < test_cnt) {
			int64_t score = (p->span.u > uspan_thresh
				? tm_filter_extend_inward(profile, ref, query, p)		/* inward extension for long-spanning chain */
				: tm_filter_extend_outward(profile, ref, query, p)
			);
			if(score < min_score) { continue; }
		}

		/* copy and fold in reference id */
		dst += tm_filter_save_chain(rid, dst, p);
	}
	return(ccnt);
}

static _force_inline
size_t tm_seed_and_sort(tm_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, tm_seed_v *seed, tm_sqiv_v *sqiv)
{
	/* enumerate seeds */
	kv_clear(*seed);
	kv_clear(*sqiv);
	if(tm_collect_seed(ref, query, qlen, seed, sqiv) == 0) {
		return(0);
	}

	/* squash overlapping seeds */
	if(tm_squash_seed(kv_ptr(*sqiv), kv_cnt(*sqiv), seed) == 0) {
		return(0);
	}

	/* sort seeds by u-coordinates for chaining */
	radix_sort_seed(kv_ptr(*seed), kv_cnt(*seed));
	return(kv_cnt(*seed));
}

static _force_inline
size_t tm_chain_and_filter(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, tm_seed_t *seed, size_t scnt, tm_chain_v *chain)
{
	/* chaining; return if no chain obtained */
	size_t const cbase = kv_cnt(*chain);		/* save chain count */
	size_t const ccnt  = tm_chain_seed(profile, seed, scnt, chain);
	if(ccnt == 0) { return(0); }

	/* filter (try small extension with simple score matrix) */
	size_t const cfilt = tm_filter_chain(si, profile, query, &kv_ptr(*chain)[cbase], ccnt);
	if(cfilt == 0) { return(0); }

	/* finalize filtered chains */
	kv_cnt(*chain) = cbase + cfilt; 
	return(kv_cnt(*chain));
}

static _force_inline
size_t tm_collect_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *seq, size_t slen)
{
	tm_idx_sketch_t const **si = (tm_idx_sketch_t const **)idx->sketch.arr;
	tm_idx_profile_t const **pf = (tm_idx_profile_t const **)idx->profile.arr;

	/* clear result bin */
	kv_clear(self->chain.arr);

	/* collect chain for each reference sequence */
	for(size_t i = 0; i < idx->sketch.cnt; i++) {
		/* collect seeds */
		if(tm_seed_and_sort(&si[i]->ref, seq, slen, &self->seed.arr, &self->seed.sqiv) == 0) { continue; }

		/* chain and filter */
		tm_seed_t *seed = kv_ptr(self->seed.arr);
		size_t scnt = kv_cnt(self->seed.arr);
		tm_chain_and_filter(si[i], pf[si[i]->h.pid], seq, seed, scnt, &self->chain.arr);
	}
	return(kv_cnt(self->chain.arr));
}



/* X-drop DP extension */

static _force_inline
tm_aln_t tm_chain_as_aln(tm_scan_t const *self, tm_chain_t const *chain)
{
	_unused(self);
	_unused(chain);

	/* FIXME */
	return((tm_aln_t){
		.pos  = { .r = 0, .q = 0 },
		.span = { .r = 0, .q = 0 },
		.qmax = 0
	});
}


static _force_inline
uint64_t tm_extend_is_complete(tm_scan_t *self)
{
	_unused(self);

	/* FIXME */
	return(0);
}

static _force_inline
uint64_t tm_extend_is_covered(tm_scan_t *self, tm_chain_t const *chain)
{
	/* convert chain position to (pseudo) alignment range */
	tm_aln_t const caln = tm_chain_as_aln(self, chain);
	uint32_t const weight = chain->attr.stat.weight;

	/* for each overlapping alignment */
	tm_aln_t const *v = kv_ptr(self->extend.arr);

	rbt_iter_t it;
	rbt_init_iter_aln(&it, v, &caln);

	tm_aln_t const *p = rbt_fetch_head_aln(&it, v, &caln);
	while(p != NULL) {
		if(p->min_weight > weight) {
			return(1);		/* already covered */
		}

		p = rbt_fetch_next_aln(&it, v, &caln);
	}
	return(0);				/* possiblilty remains for better alignment */
}

static _force_inline
void tm_extend_record(tm_scan_t *self, tm_chain_t const *chain, tm_aln_t aln)
{
	_unused(chain);

	kv_push(tm_aln_t, self->extend.arr, aln);
	rbt_insert_aln(kv_ptr(self->extend.arr), kv_cnt(self->extend.arr) - 1);
	return;
}

static _force_inline
void tm_extend_clear(tm_scan_t *self)
{
	kv_clear(self->extend.arr);
	dz_arena_flush(self->extend.trace);
	return;
}

typedef struct {
	uint8_t const *ptr;
	uint32_t len;
	uint32_t dir;
} tm_extend_ref_t;

static _force_inline
dz_state_t const *tm_extend_seed_core(dz_arena_t *mem, dz_profile_t const *profile, tm_extend_ref_t const *ref, uint8_t const *query, size_t qlen, size_t qspos)
{
	/* build query-side profile */
	dz_query_conv_t qconv = {
		.single = dz_query_conv_single,
		.bulk   = dz_query_conv_bulk
	};
	dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, (char const *)query, qlen);
	if(qspos & 0x80000000) {
		dz_pack_query_reverse(q, profile, qconv, (char const *)&query[qspos], qspos);
	} else {
		dz_pack_query_forward(q, profile, qconv, (char const *)&query[qspos], qlen - qspos);
	}

	/* reference-side */
	dz_ref_t r = {
		.ptr = (char const *)ref->ptr,
		.len = ref->len,
		.dir = ref->dir,

		.id  = 0,	/* rid is always zero */
		.init_s = INT16_MIN
	};

	dz_state_t const *ff = dz_extend_core(
		mem, profile, q, &r,
		dz_root(profile), 1		/* always use root */
	);
	return(ff);
}

static _force_inline
tm_aln_t tm_extend_seed(tm_scan_t *self, tm_idx_profile_t const *pf, tm_idx_sketch_t const *sk, size_t rrpos, uint8_t const *query, size_t qlen, size_t qrpos)
{
	dz_arena_flush(self->extend.fill);

	/* get reference sequence */
	uint8_t const *ref = tm_idx_ref_seq_ptr(sk);
	size_t const rlen = tm_idx_ref_seq_len(sk);

	/* forward */
	tm_extend_ref_t rf = {
		.ptr = &ref[rrpos],
		.len = rlen - rrpos,
		.dir = 0
	};
	dz_state_t const *r = tm_extend_seed_core(self->extend.fill, pf->extend.dz, &rf, query, qlen, qrpos);

	/* get downward max */
	uint64_t const spos = dz_calc_max_pos_core(r);
	uint64_t const rspos = spos>>32;
	uint64_t const qspos = spos & 0xffffffff;

	/* backward */
	tm_extend_ref_t rr = {
		.ptr = &ref[rspos],
		.len = rspos,
		.dir = 1
	};
	dz_state_t const *f = tm_extend_seed_core(self->extend.fill, pf->extend.dz, &rr, query, qlen, qspos);

	/* traceback */
	dz_alignment_t *aln = dz_trace_core(self->extend.trace, pf->extend.dz, f);

	uint32_t const rspan = aln->ref_length, qspan = aln->query_length;
	return((tm_aln_t){
		/* coordinates */
		.qmax = qspos + qspan,
		.pos  = { .r = rspos, .q = qspos },
		.span = { .r = rspan, .q = qspan },

		/* path */
		.path = { .ptr = aln->path, .len = aln->path_length },

		/* stats */
		.min_weight = 0,
		.score = aln->score,
		.aln = aln				/* save original */
	});
}

static _force_inline
size_t tm_extend_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *query, size_t qlen, tm_chain_t const *chain, size_t ccnt)
{
	tm_idx_sketch_t const **sk = (tm_idx_sketch_t const **)idx->sketch.arr;
	tm_idx_profile_t const **pf = (tm_idx_profile_t const **)idx->profile.arr;

	/* clear bin */
	tm_extend_clear(self);

	/* for each chain try extension */
	for(size_t i = 0; i < ccnt; i++) {
		tm_chain_t const *q = &chain[i];
		if(q->attr.stat.weight < self->extend.min_weight) { break; }

		/* skip if already covered (and there is no possibility that this chain surpasses the previous ones) */
		if(tm_extend_is_covered(self, q)) { continue; }

		tm_aln_t const caln = tm_chain_as_aln(self, q);
		tm_idx_sketch_t const *s = sk[q->attr.stat.rid];
		tm_idx_profile_t const *p = pf[s->h.pid];

		/* extend */
		tm_aln_t aln = tm_extend_seed(self,
			p, s, caln.pos.r,			/* reference */
			query, qlen, caln.pos.q		/* query */
		);
		if(aln.score <= p->extend.min_score) { continue; }

		tm_extend_record(self, q, aln);
		if(tm_extend_is_complete(self)) { break; }
	}

	/* anything else? */
	return(kv_cnt(self->extend.arr));
}


/* evaluate all; query sequence be shorter than 2Gbp */
static _force_inline
size_t tm_scan_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *seq, size_t slen)
{
	tm_scan_clear(self);

	/* seed and chain */
	if(tm_collect_all(self, idx, seq, slen) == 0) {
		return(0);
	}

	/* sort by estimated identity */
	tm_chain_t *cptr = kv_ptr(self->chain.arr);
	size_t const ccnt = kv_cnt(self->chain.arr);
	radix_sort_chain(cptr, ccnt);

	/* convert chain to alignment */
	size_t acnt = tm_extend_all(self, idx, seq, slen, cptr, ccnt);
	return(acnt);
}




/* printer */
typedef struct {
	uint32_t dummy;
} tm_print_conf_t;

typedef struct {
	uint32_t dummy;
} tm_print_t;


static _force_inline
void tm_print_destory_static(tm_print_t *self)
{
	_unused(self);

	return;
}

static _force_inline
void tm_print_init_static(tm_print_t *self, tm_print_conf_t const *conf, char const *args)
{
	_unused(self);
	_unused(conf);
	_unused(args);

	return;
}





/* multithreaded scan-and-mask */
typedef struct {
	uint32_t dummy;
} tm_mtscan_t;

static _force_inline
tm_mtscan_t *tm_mtscan_init(tm_idx_t const *mi, tm_print_t *printer, pt_t *pt)
{
	_unused(mi);
	_unused(printer);
	_unused(pt);

	return(NULL);
}

static _force_inline
void tm_mtscan_destroy(tm_mtscan_t *self)
{
	free(self);
	return;
}

static _force_inline
int tm_mtscan_file(tm_mtscan_t *self, char const *fn)
{
	_unused(self);
	_unused(fn);

	return(0);
}


/* return codes */
enum main_error_codes {
	ERROR_INTERNAL = 255,
	ERROR_NO_ARG = 1,

	ERROR_OPEN_IDX = 2,
	ERROR_LOAD_IDX = 3,

	ERROR_OPEN_RSEQ = 4,
	ERROR_OPEN_QSEQ = 5
};

typedef struct {
	/* index dump mode if not NULL */
	char const *idxdump;

	/* global args */
	char const *args;
	uint32_t verbose, help;
	size_t nth;

	/* scan-and-mask params */
	char const *profile;
	tm_idx_conf_t fallback;
	tm_print_conf_t print;

	/* option parser */
	FILE *log;
	opt_t opt;
} tm_conf_t;




#define tm_conf_append(_type, _ptr, _cnt, _body) { \
	kvec_t(_type) b; \
	kv_build(b, (void *)_ptr, _cnt, _cnt); \
	{ _body; } \
	_ptr = kv_ptr(b); \
	_cnt = kv_cnt(b); \
}


static void tm_conf_preset(opt_t *opt, tm_conf_t *conf, char const *arg)
{
	struct tm_conf_preset_s {
		/* key-value pair */
		char const *key;
		char const *val;

		/* child */
		struct tm_conf_preset_s const *children[6];
	};

	/* preset is defined recursively */
	#define _n(_k, ...)		&((struct tm_conf_preset_s const){ (_k), __VA_ARGS__ })
	struct tm_conf_preset_s const *presets[] = { 0 };		/* FIXME */
	#undef _n

	struct tm_conf_preset_s const *const *q = presets;
	split_foreach(arg, 0, ".:", {		/* traverse preset param tree along with parsing */
		while(*q != NULL && strncmp(p, (*q)->key, l) != 0) { q++; }
		if(*q == NULL) {				/* terminate if not matched, not loaded from file */
			int ret = opt_load_conf(opt, conf, p);
			oassert(opt, ret, "no preset params found for `%.*s'.", (int)l, p);
			break;
		}
		opt_parse_line(opt, conf, (*q)->val);/* apply recursively */
		q = (*q)->children;				/* visit child nodes */
	});
	return;
}

static void tm_conf_verbose(opt_t *opt, tm_conf_t *conf, char const *arg) {
	if(arg == NULL || *arg == '\0') { conf->verbose = 1; return; }
	if(isdigit(*arg)) { conf->verbose = (size_t)mm_atoi(arg, 0); return; }
	if(*arg != 'v') { return; }
	tm_conf_verbose(opt, conf, arg + 1);
	conf->verbose++;
}
static void tm_conf_threads(opt_t *opt, tm_conf_t *conf, char const *arg) {
	conf->nth = mm_atoi(arg, 0);
	oassert(opt, conf->nth < MAX_THREADS, "#threads must be less than %d.", MAX_THREADS);
}
static void tm_conf_help(opt_t *opt, tm_conf_t *conf, char const *arg) {
	_unused(opt);
	_unused(arg);
	conf->verbose++;
	conf->help++;
}

/* index filename */
static void tm_conf_idxdump(opt_t *opt, tm_conf_t *conf, char const *arg) {
	_unused(opt);
	conf->idxdump = opt_strdup(opt, arg, strlen(arg));
}
static void tm_conf_profile(opt_t *opt, tm_conf_t *conf, char const *arg) {
	_unused(opt);
	conf->profile = opt_strdup(opt, arg, strlen(arg));
}

/* indexing */
static void tm_conf_kmer(opt_t *opt, tm_conf_t *conf, char const *arg) {
	size_t k = mm_atoi(arg, 0);
	oassert(opt, k >= 3 && k <= 10, "k must be inside [3,10].");

	conf->fallback.kmer = k;
}
static void tm_conf_chain(opt_t *opt, tm_conf_t *conf, char const *arg) {
	size_t w = mm_atoi(arg, 0);
	oassert(opt, w >= 1 && w <= 64, "w must be inside [2,64].");

	conf->fallback.window = w;
}
static void tm_conf_match(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t m = mm_atoi(arg, 0);
	oassert(opt, m >= 1 && m <= 7, "match award (-a) must be inside [1,7] (%s).", arg);

	conf->fallback.match = m;
}
static void tm_conf_mismatch(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t x = mm_atoi(arg, 0);
	oassert(opt, x >= 1 && x <= 7, "mismatch penalty (-b) must be inside [1,7] (%s).", arg);

	conf->fallback.mismatch = x;
}
static void tm_conf_gap_open(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t gi = mm_atoi(arg, 0);
	oassert(opt, gi >= 1 && gi <= 7, "gap open penalty (-p) must be inside [1,7] (%s).", arg);

	conf->fallback.gap_open = gi;
}
static void tm_conf_gap_extend(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t ge = mm_atoi(arg, 0);
	oassert(opt, ge >= 1 && ge <= 7, "gap extension penalty (-q) must be inside [1,7] (%s).", arg);

	conf->fallback.gap_extend = ge;
}
static void tm_conf_min_score(opt_t *opt, tm_conf_t *conf, char const *arg) {
	_unused(opt);

	conf->fallback.min_score = mm_atoi(arg, 0);
}



static _force_inline
uint64_t tm_conf_check_sanity(tm_conf_t *conf)
{
	return(opt_ecnt(&conf->opt));
}

static _force_inline
uint64_t tm_conf_restore_default(tm_conf_t *conf)
{
	tm_conf_t defaults = {
		.fallback = {
			.kmer = 5,
			.window = 32,
			.match = 2,
			.mismatch = -3,
			.gap_open = -5,
			.gap_extend = -1,
			.min_score = 0
		},
		.print = {
			0
		}
	};

	/* we assume all element is sized 64bit */
	uint64_t *q = (uint64_t *)conf;
	uint64_t const *p = (uint64_t const *)&defaults;
	for(size_t i = 0; i < offsetof(tm_conf_t, opt) / sizeof(uint64_t); i++) {
		if(q[i] != 0) { continue; }		/* parameter set */
		if(p[i] == 0) { continue; }		/* default value not found */
		q[i] = p[i];					/* load default */
	}
	return(0);
}

static _force_inline
void tm_conf_destroy_static(tm_conf_t *conf)
{
	if(conf == NULL) { return; }
	opt_destroy_static(&conf->opt);
	return;
}

static _force_inline
uint64_t tm_conf_init_static(tm_conf_t *conf, char const *const *argv, FILE *fp)
{
	*conf = (tm_conf_t){
		/* logger */
		.log = stderr,

		#define _c(_x)			( (opt_callback_t)(_x) )
		.opt.t = {
			['\0'] = { 0, NULL },
			['d'] = { OPT_REQ,  _c(tm_conf_idxdump) },

			['v'] = { OPT_OPT,  _c(tm_conf_verbose) },
			['h'] = { OPT_BOOL, _c(tm_conf_help) },
			['t'] = { OPT_REQ,  _c(tm_conf_threads) },

			/* preset and configuration file */
			['x'] = { OPT_REQ,  _c(tm_conf_preset) },

			/* fallback parameters */
			['k'] = { OPT_REQ,  _c(tm_conf_kmer) },
			['w'] = { OPT_REQ,  _c(tm_conf_chain) },
			['a'] = { OPT_REQ,  _c(tm_conf_match) },
			['b'] = { OPT_REQ,  _c(tm_conf_mismatch) },
			['p'] = { OPT_REQ,  _c(tm_conf_gap_open) },
			['q'] = { OPT_REQ,  _c(tm_conf_gap_extend) },
			['m'] = { OPT_REQ,  _c(tm_conf_min_score) }
		}
		#undef _c
	};

	opt_init_static(&conf->opt, fp);
	if(opt_parse_argv(&conf->opt, conf, argv + 1)) { goto _tm_conf_init_fail; }
	if(tm_conf_restore_default(conf)) { goto _tm_conf_init_fail; }
	if(tm_conf_check_sanity(conf)) { goto _tm_conf_init_fail; }

	/* parsed without error */
	conf->args = opt_join(&conf->opt, argv, ' ');
	return(0);


_tm_conf_init_fail:;
	opt_destroy_static(&conf->opt);
	return(1);
}

/* determine help and verbose level */
typedef struct {
	FILE *fp;
	uint64_t help, quit;
} tm_conf_outfp_t;

static _force_inline
tm_conf_outfp_t tm_conf_get_outfp(tm_conf_t *conf)
{
	/* use stdout for explicit version (-v) and help (-h) options */
	if(conf->verbose == 1 || conf->help > 0) {
		return((tm_conf_outfp_t){
			.fp = stdout,
			.help = conf->help,
			.quit = 1
		});
	}

	/* restore default verbose level */
	if(conf->verbose == 0) {
		conf->verbose = 1;
	}

	/* implicit help invoked when no input file found; always redirected to stderr */
	return((tm_conf_outfp_t){
		.fp = conf->log,
		.help = opt_parg_cnt(&conf->opt) == 0,
		.quit = 0
	});
}

/* print help */
static _force_inline
void tm_conf_print_help(tm_conf_t const *conf, FILE *log)
{
	if(conf->verbose == 0) { return; }

	#define _msg_impl(_level, _fmt, ...) { logger_printf(log, _fmt "%s\n", __VA_ARGS__); }
	#define _msg(_level, ...) { if((_level) <= conf->verbose + 1) { _msg_impl(_level, __VA_ARGS__, ""); } }

	_msg(2, "\n"
			"  tinymasker - fast repeat masking tool\n"
			"");
	_msg(2, "Usage:\n"
			"  index construction:\n"
			"    $ tinymasker -t4 -d index.tmi repeats.fa\n"
			"  mask:\n"
			"    $ tinymasker -t4 index.tmi contigs.fa\n"
			"")
	_msg(2, "General options:");
	_msg(2, "  -x STR/FILE  load preset params or load config file []");
	_msg(2, "  -t INT       number of threads [%zu]", conf->nth);
	_msg(2, "  -d FILE      dump index to FILE (index construction mode)");
	_msg(2, "  -v [INT]     show version number or set verbose level (when number passed)");
	_msg(2, "");
	_msg(2, "Indexing options:");
	_msg(2, "  -c FILE      load profile configurations []");
	_msg(2, "");
	_msg(2, "Indexing and mapping fallbacks (default params):");
	_msg(2, "  -k INT       k-mer length [%zu]", conf->fallback.kmer);
	_msg(2, "  -w INT       chaining window size [%zu]", conf->fallback.window);
	_msg(2, "  -a INT       match award []");
	_msg(2, "  -b INT       mismatch penalty []");
	_msg(2, "  -p INT       gap-open penalty []");
	_msg(2, "  -q INT       gap-extension penalty []");
	_msg(2, "  -m INT       minimum score threshold []");
	_msg(2, "");

	#undef _msg_impl
	#undef _msg

	return;
}


/* index construction */

static _force_inline
int main_index_error(tm_conf_t *conf, int error_code, char const *filename)
{
	_unused(conf);
	switch(error_code) {
	/* argument missing */
	case ERROR_NO_ARG: error("argument is not enough. at least one reference file is required."); break;

	/* opening files */
	case ERROR_OPEN_IDX: error("failed to open index file `%s' in write mode. Please check file path and its permission.", filename); break;
	case ERROR_OPEN_RSEQ: error("failed to open sequence file `%s' for building index. Please check file path and its format.", filename); break;
	}
	return(error_code);
}

static _force_inline
int main_index_intl(tm_conf_t *conf, pg_t *pg, pt_t *pt)
{
	if(opt_parg_cnt(&conf->opt) == 0) {
		return(main_index_error(conf, ERROR_NO_ARG, NULL));
	}

	/* iterate over index blocks */
	kv_foreach(void *, opt_pargv(&conf->opt), ({

		/* generate index */
		tm_idx_t *mi = tm_idx_gen(&conf->fallback, pt, *p, conf->profile, stderr);
		if(mi == NULL) {
			/* failed to open file */
			return(main_index_error(conf, ERROR_OPEN_RSEQ, conf->idxdump));
		}

		/* dump index */
		size_t size = tm_idx_dump(mi, pg, (write_t)pgwrite);

		/* flush output for next batch */
		pg_flush(pg);

		tm_idx_destroy(mi);
		message(conf->log, "built and dumped index for `%s', on-memory size of this chunk: %.1f MB", (char const *)*p, (double)size / (1024ULL * 1024));
	}));
	return(0);
}

static _force_inline
int main_index(tm_conf_t *conf, pt_t *pt)
{
	/* add suffix if missing and if /dev/xxx */
	if(!mm_startswith(conf->idxdump, "/dev") && !mm_endswith(conf->idxdump, ".mai")) {
		message(conf->log, "index filename does not end with `.mai' (added).");
		conf->idxdump = opt_append(&conf->opt, conf->idxdump, ".mai");
	}

	/* open file in write mode */
	FILE *fp = fopen(conf->idxdump, "wb");
	if(fp == NULL) {
		/* failed open file (locked?) */
		return(main_index_error(conf, ERROR_OPEN_IDX, conf->idxdump));
	}

	/* initialize compressor */
	pg_t pg;
	pg_init_static(&pg, fp, pt);

	/* index_intl does everything */
	int error_code = main_index_intl(conf, &pg, pt);

	/* done */
	pg_stat_t stat = pg_destroy_static(&pg);
	message(conf->log, "done. total index size (compressed) on disk: %.1f MB.", (double)stat.out / (1024ULL * 1024));

	fclose(fp);
	return(error_code);
}


/* scan-and-mask */

static _force_inline
int main_scan_error(tm_conf_t *conf, int error_code, char const *file)
{
	_unused(conf);
	switch(error_code) {
	/* unknown */
	case ERROR_INTERNAL: error("failed to instanciate alignment context."); break;

	/* in mapping */
	case ERROR_OPEN_QSEQ: error("failed to open sequence file `%s'. Please check file path and format.", file); break;

	/* argument missing */
	case ERROR_NO_ARG: error("argument is not enough. a reference and at least one query files are required."); break;

	/* index loading */
	case ERROR_OPEN_IDX: error("failed to open index file `%s'. Please check file path and permission.", file); break;
	case ERROR_LOAD_IDX: error("failed to load index block from `%s'. Please check file path and version, or rebuild the index.", file); break;

	/* index construction */
	case ERROR_OPEN_RSEQ: error("failed to open sequence file `%s' for building index. Please check file path and format.", file); break;
	}
	return(error_code);
}

/* working buffer */
typedef struct {
	/* always available */
	tm_conf_t *conf;
	size_t fcnt;					/* fetched count */

	pt_t *pt;
	tm_print_t *printer;

	char const *ref;				/* reference filename */
	char const *const *query;		/* query filename */
	size_t qcnt;


	/* for prebuilt index */
	FILE *fp;				/* != NULL if prebuilt index is available */
	pg_t pg;
} main_scan_tbuf_t;

static _force_inline
int main_scan_tbuf_destroy_static(main_scan_tbuf_t *w)
{
	if(w->fp != NULL) {
		fclose(w->fp);
		pg_destroy_static(&w->pg);
	}
	return(0);
}

static _force_inline
int main_scan_tbuf_init_static(main_scan_tbuf_t *w, tm_conf_t *conf, char const *const *parg, size_t pcnt, tm_print_t *printer, pt_t *pt)
{
	/* save args */
	*w = (main_scan_tbuf_t){
		.conf = conf,
		.fcnt = 0,

		.pt = pt,
		.printer = printer,

		/* inputs */
		.ref = parg[0],
		.query = &parg[1],
		.qcnt = pcnt - 1
	};

	/* check suffix to determine if the file is prebuilt index */
	if(!mm_endswith(parg[0], ".mai")) {
		return(0);			/* not a prebuilt index */
	}

	/* open; if fails, it would be invalid path or permission */
	if((w->fp = fopen(parg[0], "rb")) == NULL) {
		return(main_scan_error(conf, ERROR_OPEN_IDX, parg[0]));
	}
	pg_init_static(&w->pg, w->fp, pt);
	return(0);
}


static _force_inline
int main_scan_idx_gen(main_scan_tbuf_t *w, tm_idx_t **pmi)
{
	if(++w->fcnt > 1) { return(0); }

	/* first call */
	tm_idx_t *mi = tm_idx_gen(&w->conf->fallback, w->pt, w->ref, w->conf->profile, stderr);
	if(mi == NULL) {
		/* failed to open file */
		return(main_scan_error(w->conf, ERROR_OPEN_RSEQ, w->ref));
	}
	message(stderr, "built index for `%s'.", w->ref);
	*pmi = mi;
	return(0);
}

static _force_inline
int main_scan_idx_load(main_scan_tbuf_t *w, tm_idx_t **pmi)
{
	if(pg_eof(&w->pg)) { return(0); }

	/* prebuilt index available; try to fetch next block */
	tm_idx_t *mi = tm_idx_load(&w->pg, (read_t)pgread);
	if(++w->fcnt == 1 && mi == NULL) {
		/* would be broken */
		return(main_scan_error(w->conf, ERROR_LOAD_IDX, w->ref));
	}

	/* NULL for tail */
	if(mi != NULL) { message(stderr, "loaded index block from `%s'.", w->ref); }
	*pmi = mi;
	return(0);
}

/* for each query file */
static _force_inline
int main_scan_foreach_qfile(main_scan_tbuf_t *w, tm_mtscan_t *mt)
{
	for(size_t i = 0; i < w->qcnt; i++) {
		if(tm_mtscan_file(mt, w->query[i])) {
			return(main_scan_error(w->conf, ERROR_OPEN_QSEQ, NULL));
		}
	}
	return(0);
}

/* for each index chunk */
static _force_inline
int main_scan_foreach_idx(main_scan_tbuf_t *w)
{
	while(1) {
		tm_idx_t *mi = NULL;

		int fetcher_error_code = (w->fp == NULL
			? main_scan_idx_gen(w, &mi)
			: main_scan_idx_load(w, &mi)
		);
		if(fetcher_error_code != 0 || mi == NULL) {
			return(fetcher_error_code);		/* error should be handled inside */
		}

		/* instanciate multithreading context for this index chunk */
		tm_mtscan_t *mt = tm_mtscan_init(mi, w->printer, w->pt);
		if(mt == NULL) {
			return(main_scan_error(w->conf, ERROR_INTERNAL, NULL));
		}

		/* do the task; for each query file */
		int scan_error_code = main_scan_foreach_qfile(w, mt);

		/* done for this index chunk */
		tm_mtscan_destroy(mt);
		tm_idx_destroy(mi);
		if(scan_error_code != 0) { return(scan_error_code); }
	}
	return(0);
}

/* for each index filename */
static _force_inline
int main_scan_intl(tm_conf_t *conf, tm_print_t *printer, pt_t *pt)
{
	char const *const *parg = (char const *const *)opt_parg(&conf->opt);
	size_t pcnt = opt_parg_cnt(&conf->opt);

	/* error if no argument is given */
	if(pcnt < 2) {
		return(main_scan_error(conf, ERROR_NO_ARG, NULL));
	}

	/* instanciate thread-local working buffers */
	main_scan_tbuf_t w;
	int init_error_code = main_scan_tbuf_init_static(&w, conf, parg, pcnt, printer, pt);
	if(init_error_code != 0) {
		return(init_error_code);
	}

	/* we consider the first argument as reference */
	int scan_error_code = main_scan_foreach_idx(&w);

	/* done */
	main_scan_tbuf_destroy_static(&w);
	return(scan_error_code);
}

/* printer */
static _force_inline
int main_scan(tm_conf_t *conf, pt_t *pt)
{
	/* instanciate alignment formatter */
	tm_print_t printer;
	tm_print_init_static(&printer, &conf->print, conf->args);

	/* dispatch */
	int error_code = main_scan_intl(conf, &printer, pt);

	tm_print_destory_static(&printer);
	return(error_code);
}

/* create worker threads for indexing and mapping */
static _force_inline
int main_dispatch(tm_conf_t *conf)
{
	pt_t *pt = pt_init(conf->nth);
	if(pt == NULL) {
		return(main_scan_error(conf, ERROR_INTERNAL, NULL));
	}

	/* dispatch either index or scan */
	int error_code = (conf->idxdump ? main_index : main_scan)(conf, pt);

	/* done */
	pt_destroy(pt);
	return(error_code);
}


/* entry */
int main(int argc, char *argv[])
{
	int error_code = ERROR_NO_ARG;

	/* unittest dispatcher */
	#if defined(UNITTEST) && UNITTEST != 0
		if(argc > 1 && strcmp(argv[1], "unittest") == 0) {
			return(unittest_main(argc, argv));
		}
	#else
		_unused(argc);
	#endif

	/* instanciate option object */
	logger_init();
	tm_conf_t conf;
	tm_conf_init_static(&conf, (char const *const *)argv, stderr);

	/* always print version */
	tm_conf_outfp_t out = tm_conf_get_outfp(&conf);
	message(out.fp, "Version: %s, Build: %s", tm_version(), MK_ARCH_NAME);

	/* when -h is passed or no input file is given, print help message */
	if(out.help) {
		if(conf.help > 0) { error_code = 0; }	/* also exit status is 0 (not an error) */
		tm_conf_print_help(&conf, out.fp);		/* we use stdout when invoked by -h option */
	}
	if(out.quit) { goto _main_final; }

	/* dispatch tasks and get return code */
	if((error_code = main_dispatch(&conf)) == 0) {
		message(conf.log, "Command: %s", conf.args);	/* print log when succeeded */
	}
	logger_destroy();

_main_final:;
	tm_conf_destroy_static(&conf);
	return(error_code);
}

/**
 * end of tinymasker.c
 */
