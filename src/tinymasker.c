// #define DEBUG
// #define DEBUG_KMER
#define _printu_v16i8(a) { \
	uint8_t _buf[16]; \
	_storeu_v16i8(_buf, a); \
	debug("(v16i8_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_buf[15] - 128, \
		_buf[14] - 128, \
		_buf[13] - 128, \
		_buf[12] - 128, \
		_buf[11] - 128, \
		_buf[10] - 128, \
		_buf[9] - 128, \
		_buf[8] - 128, \
		_buf[7] - 128, \
		_buf[6] - 128, \
		_buf[5] - 128, \
		_buf[4] - 128, \
		_buf[3] - 128, \
		_buf[2] - 128, \
		_buf[1] - 128, \
		_buf[0] - 128); \
}

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

#define DZ_WRAPPED_API			( 0 )
#define DZ_CIGAR_OP				0x44494d4d
#define DZ_REF_MAT_SIZE			( 16 )
#define DZ_QUERY_MAT_SIZE		( 4 )
#define DZ_TRANSPOSE_MATRIX		( 1 )
#include "dozeu.h"


/* toml parser and regex */
#include "toml.h"
#include "re.h"


/* misc */
#ifndef TM_VERSION
#  define TM_VERSION			"tinymasker-0.0.1"
#endif
#define TM_ARCH_NAME			ARCH_NAME			/* SSE4.1 or AVX2 */

/* version string */
static _force_inline
char const *tm_version(void)
{
	char const *prefix = "tinymasker-";

	/* remove prefix */
	uint64_t spos = mm_startswith(TM_VERSION, prefix) ? mm_strlen(prefix) : 0;
	return(&TM_VERSION[spos]);
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

/* shifted by 2bits */
enum alphabet_2bit {
	nA = 0x00,
	nC = 0x04,
	nG = 0x08,
	nT = 0x0c
};

/* FASTA/Q parser */
#include "bseq.h"


/* SIMD supplementary */
#define _mm_blendv_epi32(x, y, z) ({ \
	__m128 const _x = _mm_castsi128_ps(x); \
	__m128 const _y = _mm_castsi128_ps(y); \
	__m128 const _z = _mm_castsi128_ps(z); \
	__m128 const _w = _mm_blendv_ps(_x, _y, _z); \
	_mm_castps_si128(_w); \
})
#define _mm_shuffle2_epi32(x, y, z) ({ \
	__m128 const _x = _mm_castsi128_ps(x); \
	__m128 const _y = _mm_castsi128_ps(y); \
	__m128 const _w = _mm_shuffle_ps(_x, _y, z); \
	_mm_castps_si128(_w); \
})


/* 1-origin coorrdinates for seed positions */
#define TM_SEED_POS_OFS				( 1 )


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

/* for debugging */
static char *tm_kmer_to_str(uint32_t const *kmer)
{
	char buf[12] = { 0 };
	for(size_t i = 0; i < 11; i++) {
		buf[i] = "ACGT"[(*kmer>>(2 * (10 - i))) & 0x03];
	}
	return(strdup(buf));
}


/* construct reference object */
typedef struct {
	uint32_t next : 8;
	uint32_t pos  : 24;
	uint32_t kmer;
} tm_ref_kpos_t;
_static_assert(sizeof(tm_ref_kpos_t) == 8);

typedef struct { tm_ref_kpos_t *a; size_t n, m; } tm_ref_kpos_v;

#define tm_ref_kmer(x)			( _loadu_u64(&x) )
// KRADIX_SORT_INIT(kmer, tm_ref_kpos_t, tm_ref_kmer, 4);
KRADIX_SORT_INIT(kmer, tm_ref_kpos_t, tm_ref_kmer, 8);	/* sort all */

typedef struct {
	uint32_t f, r;
} tm_ref_kmer_pair_t;

typedef struct {
	/* constants */
	uint32_t amask, shift;
	// uint32_t kofs;

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
	size_t const prev_size = w->size;

	/* update size and branch pos */
	w->size = w->size<<1;
	w->branches |= (0x01<<w->shift);

	for(size_t i = 0; i < prev_size; i++) {
		w->kmer[prev_size + i] = w->kmer[i];
	}
	return(prev_size);
}

static _force_inline
uint64_t tm_ref_update_stack(tm_ref_work_t *w)
{
	uint64_t const shrink = (w->branches & 0x01) != 0;
	w->branches = w->branches>>2;
	return(shrink);
}

static _force_inline
void tm_ref_shrink_stack(tm_ref_work_t *w)
{
	debug("shrink");
	w->size = w->size>>1;
	for(size_t i = 0; i < w->size; i++) {
		w->kmer[i] = w->kmer[2 * i];
	}
	return;
}

static _force_inline
void tm_ref_update_kmer(tm_ref_kmer_pair_t *kmer, size_t size, uint32_t amask, uint32_t shift, uint8_t base)
{
	uint32_t const f = base<<6, r = (base ^ 0x03)<<shift;
	for(size_t i = 0; i < size; i++) {
		kmer[i].f = ((kmer[i].f<<2) | f) & amask;
		kmer[i].r = ((kmer[i].r>>2) | r) & amask;
	}
	return;
}

static _force_inline
v4i32_t tm_ref_compose_kpos(v2i32_t pos, v2i32_t kmer)
{
	static uint8_t const shuf[16] __attribute__(( aligned(16) )) = {
		0, 4, 5, 6, 1, 2, 3, 0x80, 8, 12, 13, 14, 9, 10, 11, 0x80
	};

	/* fw in lower, rv in upper */
	v16i8_t const sv = _load_v16i8(shuf);
	v16i8_t const v = (v16i8_t){
		_mm_unpacklo_epi32(kmer.v1, pos.v1)
	};
	v16i8_t const w = _shuf_v16i8(v, sv);

	// _print_v16i8(v);
	// _print_v16i8(w);
	return((v4i32_t){ w.v1 });
}

#if 1
static _force_inline
void tm_ref_push_pos(tm_ref_work_t *w, v2i32_t pos)
{
	for(size_t i = 0; i < w->size; i++) {
		v2i32_t const kmer = _loadu_v2i32(&w->kmer[i]);
		v4i32_t const kpos = tm_ref_compose_kpos(pos, kmer);
		_storeu_v4i32(w->q, kpos);
		w->q += 2;
	}
	return;
}
#else
static _force_inline
void tm_ref_push_pos(tm_ref_work_t *w, v2i32_t pos)
{
	for(size_t i = 0; i < w->size; i++) {
		*w->q++ = (tm_ref_kpos_t){
			.kmer = w->kmer[i].f,
			.pos  = pos + TM_SEED_POS_OFS				/* pos + base_offset */
		};
		*w->q++ = (tm_ref_kpos_t){
			.kmer = w->kmer[i].r,
			.pos  = w->kofs - (pos + TM_SEED_POS_OFS)	/* (k - 1) - (pos + base_offset) */
		};
		debug("pos(%u), kmer(%x, %x)", pos, w->kmer[i].f, w->kmer[i].r);
	}
	return;
}
#endif

static _force_inline
void tm_ref_push_base(tm_ref_work_t *w, uint8_t base)
{
	size_t const size = w->size;		/* save current size before duplicating stack */
	uint8_t const b2 = tm_ref_base_to_2bit(base);

	if(_unlikely(tm_ref_is_ambiguous(base))) {
		debug("ambiguous");
		tm_ref_dup_stack(w);
		tm_ref_update_kmer(&w->kmer[size], size, w->amask, w->shift, b2>>2);
	}

	/* push base */
	tm_ref_update_kmer(&w->kmer[0], size, w->amask, w->shift, b2 & 0x03);	
	// debug("branches(%lx), b2(%x)", w->branches, b2);
	return;
}

static _force_inline
tm_ref_kpos_t *tm_ref_reserve(tm_ref_kpos_v *buf, tm_ref_kpos_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);
	// debug("cnt(%zu)", kv_cnt(*buf));

	size_t const kcnt = kv_cnt(*buf);
	tm_ref_kpos_t *p = kv_reserve(tm_ref_kpos_t, *buf, kcnt + 1024);
	return(&p[kcnt]);
}

static _force_inline
size_t tm_ref_enumerate_kmer(size_t kbits, uint8_t const *seq, size_t slen, tm_ref_kpos_v *buf)
{
	tm_ref_work_t w = {
		/* place kmers at 6..kbits + 6 */
		.amask = ((0x01<<kbits) - 0x01)<<6,
		.shift = kbits - 2 + 6,
		// .kofs  = (kbits>>1),		/* k - 1 */

		/* make room in dst array */
		.q = tm_ref_reserve(buf, kv_ptr(*buf)),

		/* clear stack; keeps all the combination of ambiguous bases */
		.size = 1,
		.branches = 0,
		.kmer = { { 0 } }
	};
	// debug("kbits(%zu), seq(%p), slen(%zu)", kbits, seq, slen);

	size_t const k = kbits>>1;
	v2i32_t const inc = _seta_v2i32(-1, 1);		/* rv, fw */
	v2i32_t pos = _seta_v2i32(k - TM_SEED_POS_OFS, TM_SEED_POS_OFS);

	/* push bases at the head */
	for(size_t i = 0; i < k - 1; i++) {
		tm_ref_push_base(&w, seq[i]);
		pos = _add_v2i32(pos, inc);
	}

	/* body */
	for(size_t i = k - 1; i < slen; i++) {
		tm_ref_push_base(&w, seq[i]);			/* stack will be expanded if needed */

		/* update offset; position of the next base to compose [,) range */
		pos = _add_v2i32(pos, inc);
		tm_ref_push_pos(&w, pos);

		/* update branch stack */
 		if(_unlikely(tm_ref_update_stack(&w))) {
			tm_ref_shrink_stack(&w);			/* shrink if needed */
		}

		/*
		 * check need for expansion every 32times
		 * (Skylake's branch predictor successfully learns once-every-32times event, but fails for more)
		 */
		if((i & 0x1f) != 0) { continue; }
		w.q = tm_ref_reserve(buf, w.q);
	}

	/* tail */
	tm_ref_push_base(&w, '\0');

	pos = _add_v2i32(pos, inc);
	tm_ref_push_pos(&w, pos);

	/* done; record count */
	kv_cnt(*buf) = w.q - kv_ptr(*buf);
	return(kv_cnt(*buf));
}

static _force_inline
size_t tm_ref_dedup_kmer(tm_ref_kpos_t *kptr, size_t kcnt)
{
	tm_ref_kpos_t const *p = kptr, *t = &kptr[kcnt];	/* src */
	tm_ref_kpos_t *q = kptr + 1;						/* dst */

	uint64_t const mask = ~0xffULL;		/* clear out next field */
	uint64_t prev = _loadu_u64(p) & mask;
	while(++p < t) {
		uint64_t const kmer = _loadu_u64(p);
		uint64_t const curr = kmer & mask;

		if(curr == prev) {
			debug("dedup, i(%zu), prev(%lx), curr(%lx, %lx)", p - kptr, prev, curr, kmer);
			continue;
		}		/* skip pos and kmer are the same */

		/* pos or kmer changed; save */
		_storeu_u64(q, kmer);
		q++;

		/* update for next */
		prev = curr;
	}
	return(q - kptr);
}

static _force_inline
size_t tm_ref_collect_kmer(tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen, tm_ref_kpos_v *buf)
{
	kv_clear(*buf);

	/* slice k-mers (with additional succeeding base) */
	if(tm_ref_enumerate_kmer(conf->kbits + 2, seq, slen, buf) == 0) {
		return(0);				/* abort if no kmer found */
	}

	/* sort kmers by position */
	tm_ref_kpos_t *kptr = kv_ptr(*buf);
	size_t const kcnt = kv_cnt(*buf);

	radix_sort_kmer(kptr, kcnt);
/*
	debug("kcnt(%zu)", kcnt);
	for(size_t i = 0; i < kcnt; i++) {
		debug("i(%zu), kmer(%x), pos(%x)", i, kptr[i].kmer, kptr[i].pos);
	}
*/
	/* dedup kmers around ambiguous bases */
	debug("kcnt(%zu)", kcnt);
	size_t const dcnt = tm_ref_dedup_kmer(kptr, kcnt);
	kv_cnt(*buf) = dcnt;

/*
	debug("kcnt(%zu)", dcnt);
	for(size_t i = 0; i < dcnt; i++) {
		debug("i(%zu), kmer(%x), pos(%x)", i, kptr[i].kmer, kptr[i].pos);
	}
*/
	return(dcnt);
}

/* pack kmer-position array */
typedef struct {
	int32_t next  : 16;			/* diff to next bin */
	uint32_t tail : 12;			/* pos count for this bin */
	uint32_t plen : 4;			/* unmatching prefix length */
} tm_ref_link_t;
_static_assert(sizeof(tm_ref_link_t) == sizeof(uint32_t));

typedef struct {
	tm_ref_link_t link[4];		/* de-Bruijn-graph-like link for each nucleotide */

	#ifdef DEBUG_KMER
		uint32_t kmer;
	#endif

	uint16_t pos[];
} tm_ref_bin_t;
#ifndef DEBUG_KMER
_static_assert(sizeof(tm_ref_bin_t) == 16);		/* smallest bin size == 16 */
#endif

#define TM_REF_ALIGN_SIZE		( sizeof(uint64_t) )
#define TM_REF_BASE_OFS			( sizeof(tm_ref_sketch_t) )


typedef struct {
	uint8_t covered, exist;		/* covered flag for feeder */
	uint16_t cnt;
	uint32_t ofs;
} tm_ref_cnt_t;
typedef struct { tm_ref_cnt_t *a; size_t n, m; } tm_ref_cnt_v;

static _force_inline
size_t tm_ref_count_kmer(tm_ref_cnt_t *p, size_t ksize, tm_ref_kpos_t const *kpos, size_t kcnt)
{
	_unused(ksize);

	// debug("ksize(%zx), kcnt(%zu)", ksize, kcnt);
	for(tm_ref_kpos_t const *k = kpos, *t = &kpos[kcnt]; k < t; k++) {
		// tm_ref_cnt_t *q = &p[k->kmer>>2];
		tm_ref_cnt_t *q = &p[k->kmer];

		q->covered |= 1;
		q->exist   |= 1;
		q->cnt++;
	}
	return(kcnt);
}

static _force_inline
void tm_ref_patch_feeder(tm_ref_cnt_t *p, size_t ksize)
{
	for(size_t i = 0; i < ksize; i++) {
		// debug("i(%r), ksize(%zx)", tm_kmer_to_str, &i, ksize - 1);
		if(p[i].exist == 0) { continue; }

		for(uint64_t mask = (ksize - 1)>>2; mask != 0; mask >>= 2) {
			uint64_t krem = i & mask;
			// debug("test mask(%lx), krem(%r), covered(%u), exist(%u)", mask, tm_kmer_to_str, &krem, p[krem].covered, p[krem].exist);
			// if(p[krem].covered) { break; }

			p[krem].covered |= 1;
		}
	}

	size_t patched = 0;
	for(size_t i = 0; i < (ksize>>2); i++) {
		if(p[i].covered == 0) {
			p[i].exist |= 1;
			patched++;
		}
	}
	debug("ksize(%zu), patched(%zu), ratio(%f)", ksize, patched, 4.0 * (double)patched / (double)ksize);
	return;
}

static _force_inline
size_t tm_ref_calc_size(tm_ref_cnt_t const *p, size_t ksize)
{
	size_t const hdr = sizeof(tm_ref_bin_t);

	size_t acc = TM_REF_BASE_OFS;
	for(size_t i = 0; i < ksize; i++) {
		if(p[i].exist == 0) { continue; }
		acc += _roundup(hdr + sizeof(uint16_t) * p[i].cnt, TM_REF_ALIGN_SIZE);
		// debug("i(%zx), cnt(%zu), acc(%zu)", i, p[i].cnt, acc);
	}
	return(acc);
}

static _force_inline
size_t tm_ref_pack_kpos(tm_ref_cnt_t *p, size_t ksize, tm_ref_kpos_t const *kpos, size_t kcnt, tm_ref_sketch_t *sk)
{
	_unused(ksize);

	uint64_t const kmask = 0xffffffff000000ff;	/* extract kmer and next fields */

	/* first bin */
	size_t const hdr = sizeof(tm_ref_bin_t);
	size_t ofs = TM_REF_BASE_OFS;

	tm_ref_kpos_t const *k = kpos, *t = &kpos[kcnt];
	for(size_t i = 0; i < ksize && k < t; i++) {
		if(p[i].exist == 0) { continue; }
		while(k < t && k->kmer < i) {
			debug("something is wrong, k(%p, %p), kmer(%x)", k, t, k->kmer);
			k++;
		}

		/* save current offset */
		p[i].ofs = ofs;

		/* slice bin from current offset */
		tm_ref_bin_t *bin = _add_offset(sk, ofs);
		_storeu_v4i32(bin, _zero_v4i32());		/* clear link */

		/* pack pos array */
		#if 1
		size_t x = 0, y = 0;
		while(k < t && k->kmer == i) {
			uint64_t const next = k->next>>6;
			uint64_t const kall = _loadu_u64(k) & kmask;
			while(y < next) {
				bin->link[y++].tail = x;
			}
			while(k < t && (_loadu_u64(k) & kmask) == kall) {
				bin->pos[x++] = k++->pos;
			}
		}
		/* fill missing tail */
		while(y < 4) { bin->link[y++].tail = x; }


		#else
		size_t x = 0, y = 0;
		while(k < t && (k->kmer>>2) == i) {
			uint32_t const kall = k->kmer;
			while(y < (kall & 0x03)) {
				bin->link[y++].tail = x;
			}
			while(k < t && k->kmer == kall) {
				bin->pos[x++] = k++->pos;
			}
		}
		/* fill missing tail */
		while(y < 4) { bin->link[y++].tail = x; }
		#endif
/*
		debug("kmer(%r), ofs(%zu), cnt(%u, %zu, %zu), tail(%u, %u, %u, %u)",
			tm_kmer_to_str, &i, ofs, p[i].cnt, x, y,
			bin->link[0].tail, bin->link[1].tail, bin->link[2].tail, bin->link[3].tail
		);
*/

		/* update offset */
		_storeu_v16i8(&bin->pos[x], _zero_v16i8());		/* suppress use-of-uninitialized-value error in valgrind */
		ofs += _roundup(hdr + sizeof(uint16_t) * p[i].cnt, TM_REF_ALIGN_SIZE);
	}
	return(ofs);
}

typedef struct {
	uint32_t len, ofs;

	/* for debugging */
	#ifdef DEBUG_KMER
		uint32_t kmer;
	#endif
} tm_ref_prefix_t;

static _force_inline
tm_ref_prefix_t tm_ref_find_link(tm_ref_cnt_t const *p, uint32_t kmer, size_t max_shift)
{
	/* breadth-first search on kmer branching tree over profile table */
	size_t const ksize = 0x01<<max_shift;
	size_t shift = max_shift;		/* where (max_shift - shift) / 2 is substitution length (zero for the first iteration) */
	do {
		size_t const kpitch = 0x01<<shift;
		uint32_t const krem = kmer & (kpitch - 0x01);

		/* for each substituted kmer */
		for(size_t i = 0; i < ksize; i += kpitch) {
			uint32_t const ksub = krem + i;
			if(p[ksub].exist == 0) { continue; }

			/* hit */
			// debug("max_shift(%zu), shift(%zu), kmer(%x), ksub(%x), len(%zu)", max_shift, shift, kmer, ksub, (max_shift - shift) / 2);
			return((tm_ref_prefix_t){
				#ifdef DEBUG_KMER
					.kmer = ksub,
				#endif

				.len = (max_shift - shift) / 2,
				.ofs = p[ksub].ofs
			});
		}
	} while((shift -= 2) != 0);

	/* never reach here?? */
	return((tm_ref_prefix_t){
		.len = max_shift / 2,
		.ofs = TM_REF_BASE_OFS
	});
}

static _force_inline
uint64_t tm_ref_build_link(tm_ref_cnt_t const *p, size_t ksize, tm_ref_sketch_t *sk)
{
	size_t const max_shift = _tzcnt_u64(ksize);
	uint64_t const mask = ksize - 1;

	for(size_t i = 0; i < ksize; i++) {
		if(p[i].exist == 0) { continue; }

		tm_ref_bin_t *bin = _add_offset(sk, p[i].ofs);
		#ifdef DEBUG_KMER
			bin->kmer = i;
		#endif

		/* find next link for each base */
		for(size_t j = 0; j < 4; j++) {
			tm_ref_prefix_t n = tm_ref_find_link(p, ((i<<2) + j) & mask, max_shift);
			int32_t const next = (int32_t)(n.ofs - p[i].ofs) / (int32_t)TM_REF_ALIGN_SIZE;	/* overflows; keep sign bit */
			if(next > INT16_MAX || next < INT16_MIN) {
				error("link overflow: i(%zx), next(%d)", i, next);
				return(1);
			}

			/* save link; preserve tail */
			bin->link[j].plen = n.len;
			bin->link[j].next = next;
/*
			debugblock({
				uint32_t k = i<<2;
				debug("(%r, %u) -- (%u) -> (%r, %u), len(%u), next(%d)",
					tm_kmer_to_str, &k, p[i].ofs, j,
					tm_kmer_to_str, &n.kmer, n.ofs,
					bin->link[j].plen, bin->link[j].next
				);
				// debug("done, size(%u), head_margin(%u), kbits(%u), slen(%u)", sk->size, sk->head_margin, sk->kbits, sk->slen);
			});
*/
		}
	}
	return(0);
}

static _force_inline
tm_ref_sketch_t *tm_ref_build_index(tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen, tm_ref_kpos_t const *kpos, size_t kcnt, tm_ref_cnt_v *buf)
{
	_unused(seq);
	kv_clear(*buf);

	/* reserve profile buffer */
	size_t const ksize = 0x01ULL<<conf->kbits;
	tm_ref_cnt_t *p = kv_reserve(tm_ref_cnt_t, *buf, ksize);
	memset(p, 0, ksize * sizeof(tm_ref_cnt_t));

	/* count kmers */
	tm_ref_count_kmer(p, ksize, kpos, kcnt);
	_scan_memory(p, sizeof(tm_ref_cnt_t) * ksize);

	/* patch lead */
	tm_ref_patch_feeder(p, ksize);

	/* accumulate block size */
	size_t size = tm_ref_calc_size(p, ksize);
	if(size > UINT32_MAX) { trap(); }

	/* malloc; use entire page */
	size_t const rounded_size = _roundup(size + conf->margin.head + conf->margin.tail, 4096);
	void *base = malloc(rounded_size);
	if(base == NULL) { debug("failed"); return(NULL); }

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
	if(tm_ref_build_link(p, ksize, sk)) {
		/* failed */
		free(base);
		return(NULL);
	}
	_scan_memory(sk, sk->size);


	/* anything else? */
	debug("done, size(%u), head_margin(%u), kbits(%u), slen(%u)", sk->size, sk->head_margin, sk->kbits, sk->slen);
	return(sk);
}


/* working buffer */
typedef struct {
	tm_ref_kpos_v kpos;
	tm_ref_cnt_v prof;
} tm_ref_tbuf_t;

static _force_inline
void tm_ref_destroy_tbuf_static(tm_ref_tbuf_t *tbuf)
{
	kv_destroy(tbuf->kpos);
	kv_destroy(tbuf->prof);
	return;
}

static _force_inline
tm_ref_sketch_t *tm_ref_sketch(tm_ref_tbuf_t *self, tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen)
{
	if(slen > INT16_MAX) { return(NULL); }

	/* slice kmers from sequence (ambiguous bases are expanded here) */
	if(tm_ref_collect_kmer(conf, seq, slen, &self->kpos) == 0) {		/* with additional base */
		debug("failed");
		return(NULL);
	}

	/* build index */
	tm_ref_kpos_t *kpos = kv_ptr(self->kpos);
	size_t const kcnt = kv_cnt(self->kpos);
	return(tm_ref_build_index(conf, seq, slen, kpos, kcnt, &self->prof));
}


/* matcher */
typedef struct {
	// size_t dst, src, cnt;
	v16i8_t v;
} tm_ref_squash_t;

typedef struct {
	uint64_t prefix, unmatching;
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
		.unmatching = -1LL,
		.bin        = _add_offset(ref, TM_REF_BASE_OFS)	/* bin for AAA... */
	});
}

#if 0
static _force_inline
tm_ref_squash_t tm_ref_load_squash(tm_ref_bin_t const *bin, uint64_t unmatching, uint8_t next)
{
	static int8_t const scatter[16] __attribute__(( aligned(16) )) = {
		0, 0, 0, 0, 0x80, 0x80, 0x80, 0x80,
		0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80
	};
	static int8_t const shuf[16] __attribute__(( aligned(16) )) = {
		-2, -1, 2, 3, 14, 15, 0x80, 0x80,
		0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80
	};

	v16i8_t const tv = _load_v16i8(scatter);
	v16i8_t const sv = _load_v16i8(shuf);

	v2i64_t const nv = _seta_v2i64(0, (unmatching & 0xc0) | next);		/* 0xff when unmatching; clears src and dst columns */
	v16i8_t const mv = _shuf_v16i8(nv, tv);

	v16i8_t const w = _loadu_v16i8(bin->link);
	v16i8_t const v = _shuf_v16i8(w, _add_v16i8(mv, sv));
	return((tm_ref_squash_t){ .v = v });
}
#else
static _force_inline
tm_ref_squash_t tm_ref_load_squash(tm_ref_bin_t const *bin, uint64_t unmatching, uint8_t next)
{
	static int8_t const shuf[16] __attribute__(( aligned(16) )) = {
		-2, -1, 2, 3, 14, 15, 0x80, 0x80,
		0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80
	};
	v16i8_t const sv = _load_v16i8(shuf);

	uint64_t const mask = (0x01010101 * next) + (0xf0f0f0f0 & unmatching);
	v2i64_t const mv = _seta_v2i64(0, mask);		/* movq */
	v16i8_t const iv = _add_v16i8(mv, sv);

	v16i8_t const w = _loadu_v16i8(bin->link);
	v16i8_t const v = _shuf_v16i8(w, iv);
	return((tm_ref_squash_t){ .v = v });
}
#endif

static _force_inline
tm_ref_next_t tm_ref_match_next(tm_ref_state_t s, uint64_t keep, uint8_t next)
{
	int64_t const pitch = TM_REF_ALIGN_SIZE;

	/* get link to the next bin */
	tm_ref_link_t const *link = _add_offset(s.bin->link, next);
	tm_ref_bin_t const *bin = _add_offset(s.bin, pitch * link->next);	/* signed */

	/* update matching state */
	uint64_t const prefix     = MAX2(s.prefix + s.unmatching, link->plen);
	uint64_t const unmatching = 0ULL - (prefix != 0);

	/* calc squashable subbin of the previous node */
	tm_ref_squash_t sq = tm_ref_load_squash(s.bin,
		unmatching - keep,		/* do not squash the previous if the current bin is unmatching; also when keep == 1 */
		next
	);

	/* done */
	tm_ref_next_t r = {
		.state = {
			.prefix     = prefix,
			.unmatching = unmatching,
			.bin        = bin
		},
		.squash = sq
	};
/*
	debug("unmatching(%lu), prefix(%zu), bin(%p)", r.state.unmatching, r.state.prefix, r.state.bin);
	debug("plen(%u, %u, %u, %u), tail(%u, %u, %u, %u)",
		bin->link[0].plen, bin->link[1].plen, bin->link[2].plen, bin->link[3].plen,
		bin->link[0].tail, bin->link[1].tail, bin->link[2].tail, bin->link[3].tail
	);
*/
	return(r);
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
	size_t kmer, window, min_scnt;		/* k-mer length and chain window size */

	/* extension params */
	uint64_t match, mismatch, gap_open, gap_extend;
	int64_t min_score;
} tm_idx_conf_t;

/* matcher object; k-mer length, score matrix, and score thresholds */
typedef struct {
	uint32_t size;

	/* k-mer size */
	uint32_t kbits;

	/* chaining */
	struct {
		union {
			struct { uint32_t u, v; } sep;
			uint64_t all;
		} window;

		uint32_t kadj[2];
		uint32_t min_scnt;
		uint32_t squash_intv;
	} chain;

	/* filtering */
	struct {
		uint8_t score_matrix[16];	/* match-mismatch score matrix */
		uint8_t gap[16];			/* gap */
		uint8_t init[16];			/* p = 0 */

		/* extension test if < test_cnt, do inward extension if longer than uspan_thresh */
		uint16_t test_cnt;
		uint16_t uspan_thresh;
		int32_t min_score;			/* discard if sum of extension scores is smaller than min_score */
	} filter;

	/* metadata */
	struct {
		char *name;
	} meta;

	/* extension */
	struct {
		dz_profile_t *dz;
		int32_t min_score;
		uint8_t giv, gev, gih, geh;
		int8_t score_matrix[DZ_QUERY_MAT_SIZE * DZ_REF_MAT_SIZE];
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
#define tm_idx_ref_seq_ptr(x)		( (uint8_t const *)_add_offset((x), (x)->h.sofs) )
#define tm_idx_ref_seq_len(x)		( (x)->ref.slen )
#define tm_idx_ref_name_ptr(x)		( (uint8_t const *)_add_offset((x), sizeof(tm_idx_sketch_hdr_t) + (x)->ref.size) )
#define tm_idx_ref_name_len(x)		( strlen((char const *)tm_idx_ref_name_ptr(x)) )


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
void *tm_idx_malloc(void *unused, size_t size)
{
	_unused(unused);
	return(malloc(size));
}

static _force_inline
void tm_idx_free(void *unused, void *ptr)
{
	_unused(unused);

	free(ptr);
	return;
}

static _force_inline
void tm_idx_destroy_profile(tm_idx_profile_t **arr, size_t cnt)
{
	dz_destructor_t f = {
		.ctx = NULL,
		.fp  = (dz_free_t)tm_idx_free
	};
	for(size_t i = 0; i < cnt; i++) {
		dz_destroy_profile(&f, arr[i]->extend.dz);
		free(arr[i]);
		arr[i] = NULL;
	}
	return;
}

static _force_inline
void tm_idx_destroy_sketch(tm_idx_sketch_t **arr, size_t cnt)
{
	for(size_t i = 0; i < cnt; i++) {
		tm_ref_destroy(&arr[i]->ref);
	}
	return;
}

static _force_inline
void tm_idx_destroy_mono(tm_idx_t *mi)
{
	dz_destructor_t f = {
		.ctx = NULL,
		.fp  = (dz_free_t)tm_idx_free
	};

	/* still we need to destroy profile */
	for(size_t i = 0; i < mi->profile.cnt; i++) {
		dz_destroy_profile(&f, mi->profile.arr[i]->extend.dz);
	}
	free(mi->base);
	return;
}

static _force_inline
void tm_idx_destroy_normal(tm_idx_t *mi)
{
	if(mi->filename.ref) {
		free((void *)mi->filename.ref);
	}
	if(mi->filename.profile) {
		free((void *)mi->filename.profile);
	}

	tm_idx_destroy_profile(mi->profile.arr, mi->profile.cnt);
	tm_idx_destroy_sketch(mi->sketch.arr, mi->sketch.cnt);

	free(mi->profile.arr);
	free(mi->sketch.arr);

	free(mi);
	return;
}

static _force_inline
void tm_idx_destroy(tm_idx_t *mi)
{
	if(mi == NULL) { return; }

	if(mi->base != NULL) {
		/* monolithic index loaded by tm_idx_load */
		tm_idx_destroy_mono(mi);
		return;
	}

	/* built by tm_idx_gen; not monolithic */
	tm_idx_destroy_normal(mi);
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
	return(kv_cnt(mii->score.profile) - 1);
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
	for(size_t q = 0; q < DZ_QUERY_MAT_SIZE; q++) {
		for(size_t r = 0; r < DZ_REF_MAT_SIZE; r++) {
			fp(opaque, &profile->extend.score_matrix[q * DZ_REF_MAT_SIZE + r], q, r);
		}
	}
	return;
}

static
void tm_idx_fill_score_callback(int64_t *s, int8_t *p, size_t q, size_t r)
{
	/* query in 2bit and reference in 4bit */
	size_t const q4 = (0x8421>>(4 * q)) & 0x0f;
	*p = s[(q4 & r) == 0];
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
uint32_t tm_idx_calc_kadj(uint32_t k)
{
	return(k + 1);
}

static _force_inline
void tm_idx_fill_default(tm_idx_profile_t *profile)
{
	/* default params */
	size_t const k = 5;
	profile->kbits = 2 * k;

	profile->chain.window.sep.u = 32;
	profile->chain.window.sep.v = 32;
	profile->chain.kadj[0] = tm_idx_calc_kadj(k);
	profile->chain.kadj[1] = tm_idx_calc_kadj(k);
	profile->chain.min_scnt = 4;

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
		profile->chain.kadj[0] = tm_idx_calc_kadj(conf->kmer);
		profile->chain.kadj[1] = tm_idx_calc_kadj(conf->kmer);
	}

	/* chaining */
	if(conf->window > 0) {
		profile->chain.window.sep.u = conf->window;
		profile->chain.window.sep.v = conf->window;
	}
	if(conf->min_scnt > 0) {
		profile->chain.min_scnt = conf->min_scnt;
	}

	/* score matrix */
	if(conf->match > 0 && conf->mismatch > 0) {
		int64_t const m = conf->match;
		int64_t const x = -((int64_t)conf->mismatch);

		debug("m(%ld), x(%ld)", m, x);
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
	int64_t const mave = (int64_t)c[0].acc / (int64_t)c[0].cnt;
	int64_t const m = (c[0].min == mave
		? c[0].min
		: MAX2(1, mave * 3 / 4)
	);
	debug("(%ld, %ld, %ld, %zu)", c[0].acc, c[0].max, c[0].min, c[0].cnt);

	int64_t const xave = (int64_t)c[1].acc / (int64_t)c[1].cnt;
	int64_t const x = (c[1].max == xave
		? c[1].max
		: MIN2(-1, xave * 3 / 4)
	);
	debug("(%ld, %ld, %ld, %zu)", c[1].acc, c[1].max, c[1].min, c[1].cnt);

	debug("mave(%ld), m(%ld), xave(%ld), x(%ld)", mave, m, xave, x);

	for(size_t i = 0; i < 16; i++) {
		profile->filter.score_matrix[i] = (i == 0) ? x : m;		/* signed */
	}
	return;
}

static _force_inline
void tm_idx_calc_filter_gap(tm_idx_profile_t *profile)
{
	int32_t const ge = (profile->extend.gev + profile->extend.geh) / 2;
	v16i8_t const gv = _set_v16i8(-ge);			/* signed; negated without offset */
	_storeu_v16i8(profile->filter.gap, gv);

	debug("ge(%d)", (int8_t)ge);
	_print_v16i8(gv);
	return;
}

static _force_inline
void tm_idx_calc_filter_ivec(tm_idx_profile_t *profile)
{
	int32_t const s = 2 * profile->filter.gap[1] - profile->filter.score_matrix[1];
	for(size_t i = 0; i < 8; i++) {
		profile->filter.init[8 + i] = 128 + s *  i;
		profile->filter.init[7 - i] = 128 + s * (i + 1);
	}
	return;
}

static _force_inline
void tm_idx_calc_filter_thresh(tm_idx_profile_t *profile)
{
	uint32_t const test_cnt = profile->chain.min_scnt * 3;
	int32_t const min_score = profile->extend.min_score / 4;
	uint32_t const uspan_thresh = profile->chain.window.sep.u * 2;

	profile->filter.test_cnt  = test_cnt;
	profile->filter.min_score = MAX2(0, min_score);
	profile->filter.uspan_thresh = uspan_thresh;
	return;
}

static _force_inline
void tm_idx_calc_filter_params(tm_idx_profile_t *profile)
{
	tm_idx_calc_filter_score(profile);
	tm_idx_calc_filter_gap(profile);
	tm_idx_calc_filter_ivec(profile);
	tm_idx_calc_filter_thresh(profile);
	return;
}

static _force_inline
void tm_idx_finalize_profile(tm_idx_profile_t *profile)
{
	/* calc squash interval */
	uint32_t const u = profile->chain.window.sep.u;
	profile->chain.squash_intv = 0x40000000>>_lzcnt_u32(u);		/* divide by 2 */
	debug("wu(%u), intv(%u)", u, profile->chain.squash_intv);

	/* instanciate dz */
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
	size_t const plen = mm_strlen(pname);

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
	size_t const plen = mm_strlen(pname);	

	/* FIXME: save profile name */
	_unused(template);
	_unused(table);
	_unused(plen);

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
	tm_idx_profile_t const *profile = kv_ptr(self->score.profile)[pid];

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
tm_idx_sketch_t *tm_idx_save_seq(tm_ref_sketch_t *sr, uint32_t pid, char const *name, size_t nlen, uint8_t const *seq, size_t slen)
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
	for(size_t i = 0; i < slen; i += 16) {
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
		.rid  = 0,
		.pid  = pid,
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

	for(size_t i = 0; i < bseq_meta_cnt(&batch->bin); i++) {
		bseq_meta_t *p = &meta[i];

		p->u.ptr = NULL;		/* clear first */
		if(bseq_seq_len(p) > INT16_MAX) {
			debugblock({ fprintf(stderr, "i(%zu), len(%zu), seq(%s)\n", batch->bin.base_id + i, bseq_seq_len(p), bseq_name(p)); });
			continue;
		}

		/*
		debug("name(%s)", bseq_name(p));
		for(size_t j = 0; j < bseq_seq_len(p); j++) {
			fprintf(stderr, "%c", "NACMGRSVTWYHKDBN"[bseq_seq(p)[j]]);
		}
		fprintf(stderr, "\n");
		*/

		/* get score matrix for this sequence from its name with regex matching */
		size_t pid = tm_idx_find_profile(self,
			bseq_name(p),    bseq_name_len(p),
			bseq_comment(p), bseq_comment_len(p)
		);
		// debug("pid(%zu)", pid);

		/* build hash table */
		tm_ref_sketch_t *sr = tm_idx_build_sketch(self,
			ref, pid,
			bseq_name(p), bseq_name_len(p),
			bseq_seq(p),  bseq_seq_len(p)
		);
		debug("done, sr(%p)", sr);
		if(sr == NULL) {
			debugblock({ fprintf(stderr, "i(%zu), len(%zu), seq(%s)\n", batch->bin.base_id + i, bseq_seq_len(p), bseq_name(p)); });
			continue;
		}

		/* copy sequence */
		tm_idx_sketch_t *si = tm_idx_save_seq(
			sr, pid,	/* save profile id */
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
		bseq_meta_t const *p = &meta[i];
		if(_unlikely(p->u.ptr == NULL && bseq_seq_len(p))) {
			error("sequence too long (%.*s: length = %zu). removed.", (int)bseq_name_len(p), bseq_name(p), bseq_seq_len(p));
		}
		arr[base_rid + i] = p->u.ptr;		/* just copy */
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
		(pt_drain_t)tm_idx_record		/* no need for sorting */
	);
	for(size_t i = 0; i < pt_nth(pt); i++) {
		tm_ref_destroy_tbuf_static(&mii->ref[i]);
	}

	/* check error */
	bseq_close_t c = bseq_close(mii->col.fp);
	if(c.status != 0) {
		warn("broken file format detected for `%s'.", fn);
	}

	/* sanity check */
	// debug("read(%zu, %zu)", c.cnt, kv_cnt(mii->col.bin));
	return(c.cnt != kv_cnt(mii->col.bin));
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

size_t tm_idx_squash_invalid(tm_idx_sketch_t **sk, size_t scnt)
{
	size_t cnt = 0;
	for(size_t i = 0; i < scnt; i++) {
		if(sk[i] == NULL) { continue; }

		/* copy and assign rid */
		sk[cnt] = sk[i];
		sk[cnt]->h.rid = cnt;
		cnt++;
	}
	return(cnt);
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
	mi->sketch.cnt = tm_idx_squash_invalid(		/* remove NULL elements */
		kv_ptr(mii->col.bin),
		kv_cnt(mii->col.bin)
	);

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

static _force_inline
size_t tm_idx_string_size(char const *s)
{
	return(mm_strlen(s) + 1);
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
void *tm_idx_dump_string(tm_idx_dump_t *w, char const *s)
{
	char c = '\0';
	char const *p = s == 0 ? &c : s;
	size_t const l = tm_idx_string_size(s);
	return(tm_idx_dump_block(w, p, l));
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
	debug("pcnt(%zu)", mi->profile.cnt);
	for(size_t i = 0; i < mi->profile.cnt; i++) {
		buf[i] = tm_idx_dump_profile_core(w, buf[i]);
		debug("i(%zu), ofs(%zu)", i, (size_t)buf[i]);
		/* we don't add padding between profiles */
	}
	tm_idx_dump_pad(w);

	/* dump offsets (pointers) */
	void *ofs = tm_idx_dump_block(w, buf, sizeof(tm_idx_profile_t *) * mi->profile.cnt);
	tm_idx_dump_pad(w);

	/* done */
	free(buf);
	return(ofs);
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
	void *ofs = tm_idx_dump_block(w, buf, sizeof(tm_idx_sketch_t *) * mi->sketch.cnt);
	tm_idx_dump_pad(w);

	/* done */
	free(buf);
	return(ofs);
}

static _force_inline
size_t tm_idx_calc_size(tm_idx_t const *mi)
{
	size_t size = 0;

	size += tm_idx_scan_profile(mi);
	size += tm_idx_scan_sketch(mi);

	size += tm_idx_string_size(mi->filename.ref);
	size += tm_idx_string_size(mi->filename.profile);
	size += tm_idx_padding_size(size);		/* make aligned */

	size += sizeof(tm_idx_t);
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
		.size  = tm_idx_calc_size(&cmi)
	};
	wfp(fp, &hdr, sizeof(tm_idx_hdr_t));

	/* init dump context */
	tm_idx_dump_t w = {
		.fp = fp,
		.wfp = wfp,
		.ofs = 0
	};
	cmi.profile.arr = tm_idx_dump_profile(&w, &cmi);
	cmi.sketch.arr = tm_idx_dump_sketch(&w, &cmi);

	/* input seqence names */
	cmi.filename.ref     = tm_idx_dump_string(&w, cmi.filename.ref);
	cmi.filename.profile = tm_idx_dump_string(&w, cmi.filename.profile);

	/* done */
	tm_idx_dump_pad(&w);
	tm_idx_dump_block(&w, &cmi, sizeof(tm_idx_t));

	debug("size(%zu, %zu)", w.ofs, hdr.size);
	return(w.ofs);
}


/* laod */
static _force_inline
tm_idx_t *tm_idx_slice_root(uint8_t *base, size_t size)
{
	debug("size(%zu)", size);
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



/* scan-and-mask */

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

/* seed; see comment below for the definition */
typedef struct {
	uint32_t v, u;
} tm_seed_t;
_static_assert(sizeof(tm_seed_t) == 8);
typedef struct { tm_seed_t *a; size_t n, m; } tm_seed_v;

#define tm_seed_upos(x)			( (x).u )
KRADIX_SORT_INIT(seed, tm_seed_t, tm_seed_upos, 4);

#if 0
static _force_inline
tm_pair_t tm_seed_decode(tm_seed_t const *p)
{
	/* remove chained flag */
	uint32_t const u = p->u, v = p->v;

	/*
	 * convert to r-q corrdinates;
	 * u = 2 * q - r + 0x10000 + d
	 * v = 2 * r - q - 0x20000 - 2 * d
	 *     where d = 0 for forward and 0x40000000 for reverse
	 *
	 * 2 * u + v = 3 * q
	 * 2 * v + u = 3 * (r - 0x10000 - d)
	 */
	uint32_t const rt = (2 * v + u) & 0x7fffffff, r = (rt + 0x80030000) / 3U;
	uint32_t const q  = (2 * u + v) / 3U;

	uint32_t const vv = v & 0x3fffffff;
	uint32_t const uu = u & 0x3fffffff;
	uint32_t const rt2 = 2 * vv + uu, rr = rt2 + 0x30000, r2 = (rr - 0xc0000000) / 3U;
	uint32_t const qt2 = 2 * uu + vv, qq = qt2 & 0x3fffffff, q2 = qq / 3U;
	debug("r(%x), q(%x), v(%x), u(%x), r2(%x, %x, %x), q2(%x, %x, %x)", r, q, vv, uu, rt2, rr, r2, qt2, qq, q2);

	/* subtract offset */
	return((tm_pair_t){
		.r = r2,		/* leaves direction flag at bit 30 */
		.q = q2
	});
}
#endif


#if 0
static _force_inline
tm_chain_t tm_chain_decode(tm_chain_t const *p)
{
	return((tm_chain_t){ 0 });
}
#endif

static _force_inline
uint64_t tm_seed_dir(tm_seed_t const *s)
{
	return(1 - (s->v>>31));
}

static char *tm_seed_to_str(tm_seed_t const *s)
{
	// tm_pair_t p = tm_seed_decode(s);
	uint32_t const u = s->u & 0x7fffffff, v = s->v;		/* clear chained flag */
	uint32_t const rt = (2 * v + u), r = (rt + 0xc0030000) / 3U;
	uint32_t const qt = (2 * u + v), q = (qt & 0x7fffffff) / 3U;
	uint32_t const dir = (u & 0x80000000) != 0;
	return(xbprintf("seed(%x, %x), t(%x, %x), pos(%x, %x, %x)", s->v, s->u, qt, rt, dir, q, r));
}


typedef struct {
	uint16_t dst, src, cnt, unused;	/* upper 4bits must be cleared before use (overlaps with plen) */
} tm_sqiv_t;
_static_assert(sizeof(tm_sqiv_t) == 8);
typedef struct { tm_sqiv_t *a; size_t n, m; } tm_sqiv_v;



/* chain */

#if 0
typedef struct {
	uint16_t v, u;
} tm_span_t;

static _force_inline
tm_pair_t tm_span_decode(tm_span_t const *s)
{
	uint32_t const u = s->u, v = s->v;
	uint32_t const r = (2 * v + u) / 3U;
	uint32_t const q = (2 * u + v) / 3U;

	// debug("v(%x), u(%u), r(%u), q(%u)", v, u, r, q);
	return((tm_pair_t){
		.r = r,
		.q = q
	});
}
#endif

typedef struct {
	/*
	 * base seed position;
	 * rpos comes first because vpos comes first in seed_t.
	 * vpos is placed lower of (upos, vpos) tuple to make chaining simpler.
	 */
	uint16_t rpos;
	uint8_t dir, unused;
	uint32_t qpos;

	/* reference and query side spans */
	uint32_t rspan, qspan;
} tm_chain_raw_t;

static _force_inline
tm_pair_t tm_chain_raw_pos(tm_chain_raw_t const *chain)
{
	return((tm_pair_t){
		.q = chain->qpos,
		.r = chain->rpos		/* upper 16bits are masked out */
	});
}

static _force_inline
tm_pair_t tm_chain_raw_span(tm_chain_raw_t const *chain)
{
	return((tm_pair_t){
		.q = chain->qspan,
		.r = chain->rspan
	});
}


typedef struct {
	/* see tm_chain_compose */
	uint16_t rpos;					/* rspos + rspan for reverse */
	uint8_t dir, unused;
	uint32_t qpos;

	union {
		struct {
			uint32_t rid : 24;		/* reference id */
			uint32_t weight : 8;	/* weight; calculated from seed count and sequence length */
		} sep;
		uint32_t all;
	} attr;
	uint32_t qspan;
} tm_chain_t;
typedef struct { tm_chain_t *a; size_t n, m; } tm_chain_v;

_static_assert(sizeof(tm_chain_t) == 16);
_static_assert(sizeof(tm_chain_t) == sizeof(tm_chain_raw_t));

#define tm_chain_attr(x)		( (x).attr.all )
KRADIX_SORT_INIT(chain, tm_chain_t, tm_chain_attr, 4);


static char *tm_chain_raw_to_str(tm_chain_raw_t const *c)
{
	uint32_t const rspos = c->dir ? c->rpos + c->rspan : c->rpos;
	uint32_t const repos = c->dir ? c->rpos : c->rpos + c->rspan;
	return(xbprintf("dir(%u), (%u, %u) -- (%u, %u) --> (%u, %u)", c->dir, c->qpos, rspos, c->qspan, c->rspan, c->qpos + c->qspan, repos));
}
static char *tm_chain_to_str(tm_chain_t const *c)
{
	return(xbprintf("dir(%u), pos(%u, %u), qspan(%u), rid(%u), weight(%u)", c->dir & 0x01, c->qpos, c->rpos, c->qspan, c->attr.sep.rid, c->attr.sep.weight));
}

static _force_inline
tm_pair_t tm_chain_pos(tm_chain_t const *chain)
{
	debug("(%x, %x)", chain->qpos, chain->rpos);

	return((tm_pair_t){
		.q = chain->qpos,
		.r = chain->rpos		/* upper 16bits are masked out */
	});
}


/* alignment result */
typedef struct {
	/* rbt */
	rbt_header_t h;
	uint32_t qmax;

	/* stats */
	int32_t score;
	dz_alignment_t const *aln;

	/* positions */
	struct {
		uint32_t rid : 24;
		uint32_t max_weight : 8;
	} attr;
	uint8_t dir, unused[3];		/* 1 if reverse */
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
#define tm_aln_ivt_cmp_head(a, b)	({ \
	debug("compare head: a(%u, %u, %u), b(%u, %u, %u), cmp(%u > %u)", (a)->pos.q, (a)->pos.q + (a)->span.q, (a)->qmax, (b)->pos.q, (b)->pos.q + (b)->span.q, (b)->qmax, (a)->qmax, (b)->pos.q); \
	(a)->qmax                > (b)->pos.q; \
})
#define tm_aln_ivt_cmp_tail(a, b)	({ \
	debug("compare tail: a(%u, %u, %u), b(%u, %u, %u), cmp(%u < %u)", (a)->pos.q, (a)->pos.q + (a)->span.q, (a)->qmax, (b)->pos.q, (b)->pos.q + (b)->span.q, (b)->qmax, (a)->pos.q, (b)->pos.q + (b)->span.q); \
	(a)->pos.q               < (b)->pos.q + (b)->span.q; \
})
#define tm_aln_ivt_cmp_iter(a, b)	({ \
	debug("compare iter: a(%u, %u, %u), b(%u, %u, %u), cmp(%u > %u)", (a)->pos.q, (a)->pos.q + (a)->span.q, (a)->qmax, (b)->pos.q, (b)->pos.q + (b)->span.q, (b)->qmax, (a)->pos.q + (a)->span.q, (b)->pos.q); \
	(a)->pos.q + (a)->span.q > (b)->pos.q; \
})

static _force_inline
uint64_t tm_aln_ivt_update(tm_aln_t *parent, tm_aln_t *child)
{
	if(child == NULL) { return(1); }

	uint32_t const cepos = child->pos.q + child->span.q;
	uint32_t const qmax = parent->qmax;

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

/* array of alignment; allocated in trace stack */
typedef struct {
	size_t cnt;
	tm_aln_t arr[];
} tm_alnv_t;


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
		uint32_t max_weight;	/* working variable */
		dz_arena_t *fill, *trace;

		/* result bin */
		tm_aln_v arr;
	} extend;
} tm_scan_t;


static _force_inline
void tm_scan_init_static(tm_scan_t *self)
{
	memset(self, 0, sizeof(tm_scan_t));

	size_t const size = DZ_MEM_INIT_SIZE - _roundup(sizeof(dz_arena_t), DZ_MEM_ALIGN_SIZE);
	self->extend.fill = dz_arena_init(size);
	self->extend.trace = NULL;
	return;
}

static _force_inline
void tm_scan_destroy_static(tm_scan_t *self)
{
	_unused(self);
	kv_destroy(self->seed.arr);
	kv_destroy(self->seed.sqiv);
	kv_destroy(self->chain.arr);
	kv_destroy(self->extend.arr);

	dz_arena_destroy(self->extend.fill);
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

	size_t const scnt = kv_cnt(*buf);
	tm_seed_t *p = kv_reserve(tm_seed_t, *buf, scnt + 65536);

	// debug("scnt(%zu), q(%p), p(%p, %p)", scnt, q, p, &p[scnt]);
	return(&p[scnt]);
}

static _force_inline
size_t tm_expand_seed(tm_ref_state_t s, v4i32_t uofs, v4i32_t vofs, tm_seed_t *q)
{
	uint32_t const mask = 0xc000ffff;
	tm_ref_match_t m = tm_ref_get_arr(s);

	#if 1
	v4i32_t const wmask = _set_v4i32(mask);		/* leave lower 15bits and direction flag */
	for(size_t i = 0; i < m.cnt; i += 4) {
		/*
		 * u = 2 * r  - q'
 		 * v = 2 * q' - r
 		 */
		v4i32_t const z = ({
			__m128i const x = _mm_loadl_epi64((__m128i const *)&m.ptr[i]);	/* movq (mem), xmm */
			__m128i const y = _mm_cvtepi16_epi32(x);	/* pmovsxwd; sign expansion (negative for reverse) */
			((v4i32_t){ .v1 = y });
		});
		v4i32_t const w = _and_v4i32(z, wmask);			/* disjoin forward and reverse of the reference sequence */

		v4i32_t const u = _sub_v4i32(uofs, w);					/* 2 * q - r + 0x10000 */
		v4i32_t const v = _add_v4i32(vofs, _add_v4i32(w, w));	/* 2 * r - q - 0x20000 */

		/* save: u in upper and v in lower */
		_storeu_v4i32(&q[i],     _lo_v4i32(u, v));
		_storeu_v4i32(&q[i + 2], _hi_v4i32(u, v));
	}
	#else
	for(size_t i = 0; i < m.cnt; i++) {
		int16_t const z = m.ptr[i];
		uint32_t const w = ((int32_t)z) & mask;			/* signed expansion */
		uint32_t const u = _ext_v4i32(uofs, 0) - w;
		uint32_t const v = _ext_v4i32(vofs, 0) + 2 * w;
		q[i].v = v;
		q[i].u = u;

		_scan_memory(q, sizeof(tm_seed_t));
	}
	#endif

/*
	fprintf(stderr, "bin(%p), cnt(%zu)\n", s.bin, m.cnt);
	for(size_t i = 0; i < m.cnt; i++) {
		debug("qr(%x, %x), vu(%x, %x)", (0x80000000 - 0x20000) - _ext_v4i32(vofs, 0), m.ptr[i], q[i].v, q[i].u);
	}
*/
	return(m.cnt);
}


static _force_inline
tm_sqiv_t *tm_sqiv_reserve(tm_sqiv_v *buf, tm_sqiv_t *q)
{
	/* update count */
	kv_cnt(*buf) = q - kv_ptr(*buf);

	size_t const vcnt = kv_cnt(*buf);
	tm_sqiv_t *p = kv_reserve(tm_sqiv_t, *buf, vcnt + 64);
	return(&p[vcnt]);
}

static _force_inline
size_t tm_save_sqiv(tm_ref_squash_t sq, tm_sqiv_t *r)
{
	/* just copy */
	_storeu_v2i32(r, sq.v);
	return(1);
}

#if 0
static _force_inline
void tm_invalidate_sqiv(tm_sqiv_t *p)
{
	_storeu_u32(&p->dst, 0);		/* p->dst = p->src = 0 */
	return;
}
#endif

static _force_inline
size_t tm_save_last_sqiv(tm_ref_state_t s, tm_sqiv_t *r)
{
	tm_ref_squash_t sq = tm_ref_load_squash(s.bin, -1LL, 0);	/* unmatching = -1 to insert zeros to lower 32bit */
	_storeu_v2i32(r, sq.v);
	return(1 + s.unmatching);
}

static _force_inline
size_t tm_collect_seed(tm_ref_sketch_t const *ref, uint8_t const *query, size_t qlen, size_t intv, tm_seed_v *seed, tm_sqiv_v *sqiv)
{
	/* load coordinate constants */
	v4i32_t const uinc = _set_v4i32(2), vinc = _set_v4i32(-1);
	v4i32_t uofs = _set_v4i32(0x10000 + 2 * (TM_SEED_POS_OFS + 1));
	v4i32_t vofs = _set_v4i32(0x80000000 - 0x20000 - (TM_SEED_POS_OFS + 1));	/* query offsets */

	tm_seed_t *q = tm_seed_reserve(seed, kv_ptr(*seed));
	tm_sqiv_t *r = tm_sqiv_reserve(sqiv, kv_ptr(*sqiv));

	size_t const mask = intv - 1;

	/* initial state (of three general-purpose registers) */
	tm_ref_state_t s = tm_ref_match_init(ref);

	/* for each base... */
	uint8_t const *t = &query[qlen];
	size_t rem = qlen + 1;
	while(--rem != 0) {
		/* update coordinates */
		uofs = _add_v4i32(uofs, uinc);
		vofs = _add_v4i32(vofs, vinc);

		/* update matching status for the next bin (prefetch) */
		tm_ref_next_t n = tm_ref_match_next(s, (rem & mask) == 0, t[-rem]);

		/* save matching state of the previous bin */
		r += tm_save_sqiv(n.squash, r) + s.unmatching;		/* do not forward pointer if previous bin did not match */

		/* skip current bin if not matching */
		s = n.state;
		if(s.unmatching == 0) {
			q += tm_expand_seed(s, uofs, vofs, q);
		}

		/* test array expansion every 32 times */
		if((rem & 0x1f) != 0) { continue; }
		q = tm_seed_reserve(seed, q);
		r = tm_sqiv_reserve(sqiv, r);
	}

	/* save the last sqiv */
	r += tm_save_last_sqiv(s, r);

	/* update seed count */
	size_t const scnt = q - kv_ptr(*seed);
	size_t const vcnt = r - kv_ptr(*sqiv);
	kv_cnt(*seed) = scnt;
	kv_cnt(*sqiv) = vcnt;
	return(scnt);
}

static _force_inline
size_t tm_squash_seed(tm_sqiv_t const *sqiv, size_t icnt, tm_seed_v *seed)
{
	tm_seed_t *dp = kv_ptr(*seed);
	tm_seed_t const *sp = dp;
	// debug("sqiv(%p), icnt(%zu)", sqiv, icnt);

	/* squash succeeding match */
	for(size_t i = 0; i < icnt; i++) {
		_scan_memory(&sqiv[i], sizeof(tm_sqiv_t));
		tm_sqiv_t sq = {
			.dst = sqiv[i].dst & 0xfff,
			.src = sqiv[i].src & 0xfff,
			.cnt = sqiv[i].cnt & 0xfff
		};

		// debug("spos(%zu), dpos(%zu), cnt(%u, %u, %u)", sp - kv_ptr(*seed), dp - kv_ptr(*seed), sq.dst, sq.src, sq.cnt);

		/* former half */
		for(size_t j = 0; j < sq.dst; j += 4) {
			v32i8_t const v = _loadu_v32i8(&sp[j]);
			_storeu_v32i8(&dp[j], v);
		}

		/* latter half */
		dp -= sq.src - sq.dst;
		for(size_t j = sq.src; j < sq.cnt; j += 4) {
			v32i8_t const v = _loadu_v32i8(&sp[j]);
			_storeu_v32i8(&dp[j], v);
		}

		/* forward pointers */
		sp += sq.cnt;
		dp += sq.cnt;
		// debug("(%u, %u, %u), idx(%zu, %zu), cnt(%zu), squashed(%zu)", sq.d, sq.s, sq.c, dp - kv_ptr(*seed), sp - kv_ptr(*seed), sq.c, sq.s - sq.d);
	}

	// debug("sp(%zu), scnt(%zu)", sp - kv_ptr(*seed), kv_cnt(*seed));
	debugblock({
		if((size_t)(sp - kv_ptr(*seed)) != kv_cnt(*seed)) { trap(); }
	});

	size_t const cnt = dp - kv_ptr(*seed);
	kv_cnt(*seed) = cnt;
	return(cnt);
}

static _force_inline
int64_t tm_chain_test_ptr(tm_seed_t const *p, tm_seed_t const *t)
{
	// debug("check tail, test(%ld)", (int64_t)((ptrdiff_t)(t - p - 1)));
	return((int64_t)((ptrdiff_t)(t - p - 1)));
}

static _force_inline
tm_seed_t *tm_chain_find_first(uint64_t ruv, tm_seed_t *p, tm_seed_t *t, uint64_t lb, uint64_t window)
{
	/* constants */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const tmask = 0x8000000080000000;	/* extract two sign bit pair */

	/* load bounds */
	uint64_t const uv = ruv;					/* (ulb, vlb) */
	uint64_t ub = (uv & umask) + window;		/* (uub, vub - vlb) */
	debug("uub(%lx), vwlen(%lx)", ub>>32, ub & 0xffffffff);

	int64_t cont = 0;		/* continuous flag */
	while(1) {
		/* check if reached tail */
		if((cont | tm_chain_test_ptr(++p, t)) < 0) { return(NULL); }

		uint64_t const v = _loadu_u64(p) - lb;	/* 1: (u, v - vlb) for inclusion test */
		debug("test inclusion, pv(%lx), pu(%lx), u(%lx), {v - vlb}(%lx), test(%lx)", p->v, p->u, v>>32, v & 0xffffffff, (ub - v) & tmask);

		cont = ub - v;		/* save diff for testing MSb */
		if(((cont | v) & tmask) == 0) {			/* 2,3: break if chainable (first chainable seed found) */
			debug("chainable");
			break;
		}

		/* unchainable */
		debug("unchainable");
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
	debug("u(%lx), {v - vlb}(%lx), uub(%lx), vub(%lx)", rv>>32, rv & 0xffffffff, ub>>32, ub & 0xffffffff);

	/* keep nearest seed */
	tm_seed_t *n = p;

	int64_t cont = 0;
	while((cont | tm_chain_test_ptr(++p, t)) >= 0) {
		uint64_t const v = _loadu_u64(p) - lb;	/* 1: (u, v - vlb) for inclusion test */
		debug("u(%lx), {v - vlb}(%lx), uub(%lx), vub(%lx), test(%lx), term(%lx)", p->u, v & 0xffffffff, ub>>32, ub & 0xffffffff, (ub - v) & tmask, (int64_t)(ub - v) < 0);
		cont = ub - v;
		if((cont | v) & tmask) {		/* skip if unchainable */
			debug("unchainable");
			continue;
		}

		/* chainable; test if the seed is nearer than the previous */
		uint64_t const w = v + (v<<32);			/* 2,3: (u + v - vlb, v - vlb) */
		debug("{u + v - vlb}(%lx), {v - vlb}(%lx)", w>>32, w & 0xffffffff);
		if((ub - w) & tmask) {			/* 4,5: further than previous */
			debug("further; terminate");
			break;
		}

		/* nearer seed found */
		ub = w;			/* update bounds */
		n = p;			/* save pointer */
	}
	return(n);
}

static _force_inline
tm_seed_t *tm_chain_find_nearest(uint64_t ruv, tm_seed_t *p, tm_seed_t *t, uint64_t window)
{
	/* load root positions */
	uint64_t const umask = 0xffffffff00000000;	/* extract upper */
	uint64_t const uv = ruv;					/* (ulb, vlb) */
	uint64_t lb = (uv & ~umask);				/* (0, vlb) */
	debug("ulb(%lx), vlb(%lx)", uv>>32, lb);

	/* find first chainable seed */
	tm_seed_t *n = tm_chain_find_first(ruv, p, t, lb, window);
	if(n == NULL) { return(NULL); }				/* we expect control paths are kept isolated */

	/* then inspect alternative chainable seeds */
	return(tm_chain_find_alt(n, t, lb));
}

#if 0
static _force_inline
size_t tm_chain_record(v2i32_t root, v2i32_t tail, v2i32_t cnt, v2i32_t min_cnt, tm_chain_t *q)
{
	/* avoid use of general-purpose register to minimize spills */
	v2i32_t const span = _sub_v2i32(tail, root);

	/* FIXME: raw SSE operations */
	v4i32_t const chain = ({
		__m128i const x = _mm_shufflelo_epi16(span.v1, 0xd8);	/* (uint32_t [4]){ -, -, -, (uspan, vspan) } */
		__m128i const y = _mm_unpacklo_epi32(x, cnt.v1);		/* (uint32_t [4]){ -, -, cnt, (uspan, vspan) } */
		__m128i const z = _mm_unpacklo_epi64(root.v1, y);		/* (uint32_t [4]){ cnt, (uspan, vspan), upos, vpos } */
		((v4i32_t){ .v1 = z });
	});
	_storeu_v4i32(q, chain);

	/* forward if cnt > min_cnt */
	v2i32_t const fv = _gt_v2i32(min_cnt, cnt);
	uint32_t const fwd = _ext_v2i32(fv, 0);
	return((uint32_t)(1 + fwd));
}

static _force_inline
void tm_chain_record(v2i32_t root, v2i32_t tail, v2i32_t cnt, tm_chain_t *q)
{
	/* avoid use of general-purpose register to minimize spills */
	v2i32_t const span = _sub_v2i32(tail, root);

	/* FIXME: raw SSE operations */
	v4i32_t const chain = ({
		__m128i const x = _mm_shufflelo_epi16(span.v1, 0xd8);	/* (uint32_t [4]){ -, -, -, (uspan, vspan) } */
		__m128i const y = _mm_unpacklo_epi32(x, cnt.v1);		/* (uint32_t [4]){ -, -, cnt, (uspan, vspan) } */
		__m128i const z = _mm_unpacklo_epi64(root.v1, y);		/* (uint32_t [4]){ cnt, (uspan, vspan), upos, vpos } */
		((v4i32_t){ .v1 = z });
	});
	_storeu_v4i32(q, chain);
	return;
}
#endif

#define _print_v4i32x(a) { \
	debug("(v4i32_t) %s(%x, %x, %x, %x)", #a, _ext_v4i32(a, 3), _ext_v4i32(a, 2), _ext_v4i32(a, 1), _ext_v4i32(a, 0)); \
}

/*
 * uv -> xy coordinate conversion: epos[2] must be zero
 */
static _force_inline
v4i32_t tm_chain_compose(v2i32_t spos, v2i32_t epos)
{
	v2i64_t const coef = _seta_v2i64(0x55555556, 0x55555556);	/* 0.3333... + roundup alpha */
	v4i32_t const ofs  = _seta_v4i32(
		0, 0x80000000 - 3 * TM_SEED_POS_OFS,	/* -3 * base_offset */
		0, 0x00030000 - 3 * TM_SEED_POS_OFS		/* 3 * direction flag - 3 * base_offset */
	);

	v4i32_t const a = _lo_v4i32(epos, spos);	/* 1; (ue, us, ve, vs) */
	v4i32_t const b = _add_v4i32(				/* 3; (ye, ys, xe, xs) = (2*ue + ve, 2*us + vs, 2*ve + ue, 2*vs + us) */
		_add_v4i32(a, a),						/* 2; (2*ue, 2*us, 2*ve, 2*vs) */
		_shuf_v4i32(a, _km_v4i32(1, 0, 3, 2))	/* 2; (ve, vs, ue, us) */
	);
	// _print_v4i32x(spos);
	// _print_v4i32x(epos);
	// _print_v4i32x(a);
	// _print_v4i32x(b);

	/* decode span */
	v4i32_t const d = (v4i32_t){
		_mm_srli_epi64(b.v1, 32)				/* 4; (0, ye, 0, xe) */
	};
	v4i32_t const c = _sub_v4i32(d, b);			/* 5; (-ye, ye - ys, -xe, xe - xs) */
	v4i32_t const span = (v4i32_t){
		_mm_mul_epu32(coef.v1, c.v1)			/* 6-; (qspan, -, rspan, -) */
	};
	// _print_v4i32x(d);
	// _print_v4i32x(c);
	_print_v4i32x(span);

	/* decode pos */
	v4i32_t const e = _bsr_v4i32(c, 1);			/* 6; (-, -ye, -, -xe) */
	v4i32_t const f = (v4i32_t){
		_mm_blendv_epi32(b.v1, e.v1, epos.v1)	/* 7; (-, ys, -, dir ? -xe : xs) */
	};
	v4i32_t const g = _add_v4i32(f, ofs);		/* 8; (-, ys + 0x80000000, -, (dir ? -xe : xs) + 0xc0030000) */
	v4i32_t const pos = (v4i32_t){
		_mm_mul_epu32(coef.v1, g.v1)			/* 9-; (qpos, 0, rpos, -) */
	};
	// _print_v4i32x(e);
	// _print_v4i32x(f);
	// _print_v4i32x(g);
	_print_v4i32x(pos);

	/* fold pos and span */
	v4i32_t const h = (v4i32_t){
		_mm_shuffle2_epi32(pos.v1, span.v1, 0xdd)		/* 15; (qspan, rspan, qpos, rpos) */
	};
	_print_v4i32x(h);
	return(h);		/* (qspan, rspan, qpos, rpos) */
}

static _force_inline
tm_chain_t *tm_chain_reserve(tm_chain_v *chain, size_t scnt)
{
	size_t const used = kv_cnt(*chain);
	tm_chain_t *base = kv_reserve(tm_chain_t, *chain, used + scnt);
	return(&base[used]);
}

static _force_inline
size_t tm_chain_finalize(tm_chain_v *chain, tm_chain_t *q)
{
	tm_chain_t *base = kv_ptr(*chain);
	size_t const used = kv_cnt(*chain);
	size_t const ccnt = q - &base[used];
	// debug("cnt(%zu, %zu), q(%p, %p)", ccnt, used, q, base);
	kv_cnt(*chain) = used + ccnt;
	return(ccnt);
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
	uint64_t const chained = 0x80000000;				/* offset */
	uint64_t const window  = profile->chain.window.all;	/* window sizes */

	v4i32_t const inc = _seta_v4i32(0, 1, 0, 0);
	v4i32_t const adj = _load_v4i32(&profile->chain.kadj[0]);	/* (-, min_scnt, k, k) */

	/* src pointers */
	tm_seed_t *p = seed - 1, *t = &seed[scnt];

	/* reserve mem for chain */
	tm_chain_t *q = tm_chain_reserve(chain, scnt);
	// debug("scnt(%zu)", scnt);

	/* for each seed */
	while(++p < t) {
		/* skip chained seed to find next root */
		// debug("u(%u), chained(%lu)", p->u, p->u >= chained);
		if(p->u >= chained) { continue; }

		/* root found; keep seed on xmm register */
		v2i32_t const root = _load_v2i32(p);	/* we expect this won't be spilled */
		v4i32_t spos = _sub_v4i32((v4i32_t){ root.v1 }, adj);	/* (-, -min_scnt, vspos - k, uspos - k) */

		/* iteratively link nearest seed */
		tm_seed_t *s = p;
		uint64_t ruv = _loadu_u64(s);
		while(1) {
			tm_seed_t *n = tm_chain_find_nearest(ruv, s, t, window);
			if(n == NULL) { break; }

			/* increment seed count */
			spos = _add_v4i32(spos, inc);

			s = n;					/* save last chained seed */
			ruv = _loadu_u64(s);	/* before flag set */
			s->u += chained;		/* mark chained */
		}

		debugblock({
			if((int32_t)_ext_v2i32(spos, 2) >= 0) {
				debug("%r ---> %r, cnt(%d)", tm_seed_to_str, p, tm_seed_to_str, s, _ext_v2i32(spos, 2));
			}
		});

		v2i32_t const epos = (v2i32_t){ .v1 = _mm_cvtsi64_si128(ruv) };
		if((int32_t)_ext_v2i32(spos, 2) < 0) { continue; }		/* skip if short */

		/* save */
		v4i32_t const c = tm_chain_compose((v2i32_t){ spos.v1 }, epos);
		_store_v4i32(q, c);
		q++;
	}

	/* update chain count */
	return(tm_chain_finalize(chain, q));
}



/* filter */

typedef struct {
	v16i8_t scv, gv, cv, pv;
	uint8_t fw[64], rv[64];
} tm_filter_work_t;


static _force_inline
void tm_filter_work_init(tm_filter_work_t *w, tm_idx_profile_t const *profile)
{
	w->scv = _load_v16i8(&profile->filter.score_matrix); 
	w->gv  = _load_v16i8(&profile->filter.gap);
	w->pv  = _load_v16i8(&profile->filter.init);

	/* derive p = 1 vector */
	w->cv  = _add_v16i8(w->pv, w->gv);
	w->cv  = _maxu_v16i8(w->cv, _bsr_v16i8(w->cv, 1));
	return;
}


/* 4-bit encoding for reference side; return in reversed orientation */
static _force_inline
v32i8_t tm_filter_load_rf(uint8_t const *p) {
	v32i8_t const v = _loadu_v32i8(p);
	return(_swap_v32i8(v));
}

static _force_inline
v32i8_t tm_filter_load_rr(uint8_t const *p) {
	v32i8_t const v = _loadu_v32i8(p - 32);
	return(v);
}

/* 2-bit encoding at [3:2] for query side */
static _force_inline
v32i8_t tm_filter_load_qf(uint8_t const *p, uint64_t dir) {
	/* convert to 4bit */
	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		[nA]     = A, [nC]     = C, [nG]     = G, [nT]     = T,	/* forward */
		[nA + 1] = T, [nC + 1] = G, [nG + 1] = C, [nT + 1] = A	/* reverse-complemented */
	};

	v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(conv));
	v32i8_t const v = _loadu_v32i8(p), d = _set_v32i8(dir);
	v32i8_t const w = _shuf_v32i8(cv, _add_v32i8(v, d));
	return(w);
}

static _force_inline
v32i8_t tm_filter_load_qr(uint8_t const *p, uint64_t dir) {
	/* reverse and complement */
	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		[nA]     = A, [nC]     = C, [nG]     = G, [nT]     = T,	/* forward */
		[nA + 1] = T, [nC + 1] = G, [nG + 1] = C, [nT + 1] = A	/* reverse-complemented */
	};

	v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(conv));
	v32i8_t const v = _loadu_v32i8(p - 32), d = _set_v32i8(dir);
	v32i8_t const w = _shuf_v32i8(cv, _add_v32i8(v, d));
	return(_swap_v32i8(w));
}

static _force_inline
void tm_filter_load_seq(tm_filter_work_t *w, uint8_t const *r, uint8_t const *q, tm_chain_raw_t const *c)
{
	/* load positions */
	uint64_t const dir = c->dir;
	tm_pair_t const pos  = tm_chain_raw_pos(c);
	tm_pair_t const span = tm_chain_raw_span(c);
	// debug("rr(%u), rf(%u), qr(%u), qf(%u), dir(%x)", pos.r, pos.r + span.r, pos.q, pos.q + span.q, dir);

	/* ref */
	v32i8_t const rf = tm_filter_load_rf(&r[pos.r + span.r]);
	v32i8_t const rr = tm_filter_load_rr(&r[pos.r]);
	_store_v32i8(&w->fw[0], rf);
	_store_v32i8(&w->rv[0], rr);

	/* query */
	v32i8_t const qf = tm_filter_load_qf(&q[pos.q + span.q], dir);
	v32i8_t const qr = tm_filter_load_qr(&q[pos.q], dir);

	uint8_t *qfp = (dir & 0x01) ? &w->rv[32] : &w->fw[32];
	uint8_t *qrp = (dir & 0x01) ? &w->fw[32] : &w->rv[32];
	_store_v32i8(qfp, qf);
	_store_v32i8(qrp, qr);

	// _print_v32i8(rf);
	// _print_v32i8(rr);
	// _print_v32i8(qf);
	// _print_v32i8(qr);
	return;
}

static _force_inline
int64_t tm_filter_extend_core(tm_filter_work_t *w, uint8_t const *buf)
{
	uint8_t const *r = &buf[32 - 8], *q = &buf[32 - 8];
	v16i8_t const scv = w->scv, gv = w->gv;
	v16i8_t pv = w->pv, cv = w->cv, mv = w->pv;

	#define _calc_next(_rv, _qv, _pv, _cv, _shift) ({ \
		/* calc score profile vector with shuffling score matrix vector */ \
		v16i8_t const match = _and_v16i8(_rv, _qv);		/* cmpeq */ \
		v16i8_t const sv = _shuf_v16i8(scv, match); \
		v16i8_t const x = _add_v16i8(_pv, sv); \
		v16i8_t const y = _add_v16i8(_cv, gv); \
		v16i8_t const z = _maxu_v16i8(y, _shift(y, 1)); \
		v16i8_t const n = _maxu_v16i8(x, z); \
		n; \
	})

	/* 32 vector updates */
	v16i8_t qv = _loadu_v16i8(q);		/* load query side initial vector */
	for(size_t i = 0; i < 16; i++) {
		/* #0 (even) */
		r--;
		v16i8_t const rv = _loadu_v16i8(r);	/* move rightward */
		v16i8_t const ev = _calc_next(rv, qv, pv, cv, _bsl_v16i8);
		mv = _maxu_v16i8(mv, ev);
		// _printu_v16i8(ev);
		// _printu_v16i8(mv);

		/* #1 (odd) */
		q++;
		qv = _loadu_v16i8(q);	/* move downward */
		v16i8_t const ov = _calc_next(rv, qv, cv, ev, _bsr_v16i8);
		mv = _maxu_v16i8(mv, ov);
		// _printu_v16i8(ov);
		// _printu_v16i8(mv);

		/* alias registers for the next loop */
		pv = ev;
		cv = ov;
	}

	/* fold scores */
	return(_hmaxu_v16i8(mv) - 128LL);

	#undef _calc_next
}

static _force_inline
int64_t tm_filter_extend(tm_filter_work_t *w, tm_idx_profile_t const *profile, uint8_t const *r, uint8_t const *q, tm_chain_raw_t const *c)
{
	/* skip if long enough */
	if(c->qspan >= 128) { return(1); }

	/* load sequences */
	tm_filter_load_seq(w, r, q, c);

	/* extend */
	int64_t const fw = tm_filter_extend_core(w, w->fw);
	int64_t const rv = tm_filter_extend_core(w, w->rv);
	debug("fw(%ld), rv(%ld), score(%ld)", fw, rv, fw + rv);
	return(fw + rv >= profile->filter.min_score);
}

#if 0
static _force_inline
size_t tm_filter_save_chain(uint32_t rid, tm_chain_t *q, tm_chain_raw_t const *p)
{
	/* load */
	tm_chain_raw_t const raw = *p;
	debug("p(%p), %r", p, tm_chain_raw_to_str, p);

	/* calc weight */
	ZCNT_RESULT size_t weight = _lzcnt_u32(p->rspan);

	/* save */
	*q = (tm_chain_t){
		.rpos = raw.rpos + (raw.dir ? raw.rspan : 0),
		.qpos = raw.qpos,
		.qspan = raw.qspan,

		.dir  = raw.dir,
		.unused = raw.unused,

		.attr.sep.rid = rid,			/* overwrite */
		.attr.sep.weight = weight	/* negated and offsetted by 32; the larger for the longer chain */
	};
	debug("q(%p), %r, rid(%u), weight(%u)", q, tm_chain_to_str, q, rid, weight);
	return(1);
}
#else
static _force_inline
size_t tm_filter_save_chain(uint32_t rid, tm_chain_t *q, tm_chain_raw_t const *p)
{
	/* load */
	v2i32_t const v = _loadu_v2i32(p);
	uint8_t const dir = p->dir;
	uint32_t const rspan = p->rspan;
	debug("p(%p), %r", p, tm_chain_raw_to_str, p);

	/* calc weight */
	ZCNT_RESULT size_t weight = _lzcnt_u32(p->rspan);

	/* save */
	_storeu_v2i32(q, v);
	if(dir) { q->rpos += rspan; }

	q->attr.sep.rid = rid;			/* overwrite */
	q->attr.sep.weight = weight;	/* negated and offsetted by 32; the larger for the longer chain */
	debug("q(%p), %r, rid(%u), weight(%u)", q, tm_chain_to_str, q, rid, weight);
	return(1);
}
#endif

static _force_inline
size_t tm_filter_chain(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, tm_chain_t *chain, size_t ccnt)
{
	/* alias pointer */
	tm_chain_raw_t const *src = (tm_chain_raw_t const *)chain;
	tm_chain_t *dst = chain;

	uint32_t const rid = tm_idx_ref_rid(si);		/* FIXME: kept on xmm? */
	uint8_t const *ref = tm_idx_ref_seq_ptr(si);

	/* init working buffer (load score vectors on xmm) */
	tm_filter_work_t w __attribute__(( aligned(32) ));
	tm_filter_work_init(&w, profile);

	debug("test_cnt(%u), min_score(%d), uspan_thresh(%u), rid(%u)", profile->filter.test_cnt, profile->filter.min_score, profile->filter.uspan_thresh, rid);
	for(size_t i = 0; i < ccnt; i++) {
		tm_chain_raw_t const *p = &src[i];
		debug("i(%zu), ccnt(%zu), %r", i, ccnt, tm_chain_raw_to_str, p);

		/* try short extension if not heavy enough */
		if(!tm_filter_extend(&w, profile, ref, query, p)) { continue; }

		/* copy and fold in reference id */
		dst += tm_filter_save_chain(rid, dst, p);
	}
	debug("cfilt(%zu)", dst - chain);
	return(dst - chain);
}

static _force_inline
size_t tm_seed_and_sort(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, size_t qlen, tm_seed_v *seed, tm_sqiv_v *sqiv)
{
	/* enumerate seeds */
	kv_clear(*seed);
	kv_clear(*sqiv);
	if(tm_collect_seed(&si->ref, query, qlen, profile->chain.squash_intv, seed, sqiv) == 0) {
		return(0);
	}

	/* squash overlapping seeds */
	if(tm_squash_seed(kv_ptr(*sqiv), kv_cnt(*sqiv), seed) == 0) {
		return(0);
	}

	/* sort seeds by u-coordinates for chaining */
	tm_seed_t *sptr = kv_ptr(*seed);
	size_t const scnt = kv_cnt(*seed);

	radix_sort_seed(sptr, scnt);
/*
	debug("scnt(%zu)", scnt);
	for(size_t i = 0; i < scnt; i++) {
		debug("i(%zu), %r", i, tm_seed_to_str, &sptr[i]);
	}
*/
	return(scnt);
}

static _force_inline
size_t tm_chain_and_filter(tm_idx_sketch_t const *si, tm_idx_profile_t const *profile, uint8_t const *query, tm_seed_t *seed, size_t scnt, tm_chain_v *chain)
{
	/* chaining; return if no chain obtained */
	size_t const cbase = kv_cnt(*chain);		/* save chain count */
	size_t const ccnt  = tm_chain_seed(profile, seed, scnt, chain);
	if(ccnt == 0) { return(0); }
	debug("cbase(%zu), ccnt(%zu)", cbase, ccnt);

	/* filter (try small extension with simple score matrix) */
	size_t const cfilt = tm_filter_chain(si, profile, query, &kv_ptr(*chain)[cbase], ccnt);
	
	/* update chain count */
	kv_cnt(*chain) = cbase + cfilt;
	return(cfilt);
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
		if(tm_seed_and_sort(si[i], pf[si[i]->h.pid], seq, slen, &self->seed.arr, &self->seed.sqiv) == 0) { continue; }

		/* chain and filter */
		tm_seed_t *seed = kv_ptr(self->seed.arr);
		size_t scnt = kv_cnt(self->seed.arr);
		tm_chain_and_filter(si[i], pf[si[i]->h.pid], seq, seed, scnt, &self->chain.arr);

		// fprintf(stderr, "sketch(%zu)\n", i);
	}

	debug("ccnt(%zu)", kv_cnt(self->chain.arr));
	return(kv_cnt(self->chain.arr));
}



/* X-drop DP extension */
#define TM_WEIGHT_MARGIN			( 3U )

static _force_inline
tm_aln_t tm_chain_as_aln(tm_chain_t const *c)
{
	return((tm_aln_t){
		.pos  = tm_chain_pos(c),
		.span = { .q = c->qspan, .r = c->qspan },	/* duplicate qspan because rspan is missing in chain object */
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
uint64_t tm_extend_test_range(tm_chain_t const *c, tm_aln_t const *p)
{
	/*
	 * test q-range; skip if entirely covered
	 * c->qpos > p->qpos && c->qpos + c->qspan < p->qpos + p->qspan
	 */
	tm_pair_t const pos = tm_chain_pos(c);
	if((uint32_t)(pos.q - p->pos.q) < (uint32_t)(p->span.q - c->qspan)) {
		debug("covered, (%u, %u), (%u, %u)", pos.q, pos.q + c->qspan, p->pos.q, p->pos.q + p->span.q);

		if(c->attr.sep.weight + TM_WEIGHT_MARGIN > p->attr.max_weight + 1U) {
			debug("weight not enough, weight(%u), max_weight(%u)", c->attr.sep.weight, p->attr.max_weight);
			return(1);	/* still too shorter than overlapping alignment */
		}
	}
	return(0);
}

static _force_inline
uint64_t tm_extend_is_covered(tm_scan_t *self, tm_chain_t const *c)
{
	debug("%r", tm_chain_to_str, c);

	/* for each overlapping alignment */
	tm_aln_t const *v = kv_ptr(self->extend.arr);
	if(v == NULL) { return(0); }	/* no result found */

	/* init iterator */
	rbt_iter_t it;
	tm_aln_t const caln = tm_chain_as_aln(c);		/* convert chain position to (pseudo) alignment range */
	rbt_init_iter_aln(&it, v, &caln);

	tm_aln_t const *p = rbt_fetch_head_aln(&it, v, &caln);
	while(p != NULL) {
		debug("p(%p), weight(%u), (%u) --- %u --> (%u)", p, p->attr.max_weight, p->pos.q, p->span.q, p->pos.q + p->span.q);

		/* overlapping and heavy enough */
		if(c->attr.sep.weight > p->attr.max_weight) {	/* smaller weight is better */
			debug("weight too large, weight(%u), max_weight(%u)", c->attr.sep.weight, p->attr.max_weight);
			return(1);		/* already covered */
		}

		if(tm_extend_test_range(c, p)) {
			return(1);		/* still too shorter than overlapping alignment */
		}

		p = rbt_fetch_next_aln(&it, v, &caln);
	}
	return(0);				/* possiblilty remains for better alignment */
}


static _force_inline
tm_pair_t tm_aln_span(dz_alignment_t const *aln)
{
	return((tm_pair_t){
		.r = aln->ref_length,
		.q = aln->query_length
	});
}

static _force_inline
void tm_extend_record(tm_scan_t *self, tm_chain_t const *chain, tm_pair_t spos, dz_alignment_t const *aln)
{
	_unused(chain);

	tm_pair_t const span = tm_aln_span(aln);
	tm_aln_t a = {
		/* coordinates */
		.qmax = spos.q + span.q,	/* spos.q + span.q, */
		.dir  = chain->dir,
		.pos  = spos,
		.span = span,
		.attr = {
			.rid = chain->attr.sep.rid,
			.max_weight = chain->attr.sep.weight + TM_WEIGHT_MARGIN		/* x8 */
		},

		/* stats */
		.score = aln->score,		/* copy */
		.aln = aln,					/* save original */

		/* path */
		.path = { .ptr = aln->path, .len = aln->path_length }
	};

	kv_push(tm_aln_t, self->extend.arr, a);
	rbt_insert_aln(kv_ptr(self->extend.arr), kv_cnt(self->extend.arr) - 1);
	return;
}

static _force_inline
void tm_extend_clear(tm_scan_t *self)
{
	/* clear test threshold */
	self->extend.max_weight = 32;

	/* clear rbt */
	kv_clear(self->extend.arr);
	kv_pushp(tm_aln_t, self->extend.arr);
	rbt_init_static_aln(kv_ptr(self->extend.arr));

	/* clear working stack */
	dz_arena_flush(self->extend.trace);
	return;
}

static _force_inline
tm_alnv_t *tm_extend_finalize(tm_scan_t *self)
{
	/* results */
	tm_aln_t const *src = kv_ptr(self->extend.arr);
	size_t const cnt = kv_cnt(self->extend.arr) - 1;

	/* allocate dst array */
	tm_alnv_t *dst = dz_arena_malloc(self->extend.trace,
		  sizeof(tm_alnv_t)			/* header */
		+ sizeof(tm_aln_t) * cnt	/* payload */
	);
	debug("src(%p), cnt(%zu), dst(%p)", src, cnt, dst);

	/* sort by qpos: traverse tree */
	tm_aln_t const all = {
		.pos  = { .q = 0 },
		.span = { .q = INT32_MAX },
		// .qmax = INT32_MAX
	};

	rbt_iter_t it;
	rbt_init_iter_aln(&it, src, &all);

	tm_aln_t const *p = rbt_fetch_head_aln(&it, src, &all);
	tm_aln_t *q = dst->arr;
	while(p != NULL) {
		debug("p(%p), pos(%u, %u), span(%u, %u)", p, p->pos.q, p->pos.r, p->span.q, p->span.r);

		*q++ = *p;
		p = rbt_fetch_next_aln(&it, src, &all);
	}
	size_t const saved = q - dst->arr;
	dst->cnt = saved;

	debug("cnt(%zu, %zu)", cnt, saved);
	return(dst);
}



/* reference fetcher and query converter for dz */
typedef struct {
	__m128i mv;
	uint8_t const *p;
	int64_t inc;
	uint64_t conv;
} tm_extend_fetcher_t;

static __dz_vectorize
void tm_extend_fetcher_init(tm_extend_fetcher_t *self, uint8_t const *ref, uint32_t dir)
{
	self->p    = ref - dir;
	self->inc  = dir ? -1LL : 1LL;
	self->conv = dir ? 0xf7b3d591e6a2c480 : 0xfedcba9876543210;

/*
	uint8_t const *p = ref - dir;
	while(*p != '\0') {
		fprintf(stderr, "%c", "NACMGRSVTWYHKDBN"[*p]);
		p += self->inc;
	}
	fprintf(stderr, "\n");
*/
	return;
}

static __dz_vectorize
dz_fill_fetch_t tm_extend_fetch_next(tm_extend_fetcher_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(query);

	if(*self->p == '\0') {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.ch      = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = (self->conv>>(4 * c)) & 0x0f;
	self->mv = _mm_cvtsi64_si128(_loadu_u32(&score_matrix[e * DZ_QUERY_MAT_SIZE]));	/* <<2 */
	// debug("p(%p), c(%x, %c), e(%x), inc(%ld)", self->p, c, "NACMGRSVTWYHKDBN"[c], e, self->inc);
	// _print_v16i8((v16i8_t){ self->mv });

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.ch      = e
	});
}

static __dz_vectorize
__m128i tm_extend_get_profile(tm_extend_fetcher_t *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	dz_unused(score_matrix);

	uint8_t const *packed = dz_query_packed_array(query);
	__m128i const v = ({
		__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
		__m128i const sc = _mm_shuffle_epi8(self->mv, qv);
		// _print_v16i8((v16i8_t){ qv });
		_mm_cvtepi8_epi16(sc);
	});
	// _print_v16i8((v16i8_t){ v });
	return(v);
}

static __dz_vectorize
dz_trace_match_t tm_extend_get_match(int8_t const *score_matrix, dz_query_t const *query, size_t qidx, uint32_t ch)
{
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = score_matrix[ch * DZ_QUERY_MAT_SIZE + packed[qidx]],
		.match = packed[qidx] == ch
	});
}

static
size_t tm_extend_calc_dim(size_t qlen)
{
	dz_unused(qlen);
	return(1);
}

static __dz_vectorize
__m128i tm_extend_conv(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		[nA]     = 0x00, [nC]     = 0x01, [nG]     = 0x02, [nT]     = 0x03,	/* forward */
		[nA + 1] = 0x03, [nC + 1] = 0x02, [nG + 1] = 0x01, [nT + 1] = 0x00,	/* reverse-complemented */
		[0x0f]   = 0x04		/* invalid */
	};
	v16i8_t const cv = _load_v16i8(conv);
	v16i8_t const d = _set_v16i8(dir);
	v16i8_t const x = _or_v16i8((v16i8_t){ v }, d);
	v16i8_t const y = _shuf_v16i8(cv, x);

	// _print_v16i8((v16i8_t){ v });
	// _print_v16i8(y);
	return(y.v1);
}

static _force_inline
dz_query_t *tm_pack_query_wrap(dz_arena_t *mem, dz_profile_t const *profile, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, tm_pair_t dir, tm_pair_t pos)
{
	dz_pack_query_t const pack = {
		.dir      = 0,		/* ignored */
		.invalid  = 0x0f,
		.calc_dim = tm_extend_calc_dim,
		.conv     = tm_extend_conv
	};

	/* calc reference remaining length */
	size_t const rlen = tm_idx_ref_seq_len(sk);
	size_t const rrem = dir.r ? pos.r : rlen - pos.r;

	if(dir.q) {
		size_t const qrem = MIN2(2 * rrem, pos.q);
		debug("reverse: qpos(%u) --- qspan(%u) --> qspos(%u)", pos.q, qrem, pos.q - qrem);

		dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, &query[pos.q - qrem], qrem, &pack);
		dz_pack_query_reverse_core(q, profile, &query[pos.q - qrem], qrem, &pack);
		return(q);
	} else {
		size_t const qrem = MIN2(2 * rrem, qlen - pos.q);
		debug("forward: qpos(%u) --- qspan(%u) --> qepos(%u)", pos.q, qrem, pos.q + qrem);

		dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, &query[pos.q], qrem, &pack);
		dz_pack_query_forward_core(q, profile, &query[pos.q], qrem, &pack);
		return(q);
	}
}

static _force_inline
dz_state_t const *tm_extend_wrap(dz_arena_t *mem, dz_profile_t const *profile, tm_idx_sketch_t const *sk, dz_query_t const *q, uint32_t rdir, tm_pair_t pos)
{
	uint8_t const *ref = tm_idx_ref_seq_ptr(sk);

	/* use 4bit fetcher */
	tm_extend_fetcher_t w __attribute__(( aligned(16) ));
	tm_extend_fetcher_init(&w, &ref[pos.r], rdir);
	debug("ref(%p, %p, %zu), rlen(%zu), rdir(%u), rpos(%u)", ref, &ref[pos.r], &ref[pos.r] - ref, tm_idx_ref_seq_len(sk), rdir, pos.r);

	dz_fetcher_t fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fill_fetch_next_t)tm_extend_fetch_next,
		.get_profile = (dz_fill_get_profile_t)tm_extend_get_profile,		/* direct conversion */
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, q, &fetcher, (dz_state_t const **)&profile->root, 1));
}

static _force_inline
tm_pair_t tm_calc_max_wrap(dz_state_t const *r, uint32_t rdir, tm_pair_t rpos)
{
	if(r == NULL || r->max.cap == NULL) {
		return(rpos);
	}

	/* get downward max */
	dz_max_pos_t const s = dz_calc_max_pos_core(r);
	return((tm_pair_t){
		.r = rpos.r + (rdir ? -s.rpos : s.rpos),
		.q = rpos.q - s.qpos
	});
}

typedef struct {
	tm_pair_t spos;
	dz_alignment_t const *aln;
} tm_extend_res_t;

static _force_inline
tm_extend_res_t tm_extend_core(tm_scan_t *self, tm_idx_profile_t const *pf, tm_idx_sketch_t const *sk, uint8_t const *query, size_t qlen, tm_chain_t const *c)
{
	/* flush working memory */
	dz_arena_flush(self->extend.fill);
	tm_extend_res_t const failed = { .aln  = NULL };

	/* load root position */
	uint32_t const rdir = c->dir;
	tm_pair_t const rpos = tm_chain_pos(c);
	debug("reverse: rdir(%u), rpos(%u, %u)", rdir, rpos.q, rpos.r);

	/* reference reverse */
	dz_query_t *qr = tm_pack_query_wrap(self->extend.fill, pf->extend.dz, sk, query, qlen, (tm_pair_t){ .r = rdir ^ 0x01, .q = 1 }, rpos);
	dz_state_t const *r = tm_extend_wrap(self->extend.fill, pf->extend.dz, sk, qr, rdir ^ 0x01, rpos);
	tm_pair_t const spos = tm_calc_max_wrap(r, rdir ^ 0x01, rpos);
	debug("reverse: rpos(%u, %u) --- score(%d) --> spos(%u, %u)", rpos.q, rpos.r, r != NULL ? r->max.score : 0, spos.q, spos.r);

	/* forward */
	dz_query_t *qf = tm_pack_query_wrap(self->extend.fill, pf->extend.dz, sk, query, qlen, (tm_pair_t){ .r = rdir, .q = 0 }, spos);
	dz_state_t const *f = tm_extend_wrap(self->extend.fill, pf->extend.dz, sk, qf, rdir, spos);
	debug("forward: score(%d)", (f != NULL ? f->max.score : -1));
	if(f == NULL || f->max.score == 0) { return(failed); }

	/* traceback */
	dz_alignment_t const *aln = dz_trace_core(self->extend.trace, pf->extend.dz, f, (dz_trace_get_match_t)tm_extend_get_match);
	debug("forward: spos(%u, %u) --- score(%d) --> epos(%u, %u), f(%p), aln(%p)", spos.q, spos.r, aln->score, aln->query_length + spos.q, aln->ref_length + spos.r, f, aln);
	return((tm_extend_res_t){
		.spos = spos,
		.aln  = aln
	});
}

static _force_inline
size_t tm_extend_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *query, size_t qlen, tm_chain_t const *chain, size_t ccnt)
{
	tm_idx_sketch_t const **sk  = (tm_idx_sketch_t const **)idx->sketch.arr;
	tm_idx_profile_t const **pf = (tm_idx_profile_t const **)idx->profile.arr;

	/* clear bin */
	tm_extend_clear(self);
/*
	for(size_t i = 0; i < qlen; i++) {
		fprintf(stderr, "%c", "A   C   G   T   "[query[i]]);
	}
	fprintf(stderr, "\n");
	for(size_t j = 0; j < DZ_REF_MAT_SIZE; j++) {
		for(size_t i = 0; i < DZ_QUERY_MAT_SIZE; i++) {
			fprintf(stderr, "ref(%zu), query(%zu), score(%d)\n", j, i, pf[0]->extend.dz->matrix[j * DZ_QUERY_MAT_SIZE + i]);
		}
	}
*/

	/* for each chain try extension */
	for(size_t i = 0; i < ccnt; i++) {
		tm_chain_t const *q = &chain[i];
		if(q->attr.sep.weight >= self->extend.max_weight) { break; }

		/* skip if already covered (and there is no possibility that this chain surpasses the previous ones) */
		if(tm_extend_is_covered(self, q)) { continue; }

		// tm_aln_t const caln = tm_chain_as_aln(q);
		tm_idx_sketch_t const *s  = sk[q->attr.sep.rid];
		tm_idx_profile_t const *p = pf[s->h.pid];
		debug("%r, rid(%u), pid(%u)", tm_chain_to_str, q, q->attr.sep.rid, s->h.pid);

		/* extend */
		tm_extend_res_t r = tm_extend_core(self, p, s, query, qlen, q);
		if(r.aln == NULL || r.aln->score <= p->extend.min_score) { continue; }

		tm_extend_record(self, q, r.spos, r.aln);
		if(tm_extend_is_complete(self)) { break; }
	}

	/* anything else? */
	return(kv_cnt(self->extend.arr) - 1);
}


/* evaluate all; query sequence must be encoded in 2bit and shorter than 2Gbp */
static _force_inline
tm_alnv_t *tm_scan_all(tm_scan_t *self, tm_idx_t const *idx, uint8_t const *seq, size_t slen)
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
	if(tm_extend_all(self, idx, seq, slen, cptr, ccnt) == 0) {
		return(NULL);		/* alignment not found */
	}

	/* copy alignment to dst buffer */
	return(tm_extend_finalize(self));
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

static _force_inline
void tm_print_cigar(tm_print_t *self, uint8_t const *path, size_t len)
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
			ZCNT_RESULT size_t cnt = _tzcnt_u64(~mask);
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

static _force_inline
void tm_print_aln(tm_print_t *self, tm_idx_sketch_t const **si, bseq_meta_t const *query, tm_aln_t const *aln)
{
	_unused(self);

	tm_idx_sketch_t const *ref = si[aln->attr.rid];
	printf("%.*s\t%u\t%u\t%.*s\t%u\t%u\t%u\t",
		(int)bseq_name_len(query),
		bseq_name(query),
		aln->pos.q,
		aln->span.q,
		(int)tm_idx_ref_name_len(ref),
		tm_idx_ref_name_ptr(ref),
		aln->pos.r,
		aln->span.r,
		aln->dir
	);
	tm_print_cigar(self, aln->path.ptr, aln->path.len);
	printf("\n");
	return;
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

typedef struct {
	tm_idx_t const *mi;
	bseq_file_t *fp;
	tm_print_t *printer;
	pt_t *pt;

	/* counters */
	size_t icnt, ocnt;
	kvechq_t(mm_mtscan_drain_t) hq;

	tm_scan_t scan[];
} tm_mtscan_t;

static _force_inline
tm_mtscan_t *tm_mtscan_init(tm_idx_t const *mi, tm_print_t *printer, pt_t *pt)
{
	size_t const size = sizeof(tm_mtscan_t) + pt_nth(pt) * sizeof(tm_scan_t);
	tm_mtscan_t *self = malloc(size);
	*self = (tm_mtscan_t){
		.mi = mi,
		.printer = printer,
		.pt = pt
	};

	for(size_t i = 0; i < pt_nth(pt); i++) {
		tm_scan_init_static(&self->scan[i]);
	}
	return(self);
}

static _force_inline
void tm_mtscan_destroy(tm_mtscan_t *self)
{
	if(self == NULL) { return; }

	for(size_t i = 0; i < pt_nth(self->pt); i++) {
		tm_scan_destroy_static(&self->scan[i]);
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
	tm_scan_t *scan = &self->scan[tid];
	bseq_meta_t *meta = bseq_meta_ptr(&batch->seq_bin);

	/* setup result buffer */
	scan->extend.trace = batch->mem;

	/* for each query sequence */
	for(size_t i = 0; i < bseq_meta_cnt(&batch->seq_bin); i++) {
		bseq_meta_t *seq = &meta[i];

		/* print info */
		debugblock({
			fprintf(stderr, "begin, tid(%u), i(%zu), len(%zu), seq(%s)\n", tid, i, bseq_seq_len(seq), bseq_name(seq));
		});

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


static _force_inline
int tm_mtscan_file(tm_mtscan_t *self, char const *fn)
{
	bseq_conf_t const bseq_conf = {
		.batch_size  = BATCH_SIZE,
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
	conf->idxdump = opt_strdup(opt, arg, mm_strlen(arg));
}
static void tm_conf_profile(opt_t *opt, tm_conf_t *conf, char const *arg) {
	_unused(opt);
	conf->profile = opt_strdup(opt, arg, mm_strlen(arg));
}

/* indexing */
static void tm_conf_kmer(opt_t *opt, tm_conf_t *conf, char const *arg) {
	size_t k = mm_atoi(arg, 0);
	oassert(opt, k >= 3 && k <= 7, "k-mer size (-k) must be >= 3 and <= 7.");

	conf->fallback.kmer = k;
}
static void tm_conf_chain(opt_t *opt, tm_conf_t *conf, char const *arg) {
	size_t w = mm_atoi(arg, 0);
	oassert(opt, w >= 1 && w <= 256, "chain window size (-w) must be >= 2 and <= 256.");

	conf->fallback.window = w;
}
static void tm_conf_ccnt(opt_t *opt, tm_conf_t *conf, char const *arg) {
	size_t c = mm_atoi(arg, 0);
	oassert(opt, c >= 1 && c <= 1024, "minimum seed count (-c) must be >= 1 and <= 1024.");

	conf->fallback.min_scnt = c;
}
static void tm_conf_match(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t m = mm_atoi(arg, 0);
	oassert(opt, m >= 1 && m <= 7, "match award (-a) must be >= 1 and <= 7 (%s).", arg);

	conf->fallback.match = m;
}
static void tm_conf_mismatch(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t x = mm_atoi(arg, 0);
	oassert(opt, x >= 1 && x <= 7, "mismatch penalty (-b) must be >= 1 and <= 7 (%s).", arg);

	conf->fallback.mismatch = x;
}
static void tm_conf_gap_open(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t gi = mm_atoi(arg, 0);
	oassert(opt, gi >= 1 && gi <= 7, "gap open penalty (-p) must be >= 1 and <= 7 (%s).", arg);

	conf->fallback.gap_open = gi;
}
static void tm_conf_gap_extend(opt_t *opt, tm_conf_t *conf, char const *arg) {
	int64_t ge = mm_atoi(arg, 0);
	oassert(opt, ge >= 1 && ge <= 7, "gap extension penalty (-q) must be >= 1 and <= 7 (%s).", arg);

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
			.kmer = 3,
			.window = 32,
			.match = 2,
			.mismatch = 3,
			.gap_open = 5,
			.gap_extend = 1,
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
			['c'] = { OPT_REQ,  _c(tm_conf_ccnt) },
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
		debug("p(%s)", *p);

		/* generate index */
		tm_idx_t *mi = tm_idx_gen(&conf->fallback, pt, *p, conf->profile, stderr);
		if(mi == NULL) {
			debug("mi == NULL");
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
	if(!mm_startswith(conf->idxdump, "/dev") && !mm_endswith(conf->idxdump, ".tmi")) {
		message(conf->log, "index filename does not end with `.tmi' (added).");
		conf->idxdump = opt_append(&conf->opt, conf->idxdump, ".tmi");
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
	fclose(fp);
	if(error_code) { return(error_code); }

	message(conf->log, "done. total index size (compressed) on disk: %.1f MB.", (double)stat.out / (1024ULL * 1024));
	return(0);
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
	case ERROR_NO_ARG: error("argument is not enough. a reference and at least one query file are required."); break;

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
	if(!mm_endswith(parg[0], ".tmi")) {
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
	if(tm_conf_init_static(&conf, (char const *const *)argv, stderr)) {
		error("error in parsing arguments. abort.");
		goto _main_final;
	}

	/* always print version */
	tm_conf_outfp_t out = tm_conf_get_outfp(&conf);
	message(out.fp, "Version: %s, Build: %s", tm_version(), TM_ARCH_NAME);

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
