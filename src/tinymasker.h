
/**
 * @file tinymasker.h
 *
 * @author Hajime Suzuki
 * @license MIT
 */


#ifndef _TM_H_INCLUDED
#define _TM_H_INCLUDED

#ifndef _UTILS_H_INCLUDED
#  error "#include \"utils.h\" must be before #include \"tinymasker.h\""
#endif


/* global configurations */
#define TM_BATCH_SIZE			( 512ULL * 1024 )
#define TM_MAX_THREADS			( 256 )				/* FIXME: can be larger? */


/* index data structure */
#define TM_MIN_KMER				( 3 )
#define TM_MAX_KMER				( 9 )
#define TM_REF_ALPH_SIZE		( 16 )
#define TM_QUERY_ALPH_SIZE		( 4 )
#define TM_MAX_SKETCH_TRIAL		( 3 )
#define TM_MAX_AMB_COLUMNS		( 6 )
#define TM_SEED_POS_OFS			( 1 )				/* 1-origin coorrdinates for seed positions; zero is reserved for invalid letter */


/* dozeu configurations */
#define DZ_REF_MAT_SIZE			( TM_REF_ALPH_SIZE )
#define DZ_QUERY_MAT_SIZE		( 2 * TM_QUERY_ALPH_SIZE )	/* with / without bonus */
#define DZ_WRAPPED_API			( 0 )
#define DZ_CIGAR_OP				0x44494d4d
#define DZ_TRANSPOSE_MATRIX		( 1 )
#define DZ_NONGRAPH
#define DZ_OFFSET
// #define DZ_FULL_LENGTH_BONUS	( 1 )				/* query-side bonus disabled */


/* version and arch */
#ifndef TM_VERSION
#  define TM_VERSION			"tinymasker-0.0.1"
#endif
#ifndef TM_COMMIT
#  define TM_COMMIT				"unknown"
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

static _force_inline char const *tm_commit(void) { return(TM_COMMIT); }
static _force_inline char const *tm_arch_name(void) { return(TM_ARCH_NAME); }


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
enum alphabet_query {
	nA = 0x00,
	nC = 0x04,
	nG = 0x08,
	nT = 0x0c
};

enum alphabet_reference {
	tA = 0x0c, tT = 0x03,
	tC = 0x0a, tG = 0x05,

	tR = 0x0d, tY = 0x02,
	tS = 0x0f,	/* 0x00; does not appear in canonical sequence */
	tW = 0x0e,	/* 0x01 */
	tK = 0x0b, tM = 0x04,

	tB = 0x09, tV = 0x06,
	tD = 0x08, tH = 0x07
};
#define tm_rch_pack(x)			( (x)<<1 )
#define tm_rch_unpack(x)		( (x)>>1 )


/* encoding conversion */
static _force_inline
char tm_qch_to_ascii(uint8_t base)
{
	return("A   C   G   T   "[base & 0x0f]);
}

static _force_inline
uint8_t tm_2bit_to_qch(uint8_t q2)
{
	return(q2<<2);
}

static _force_inline
uint64_t tm_rch_is_amb(uint8_t base)
{
	uint64_t const magic = 0x1110101111010111;
	return((magic>>(2 * base)) & 0x01);
}

static _force_inline
char tm_rch_to_ascii(uint8_t base)
{
	return("SWYTMGVHDBCKARWS"[tm_rch_unpack(base) & 0x0f]);
}

static _force_inline
uint64_t tm_rch_to_2bit(uint8_t base)
{
	/* convert iupac nucleotide to pair of 2-bit encoded bases */
	uint64_t const magic = 0x9c80e1e8d4243dc9;
	return((magic>>(2 * base)) & 0x0f);
}

unittest( .name = "base.conv.magic" ) {
	struct { uint8_t base, amb, enc; } const x[16] = {
		{ A, 0, tm_rch_pack(tA) },
		{ C, 0, tm_rch_pack(tC) },
		{ G, 0, tm_rch_pack(tG) },
		{ T, 0, tm_rch_pack(tT) },
		{ R, 1, tm_rch_pack(tR) },
		{ Y, 1, tm_rch_pack(tY) },
		{ S, 1, tm_rch_pack(tS) },
		{ S, 1, tm_rch_pack(tS ^ 0x0f) },
		{ W, 1, tm_rch_pack(tW) },
		{ W, 1, tm_rch_pack(tW ^ 0x0f) },
		{ K, 1, tm_rch_pack(tK) },
		{ M, 1, tm_rch_pack(tM) },
		{ B, 1, tm_rch_pack(tB) },
		{ D, 1, tm_rch_pack(tD) },
		{ H, 1, tm_rch_pack(tH) },
		{ V, 1, tm_rch_pack(tV) }
	};
	uint8_t const conv[4] = { A, C, G, T };

	for(size_t i = 0; i < 16; i++) {
		ut_assert(tm_rch_is_amb(x[i].enc) == x[i].amb);
		ut_assert((conv[tm_rch_to_2bit(x[i].enc) & 0x03] & x[i].base) != 0);
		ut_assert((conv[tm_rch_to_2bit(x[i].enc) & 0x03] & ~x[i].base) == 0);

		if(x[i].amb) {
			ut_assert((conv[tm_rch_to_2bit(x[i].enc)>>2] & x[i].base) != 0);
			ut_assert((conv[tm_rch_to_2bit(x[i].enc)>>2] & ~x[i].base) == 0);
		}
	}
}

/* calc matching state between query and reference bases */
static _force_inline
uint64_t tm_qrch_is_match(uint8_t q, uint8_t r)
{
	/* query in 2bit and reference in 4bit; faster? */
	#define _q(x)			( (x)<<(n##x) )
	#define _rx(x, y)		( ((uint64_t)(x))<<(4 * (y)) )
	#define _r(x)			( _rx((x), (t##x)) )
	uint64_t const qconv = _q(A) | _q(C) | _q(G) | _q(T);
	uint64_t const rconv = (
		  _r(A) | _r(C) | _r(G) | _r(T)
		| _r(R) | _r(Y) | _r(S) | _rx(S, tS ^ 0x0f)
		| _r(K) | _r(M) | _r(W) | _rx(W, tW ^ 0x0f)
		| _r(B) | _r(D) | _r(H) | _r(V)
	);
	return(((qconv>>q) & (rconv>>(r<<1)) & 0x0f) != 0);
}

#if defined(UNITTEST) && (UNITTEST != 0)
static _force_inline
uint64_t tm_qrch_is_match_naive(uint8_t q, uint8_t r)
{
	/* standard impl */
	static uint8_t const qconv[TM_QUERY_ALPH_SIZE] = {
		[nA>>2] = A,
		[nC>>2] = C,
		[nG>>2] = G,
		[nT>>2] = T
	};
	static uint8_t const rconv[TM_REF_ALPH_SIZE] = {
		[tA] = A, [tC] = C, [tG] = G, [tT] = T,
		[tR] = R, [tY] = Y, [tS] = S, [tS ^ 0x0f] = S,
		[tK] = K, [tM] = M, [tW] = W, [tW ^ 0x0f] = W,
		[tB] = B, [tD] = D, [tH] = H, [tV] = V
	};
	return((qconv[q>>2] & rconv[tm_rch_unpack(r)]) != 0);
}

unittest( .name = "base.match" ) {
	uint8_t const qconv[4] = { nA, nC, nG, nT };
	for(size_t qidx = 0; qidx < 4; qidx++) {
		uint8_t const q = qconv[qidx];

		for(size_t ridx = 0; ridx < 16; ridx++) {
			uint8_t const r = tm_rch_pack(ridx);
			ut_assert(tm_qrch_is_match(q, r) == tm_qrch_is_match_naive(q, r), "q(%x), r(%x), match(%u, %u)", q, r, tm_qrch_is_match(q, r), tm_qrch_is_match_naive(q, r));
		}
	}
}
#endif


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


/* debug supplementary */
#define _print_v16u8(a) { \
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

#define _print_v32i8x(a) { \
	uint8_t _buf[32]; \
	_storeu_v32i8(_buf, a); \
	debug("(v32i8_t) %s(%x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x)", \
		#a, \
		_buf[16 + 15], \
		_buf[16 + 14], \
		_buf[16 + 13], \
		_buf[16 + 12], \
		_buf[16 + 11], \
		_buf[16 + 10], \
		_buf[16 + 9], \
		_buf[16 + 8], \
		_buf[16 + 7], \
		_buf[16 + 6], \
		_buf[16 + 5], \
		_buf[16 + 4], \
		_buf[16 + 3], \
		_buf[16 + 2], \
		_buf[16 + 1], \
		_buf[16 + 0], \
		_buf[15], \
		_buf[14], \
		_buf[13], \
		_buf[12], \
		_buf[11], \
		_buf[10], \
		_buf[9], \
		_buf[8], \
		_buf[7], \
		_buf[6], \
		_buf[5], \
		_buf[4], \
		_buf[3], \
		_buf[2], \
		_buf[1], \
		_buf[0]); \
}


#endif	/* _TM_H_INCLUDED */
/**
 * end of tinymasker.h
 */
