
/**
 * @file v32i8.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V32I8_H_INCLUDED
#define _V32I8_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 32cell */
typedef struct v32i8_s {
	__m128i v1;
	__m128i v2;
} v32i8_t;

/* expanders (without argument) */
#define _e_x_v32i8_1(u)
#define _e_x_v32i8_2(u)

/* expanders (without immediate) */
#define _e_v_v32i8_1(a)				(a).v1
#define _e_v_v32i8_2(a)				(a).v2
#define _e_vv_v32i8_1(a, b)			(a).v1, (b).v1
#define _e_vv_v32i8_2(a, b)			(a).v2, (b).v2
#define _e_vvv_v32i8_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v32i8_2(a, b, c)		(a).v2, (b).v2, (c).v2

/* expanders with immediate */
#define _e_i_v32i8_1(imm)			(imm)
#define _e_i_v32i8_2(imm)			(imm)
#define _e_vi_v32i8_1(a, imm)		(a).v1, (imm)
#define _e_vi_v32i8_2(a, imm)		(a).v2, (imm)
#define _e_vvi_v32i8_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v32i8_2(a, b, imm)	(a).v2, (b).v2, (imm)

/* address calculation macros */
#define _addr_v32i8_1(imm)			( (__m128i *)(imm) )
#define _addr_v32i8_2(imm)			( (__m128i *)(imm) + 1 )
#define _pv_v32i8(ptr)				( _addr_v32i8_1(ptr) )
/* expanders with pointers */
#define _e_p_v32i8_1(ptr)			_addr_v32i8_1(ptr)
#define _e_p_v32i8_2(ptr)			_addr_v32i8_2(ptr)
#define _e_pv_v32i8_1(ptr, a)		_addr_v32i8_1(ptr), (a).v1
#define _e_pv_v32i8_2(ptr, a)		_addr_v32i8_2(ptr), (a).v2

/* expand intrinsic name */
#define _i_v32i8(intrin) 			_mm_##intrin##_epi8
#define _i_v32u8(intrin) 			_mm_##intrin##_epu8
#define _i_v32i8x(intrin)			_mm_##intrin##_si128

/* apply */
#define _a_v32i8(intrin, expander, ...) ( \
	(v32i8_t) { \
		_i_v32i8(intrin)(expander##_v32i8_1(__VA_ARGS__)), \
		_i_v32i8(intrin)(expander##_v32i8_2(__VA_ARGS__)) \
	} \
)
#define _a_v32u8(intrin, expander, ...) ( \
	(v32i8_t) { \
		_i_v32u8(intrin)(expander##_v32i8_1(__VA_ARGS__)), \
		_i_v32u8(intrin)(expander##_v32i8_2(__VA_ARGS__)) \
	} \
)
#define _a_v32i8x(intrin, expander, ...) ( \
	(v32i8_t) { \
		_i_v32i8x(intrin)(expander##_v32i8_1(__VA_ARGS__)), \
		_i_v32i8x(intrin)(expander##_v32i8_2(__VA_ARGS__)) \
	} \
)
#define _a_v32i8xv(intrin, expander, ...) { \
	_i_v32i8x(intrin)(expander##_v32i8_1(__VA_ARGS__)); \
	_i_v32i8x(intrin)(expander##_v32i8_2(__VA_ARGS__)); \
}

/* load and store */
#define _load_v32i8(...)	_a_v32i8x(load, _e_p, __VA_ARGS__)
#define _loadu_v32i8(...)	_a_v32i8x(loadu, _e_p, __VA_ARGS__)
#define _store_v32i8(...)	_a_v32i8xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v32i8(...)	_a_v32i8xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v32i8(...)		_a_v32i8(set1, _e_i, __VA_ARGS__)
#define _zero_v32i8()		_a_v32i8x(setzero, _e_x, _unused)

/* swap (reverse) */
#define _swap_idx_v32i8() ( \
	_mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15) \
)
#define _swap_v32i8(a) ( \
	(v32i8_t) { \
		_mm_shuffle_epi8((a).v2, _swap_idx_v32i8()), \
		_mm_shuffle_epi8((a).v1, _swap_idx_v32i8()) \
	} \
)

/* logics */
#define _not_v32i8(...)		_a_v32i8x(not, _e_v, __VA_ARGS__)
#define _and_v32i8(...)		_a_v32i8x(and, _e_vv, __VA_ARGS__)
#define _or_v32i8(...)		_a_v32i8x(or, _e_vv, __VA_ARGS__)
#define _xor_v32i8(...)		_a_v32i8x(xor, _e_vv, __VA_ARGS__)
#define _andn_v32i8(...)	_a_v32i8x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v32i8(...)		_a_v32i8(add, _e_vv, __VA_ARGS__)
#define _sub_v32i8(...)		_a_v32i8(sub, _e_vv, __VA_ARGS__)
#define _adds_v32i8(...)	_a_v32i8(adds, _e_vv, __VA_ARGS__)
#define _subs_v32i8(...)	_a_v32i8(subs, _e_vv, __VA_ARGS__)
#define _addus_v32i8(...)	_a_v32u8(adds, _e_vv, __VA_ARGS__)
#define _subus_v32i8(...)	_a_v32u8(subs, _e_vv, __VA_ARGS__)
#define _max_v32i8(...)		_a_v32i8(max, _e_vv, __VA_ARGS__)
#define _min_v32i8(...)		_a_v32i8(min, _e_vv, __VA_ARGS__)
#define _maxu_v32i8(...)	_a_v32u8(max, _e_vv, __VA_ARGS__)
#define _minu_v32i8(...)	_a_v32u8(min, _e_vv, __VA_ARGS__)

/* shuffle */
#define _shuf_v32i8(...)	_a_v32i8(shuffle, _e_vv, __VA_ARGS__)

/* blend */
#define _sel_v32i8(mask, a, b) ( \
	(v32i8_t) { \
		_mm_blendv_epi8((b).v1, (a).v1, (mask).v1), \
		_mm_blendv_epi8((b).v2, (a).v2, (mask).v2) \
	} \
)

/* compare */
#define _eq_v32i8(...)		_a_v32i8(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v32i8(...)		_a_v32i8(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _ins_v32i8(a, val, imm) { \
	if((imm) < sizeof(__m128i)) { \
		(a).v1 = _i_v32i8(insert)((a).v1, (val), (imm)); \
	} else { \
		(a).v2 = _i_v32i8(insert)((a).v2, (val), (imm) - sizeof(__m128i)); \
	} \
}
#define _ext_v32i8(a, imm) ( \
	(int8_t)(((imm) < sizeof(__m128i)) ? ( \
		_i_v32i8(extract)((a).v1, (imm)) \
	) : ( \
		_i_v32i8(extract)((a).v2, (imm) - sizeof(__m128i)) \
	)) \
)

/* shift */
#define _bsl_v32i8(a, imm) ( \
	(v32i8_t) { \
		_i_v32i8x(slli)((a).v1, (imm)), \
		_i_v32i8(alignr)((a).v2, (a).v1, sizeof(__m128i) - (imm)) \
	} \
)
#define _bsr_v32i8(a, imm) ( \
	(v32i8_t) { \
		_i_v32i8(alignr)((a).v2, (a).v1, (imm)), \
		_i_v32i8x(srli)((a).v2, (imm)) \
	} \
)

/* double shift (palignr) */
#define _bsld_v32i8(a, b, imm) ( \
	(v32i8_t) { \
		_i_v32i8(alignr)((a).v1, (b).v2, sizeof(__m128i) - (imm)), \
		_i_v32i8(alignr)((a).v2, (a).v1, sizeof(__m128i) - (imm)) \
	} \
)
#define _bsrd_v32i8(a, b, imm) ( \
	(v32i8_t) { \
		_i_v32i8(alignr)((b).v2, (b).v1, (imm)), \
		_i_v32i8(alignr)((a).v1, (b).v2, (imm)) \
	} \
)

/* bit shift */
#define _shl_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm_slli_epi32((a).v1, (imm)), \
		_mm_slli_epi32((a).v2, (imm)) \
	} \
)
#define _shr_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm_srli_epi32((a).v1, (imm)), \
		_mm_srli_epi32((a).v2, (imm)) \
	} \
)
#define _sal_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm_slai_epi32((a).v1, (imm)), \
		_mm_slai_epi32((a).v2, (imm)) \
	} \
)
#define _sar_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm_srai_epi32((a).v1, (imm)), \
		_mm_srai_epi32((a).v2, (imm)) \
	} \
)

/* mask */
#define _mask_v32i8(a) ( \
	(v32_mask_t) { \
		.m1 = _i_v32i8(movemask)((a).v1), \
		.m2 = _i_v32i8(movemask)((a).v2) \
	} \
)

/* convert */
#define _cvt_v32i16_v32i8(a) ( \
	(v32i8_t) { \
		_mm_packs_epi16((a).v1, (a).v2), \
		_mm_packs_epi16((a).v3, (a).v4) \
	} \
)
#define _cvt_v32i16_v32u8(a) ( \
	(v32i8_t) { \
		_mm_packus_epi16((a).v1, (a).v2), \
		_mm_packus_epi16((a).v3, (a).v4) \
	} \
)

/* debug print */
#ifdef DEBUG
#define _print_v32i8(a) { \
	int8_t _buf[32]; \
	_storeu_v32i8(_buf, a); \
	debug("(v32i8_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_buf[31], \
		_buf[30], \
		_buf[29], \
		_buf[28], \
		_buf[27], \
		_buf[26], \
		_buf[25], \
		_buf[24], \
		_buf[23], \
		_buf[22], \
		_buf[21], \
		_buf[20], \
		_buf[19], \
		_buf[18], \
		_buf[17], \
		_buf[16], \
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
#else
#define _print_v32i8(x)		;
#endif

#endif /* _V32I8_H_INCLUDED */
/**
 * end of v32i8.h
 */
