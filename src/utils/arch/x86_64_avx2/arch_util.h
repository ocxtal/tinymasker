
/**
 * @file arch_util.h
 *
 * @brief architecture-dependent utilities devided from util.h
 */
#ifndef _ARCH_UTIL_H_INCLUDED
#define _ARCH_UTIL_H_INCLUDED
#define ARCH_NAME			"x86-64::AVX2"

#include "vector.h"
#include <x86intrin.h>
#include <stdlib.h>
#include <stdint.h>

/**
 * misc bit operations (popcnt, tzcnt, and lzcnt)
 */

/**
 * @macro popcnt
 */
#define _popc_u64(x)		( (size_t)_mm_popcnt_u64(x) )
#define _popc_u32(x)		( (size_t)_mm_popcnt_u32(x) )

/**
 * @macro ZCNT_RESULT
 * @brief workaround for a bug in gcc (<= 5), all the results of tzcnt / lzcnt macros must be modified by this label
 */
#ifndef ZCNT_RESULT
#  if defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 600
#    define ZCNT_RESULT		volatile
#  else
#    define ZCNT_RESULT
#  endif
#endif

/**
 * @macro tzcnt
 * @brief trailing zero count (count #continuous zeros from LSb)
 */
/** immintrin.h is already included */
#if defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 490
#  define _tzc_u64(x)		( (size_t)__tzcnt_u64(x) )
#  define _tzc_u32(x)		( (size_t)__tzcnt_u32(x) )
#else
#  define _tzc_u64(x)		( (size_t)_tzcnt_u64(x) )
#  define _tzc_u32(x)		( (size_t)_tzcnt_u32(x) )
#endif

/**
 * @macro lzcnt
 * @brief leading zero count (count #continuous zeros from MSb)
 */
/* __lzcnt_u64 in bmiintrin.h gcc-4.6, _lzcnt_u64 in lzcntintrin.h from gcc-4.7 */
#if defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 470
#  define _lzc_u64(x)		( (size_t)__lzcnt_u64(x) )
#  define _lzc_u32(x)		( (size_t)__lzcnt_u32(x) )
#else
#  define _lzc_u64(x)		( (size_t)_lzcnt_u64(x) )
#  define _lzc_u32(x)		( (size_t)_lzcnt_u32(x) )
#endif

/**
 * @macro _swap_u64
 */
#if defined(__clang__) || (defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 470)
#  define _swap_u64(x)		({ uint64_t _x = (x); __asm__( "bswapq %0" : "+r"(_x) ); _x; })
#  define _swap_u32(x)		({ uint32_t _x = (x); __asm__( "bswapl %0" : "+r"(_x) ); _x; })
#else
#  define _swap_u64(x)		( (uint64_t)_bswap64(x) )
#  define _swap_u32(x)		( (uint32_t)_bswap32(x) )
#endif

/**
 * @macro _loadu_u64, _storeu_u64
 */
#define _loadu_u64(p)		({ uint8_t const *_p = (uint8_t const *)(p); *((uint64_t const *)_p); })
#define _storeu_u64(p, e)	({ uint8_t *_p = (uint8_t *)(p); *((uint64_t *)(_p)) = (e); sizeof(uint64_t); })
#define _loadu_u32(p)		({ uint8_t const *_p = (uint8_t const *)(p); *((uint32_t const *)_p); })
#define _storeu_u32(p, e)	({ uint8_t *_p = (uint8_t *)(p); *((uint32_t *)(_p)) = (e); sizeof(uint32_t); })
#define _loadu_u16(p)		({ uint8_t const *_p = (uint8_t const *)(p); *((uint16_t const *)_p); })
#define _storeu_u16(p, e)	({ uint8_t *_p = (uint8_t *)(p); *((uint16_t *)(_p)) = (e); sizeof(uint16_t); })

/**
 * @macro _unpack_v16i8, _unpack_v32i8
 * @brief 4bit packed -> uint8_t array
 */
#define _unpack_v16i8(_v) ({ \
	__m128i t0 = (_v).v1, t1 = _mm_cvtepu8_epi16(t0); \
	t1 = _mm_and_si128(_mm_or_si128(t1, _mm_slli_epi16(t1, 4)), _mm_set1_epi8(0x0f)); \
	((v16i8_t){ .v1 = t1 }); \
})
#define _unpack_v32i8(_v) ({ \
	__m128i t0 = (_v).v1; \
	__m256i t1 = _mm256_cvtepu8_epi16(t0); \
	t1 = _mm256_and_si256(_mm256_or_si256(t1, _mm256_slli_epi16(t1, 4)), _mm256_set1_epi8(0x0f)); \
	((v32i8_t){ .v1 = t1 }); \
})

/**
 * @macro _aligned_block_memcpy
 *
 * @brief copy size bytes from src to dst.
 *
 * @detail
 * src and dst must be aligned to 16-byte boundary.
 * copy must be multipe of 16.
 */
#define _ymm_rd_a(src, n) (ymm##n) = _mm256_load_si256((__m256i *)(src) + (n))
#define _ymm_rd_u(src, n) (ymm##n) = _mm256_loadu_si256((__m256i *)(src) + (n))
#define _ymm_wr_a(dst, n) _mm256_store_si256((__m256i *)(dst) + (n), (ymm##n))
#define _ymm_wr_u(dst, n) _mm256_storeu_si256((__m256i *)(dst) + (n), (ymm##n))
#define _memcpy_blk_intl(dst, src, size, _wr, _rd) { \
	/** duff's device */ \
	uint8_t *_src = (uint8_t *)(src), *_dst = (uint8_t *)(dst); \
	uint64_t const _nreg = 16;		/** #ymm registers == 16 */ \
	uint64_t const _tcnt = (size) / sizeof(__m256i); \
	uint64_t const _offset = ((_tcnt - 1) & (_nreg - 1)) - (_nreg - 1); \
	uint64_t _jmp = _tcnt & (_nreg - 1); \
	uint64_t _lcnt = (_tcnt + _nreg - 1) / _nreg; \
	register __m256i ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7; \
	register __m256i ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15; \
	_src += _offset * sizeof(__m256i); \
	_dst += _offset * sizeof(__m256i); \
	switch(_jmp) { \
		case 0: do { _rd(_src, 0); \
		case 15:     _rd(_src, 1); \
		case 14:     _rd(_src, 2); \
		case 13:     _rd(_src, 3); \
		case 12:     _rd(_src, 4); \
		case 11:     _rd(_src, 5); \
		case 10:     _rd(_src, 6); \
		case 9:      _rd(_src, 7); \
		case 8:      _rd(_src, 8); \
		case 7:      _rd(_src, 9); \
		case 6:      _rd(_src, 10); \
		case 5:      _rd(_src, 11); \
		case 4:      _rd(_src, 12); \
		case 3:      _rd(_src, 13); \
		case 2:      _rd(_src, 14); \
		case 1:      _rd(_src, 15); \
		switch(_jmp) { \
			case 0:  _wr(_dst, 0); \
			case 15: _wr(_dst, 1); \
			case 14: _wr(_dst, 2); \
			case 13: _wr(_dst, 3); \
			case 12: _wr(_dst, 4); \
			case 11: _wr(_dst, 5); \
			case 10: _wr(_dst, 6); \
			case 9:  _wr(_dst, 7); \
			case 8:  _wr(_dst, 8); \
			case 7:  _wr(_dst, 9); \
			case 6:  _wr(_dst, 10); \
			case 5:  _wr(_dst, 11); \
			case 4:  _wr(_dst, 12); \
			case 3:  _wr(_dst, 13); \
			case 2:  _wr(_dst, 14); \
			case 1:  _wr(_dst, 15); \
		} \
				     _src += _nreg * sizeof(__m256i); \
				     _dst += _nreg * sizeof(__m256i); \
				     _jmp = 0; \
			    } while(--_lcnt > 0); \
	} \
}
#define _memcpy_blk_aa(dst, src, len)		_memcpy_blk_intl(dst, src, len, _ymm_wr_a, _ymm_rd_a)
#define _memcpy_blk_au(dst, src, len)		_memcpy_blk_intl(dst, src, len, _ymm_wr_a, _ymm_rd_u)
#define _memcpy_blk_ua(dst, src, len)		_memcpy_blk_intl(dst, src, len, _ymm_wr_u, _ymm_rd_a)
#define _memcpy_blk_uu(dst, src, len)		_memcpy_blk_intl(dst, src, len, _ymm_wr_u, _ymm_rd_u)
#define _memset_blk_intl(dst, a, size, _wr) { \
	uint8_t *_dst = (uint8_t *)(dst); \
	__m256i const ymm0 = _mm256_set1_epi8((int8_t)a); \
	uint64_t i; \
	for(i = 0; i < size / sizeof(__m256i); i++) { \
		_wr(_dst, 0); _dst += sizeof(__m256i); \
	} \
}
#define _memset_blk_a(dst, a, size)			_memset_blk_intl(dst, a, size, _ymm_wr_a)
#define _memset_blk_u(dst, a, size)			_memset_blk_intl(dst, a, size, _ymm_wr_u)


/**
 * substitution matrix abstraction
 */
/* store */
#define _store_sb(_scv, sv16)				{ _store_v32i8((_scv).v1, _from_v16i8_v32i8(sv16)); }

/* load */
#define _load_sb(scv)						( _from_v32i8_n(_load_v32i8((scv).v1)) )

/**
 * gap penalty vector abstraction macros
 */
/* store */
#define _make_gap(_e1, _e2, _e3, _e4) ( \
	(v16i8_t){ _mm_set_epi8( \
		(_e4), (_e4), (_e4), (_e4), \
		(_e3), (_e3), (_e3), (_e3), \
		(_e2), (_e2), (_e2), (_e2), \
		(_e1), (_e1), (_e1), (_e1)) \
	} \
)
#define _store_adjh(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	_store_v32i8((_scv).v3, _from_v16i8_v32i8(_make_gap(_adjh, _adjv, _ofsh, _ofsv))) \
}
#define _store_adjv(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	/* nothing to do */ \
	/*_store_v32i8((_scv).v3, _from_v16i8_v32i8(_make_gap(_adjh, _adjv, _ofsh, _ofsv)))*/ \
}
#define _store_ofsh(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	/* nothing to do */ \
	/* _store_v32i8((_scv).v5, _from_v16i8_v32i8(_make_gap(_adjh, _adjv, _ofsh, _ofsv)))*/ \
}
#define _store_ofsv(_scv, _adjh, _adjv, _ofsh, _ofsv) { \
	/* nothing to do */ \
	/*_store_v32i8((_scv).v5, _from_v16i8_v32i8(_make_gap(_adjh, _adjv, _ofsh, _ofsv)))*/ \
}

/* load */
#define _load_gap(_ptr, _idx) ( \
	(v32i8_t){ _mm256_shuffle_epi32(_mm256_load_si256((__m256i const *)(_ptr)), (_idx)) } \
)

#define _load_adjh(_scv)					( _from_v32i8_n(_load_gap((_scv).v3, 0x00)) )
#define _load_adjv(_scv)					( _from_v32i8_n(_load_gap((_scv).v3, 0x00)) )
#define _load_ofsh(_scv)					( _from_v32i8_n(_load_gap((_scv).v3, 0x55)) )
#define _load_ofsv(_scv)					( _from_v32i8_n(_load_gap((_scv).v3, 0x55)) )
#define _load_gfh(_scv)						( _from_v32i8_n(_load_gap((_scv).v3, 0xaa)) )
#define _load_gfv(_scv)						( _from_v32i8_n(_load_gap((_scv).v3, 0xff)) )
/*
#define _load_adjv(_scv)					( _from_v32i8_n(_load_gap((_scv).v3, 0x55)) )
#define _load_ofsv(_scv)					( _from_v32i8_n(_load_gap((_scv).v3, 0xff)) )
*/


/* compare and swap (cas) */
#if defined(__GNUC__)
#  if (defined(_ARCH_GCC_VERSION) && _ARCH_GCC_VERSION < 470) || (defined(__INTEL_COMPILER) && _ARCH_GCC_COMPAT < 470)
#    define cas(ptr, cmp, val) ({ \
		uint8_t _res; \
		__asm__ volatile ("lock cmpxchg %[src], %[dst]\n\tsete %[res]" \
			: [dst]"+m"(*ptr), [res]"=a"(_res) \
			: [src]"r"(val), "a"(*cmp) \
			: "memory", "cc"); \
		_res; \
	})
#    define fence() ({ \
		__asm__ volatile ("mfence"); \
	})
#  else											/* > 4.7 */
#    define cas(ptr, cmp, val)	__atomic_compare_exchange_n(ptr, cmp, val, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)
#    define fence()				__sync_synchronize()
#  endif
#else
#  error "atomic compare-and-exchange is not supported in this version of compiler."
#endif

#endif /* #ifndef _ARCH_UTIL_H_INCLUDED */
/**
 * end of arch_util.h
 */
