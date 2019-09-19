
/**
 * @rbread.h
 * @brief gzread-compatible stream reader
 */

#ifndef _RBREAD_H_INCLUDED
#define _RBREAD_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"
#include "arch.h"
#include "wmalloc.h"


#if !defined(RB_USE_ZLIB) && !defined(RB_DONTUSE_ZLIB)
#  define RB_USE_ZLIB
#endif
#if !defined(RB_USE_BZLIB) && !defined(RB_DONTUSE_BZLIB)
#  define RB_USE_BZLIB
#endif
#if !defined(RB_USE_LZMA) && !defined(RB_DONTUSE_LZMA)
#  define RB_USE_LZMA
#endif


#ifdef RB_USE_ZLIB
#  include <zlib.h>
#endif
#ifdef RB_USE_BZLIB
#  include <bzlib.h>
#endif
#ifdef RB_USE_LZMA
#  include <lzma.h>
#endif


#if 1
#define RB_BUF_SIZE				( 2ULL * 1024 * 1024 )
#define RB_BUF_ALIGN_SIZE		( ARCH_HUGEPAGE_SIZE )
#else
#define RB_BUF_SIZE				( 1024 )
#define RB_BUF_ALIGN_SIZE		( 1024 )
#endif


#define rbdebug(_x, ...)		;


/**
 * @struct rbread_s
 * @brief context (gzFile)
 */
typedef struct {
	/* (%rax, %rdx) pair */
	size_t obtained;
	size_t remaining;
} rb_reader_res_t;

typedef struct rbread_s {
	FILE *fp;
	uint8_t *ibuf, *buf;	/* input buffer and output buffer; ibuf is unused for transparent mode */
	size_t head, tail;
	size_t bulk_size, buf_size;

	/*
	 * end-of-file status:
	 *   0: stream continues
	 *   1: input stream starved, input buffer remains
	 *   2: input buffer starved, output buffer remains
	 *   3: output buffer starved, no more bytes
	 */
	uint32_t eof, _pad;

	rb_reader_res_t (*read)(struct rbread_s *, void *, size_t);
	union {					/* exclusive */
		#ifdef RB_USE_ZLIB
			z_stream z;		/* .gz */
		#endif
		#ifdef RB_USE_BZLIB
			bz_stream b;	/* .bz2 */
		#endif
		#ifdef RB_USE_LZMA
			lzma_stream x;	/* .xz */
		#endif

		uint64_t _pad;
	} s;
} rbread_t;
#define rbeof(_rb)				( (_rb)->eof >= 3 )


/**
 * @fn rb_read_gzip, rb_read_bz2, rb_read_xz, rb_read_transparent
 * @brief block fetcher
 */
#define _rb_reader(_sfx, _ctx, _type_t, _body) \
	static inline rb_reader_res_t rb_read_##_sfx(rbread_t *rb, void *dst, size_t len) \
	{ \
		/* eof == 0 indicates input stream remains */ \
		if(rb->s._ctx.avail_in == 0 && rb->eof == 0) { \
			/* copy to the head then fill */ \
			/* uint64_t const _fragment = _loadu_u64(rb->s._ctx.next_in); _storeu_u64(rb->ibuf, _fragment); */ \
			rb->s._ctx.next_in   = (_type_t *)&rb->ibuf[rb->s._ctx.avail_in]; \
			rb->s._ctx.avail_in += fread((void *)rb->s._ctx.next_in, sizeof(char), rb->bulk_size, rb->fp); \
			/* input stream starved, mark eof */ \
			if(feof(rb->fp)) { rb->eof = 1; } \
			rbdebug("refill, avail_in(%d), eof(%u, %u)", rb->s._ctx.avail_in, feof(rb->fp), rb->eof); \
		} \
		rb->s._ctx.next_out  = (_type_t *)dst; \
		rb->s._ctx.avail_out = len; \
		rbdebug("len(%zu), remaining(%zu)", len, (size_t)rb->s._ctx.avail_in); \
		_body; \
		rbdebug("len(%zu), obtained(%zu), remaining(%zu)", len, (size_t)(len - rb->s._ctx.avail_out), (size_t)rb->s._ctx.avail_in); \
		return((rb_reader_res_t){ \
			.obtained  = len - rb->s._ctx.avail_out, \
			.remaining = rb->s._ctx.avail_in \
		}); \
	}

#ifdef RB_USE_ZLIB
_rb_reader(gzip, z, uint8_t, {
	switch(inflate(&rb->s.z, Z_NO_FLUSH)) {
		default:
		rbdebug("unhandled error");
		case Z_STREAM_END: rbdebug("reset stream"); inflateReset(&rb->s.z);
		if(rb->eof == 1) { rbdebug("input buffer starved"); rb->eof = 2; }
		case Z_OK: break;
	}
});
#endif

#ifdef RB_USE_BZLIB
_rb_reader(bz2, b, char, {
	switch(BZ2_bzDecompress(&rb->s.b)) {
		default:
		case BZ_STREAM_END:
		if(rb->eof == 1) { rbdebug("input buffer starved"); rb->eof = 2; }
		case BZ_OK: break;
	}
});
#endif

#ifdef RB_USE_LZMA
_rb_reader(xz, x, uint8_t, {
	switch(lzma_code(&rb->s.x, rb->eof == 1 ? LZMA_FINISH : LZMA_RUN)) {
		default:
		case LZMA_STREAM_END:
		if(rb->eof == 1) { rbdebug("input buffer starved"); rb->eof = 2; }
		case LZMA_OK: break;
	}
});
#endif
#undef _rb_reader

static inline
rb_reader_res_t rb_read_transparent(rbread_t *rb, void *dst, size_t len)
{
	size_t const bytes = fread(dst, sizeof(char), len, rb->fp);
	rb->eof = feof(rb->fp)<<1;		/* input stream starved == input buffer starved for transparent */
	return((rb_reader_res_t){
		.obtained  = bytes,
		.remaining = 0
	});
}

/**
 * @fn rbread
 * @brief gzread compatible
 */
static inline
size_t rbread_bulk(rbread_t *rb, void *dst, size_t len)
{
	/* output buffer must be larger than (or equal to) bulk_size */
	if(rb->eof > 2) { return(0); }	/* if output buffer starved */
	size_t rem = len;
	uint8_t *p = dst;
	rbdebug("len(%zu), rem(%zu), head(%zu), tail(%zu), eof(%u)", len, rem, rb->head, rb->tail, rb->eof);

	/* flush remaining content in buffer */
	size_t const hlen = MIN2(rem, rb->tail - rb->head);
	if(hlen > 0) { memcpy(&p[len - rem], &rb->buf[rb->head], hlen); }
	rem      -= hlen;
	rb->head += hlen;
	rbdebug("len(%zu), rem(%zu), head(%zu), tail(%zu), eof(%u)", len, rem, rb->head, rb->tail, rb->eof);

	/*
	 * EOF state fixup before continue to bulk loop;
	 * we only have input buffer in transparent mode so starvation of the buffer is the end of file.
	 */
	if(rb->head == rb->tail) {
		rbdebug("output buffer starved");
		rb->eof += rb->eof == 2;
	}

	/*
	 * rb->eof == 1 indicates there is only small chunk remaining in the decoding input buffer,
	 * so try get everything remaining in this case.
	 *
	 * otherwise the stream might continue longer than the output buffer,
	 * so check remaining buffer size first and postpone reading if it is smaller than bulk_size / 8.
	 */

	/* if not EOF, try bulk read */
	while(rb->eof < 3 && rem > rb->bulk_size / 8) {
		rb_reader_res_t const r = rb->read(rb, &p[len - rem], rem);
		rem -= r.obtained;
		if(rb->eof == 2 && r.remaining == 0) { rbdebug("output buffer starved"); rb->eof = 3; }
		rbdebug("len(%zu), rem(%zu), head(%zu), tail(%zu), eof(%u)", len, rem, rb->head, rb->tail, rb->eof);
	}
	rbdebug("len(%zu), rem(%zu), head(%zu), tail(%zu), eof(%u)", len, rem, rb->head, rb->tail, rb->eof);
	return(len - rem);
}
static inline
size_t rbread(rbread_t *rb, void *dst, size_t len)
{
	/* if output buffer starved it's the end */
	if(rb->eof > 2) { return(0); }

	size_t rem = len - rbread_bulk(rb, dst, len);
	rbdebug("len(%zu), rem(%zu), head(%zu), tail(%zu), eof(%u)", len, rem, rb->head, rb->tail, rb->eof);

	uint8_t *p = dst;
	while(rb->eof < 3 && rem > 0) {
		rb_reader_res_t const r = rb->read(rb, rb->buf, rb->bulk_size);
		rb->tail = r.obtained;
		rb->head = 0;
		size_t const tlen = MIN2(rem, rb->tail);
		memcpy(&p[len - rem], &rb->buf[rb->head], tlen);

		rem      -= tlen;
		rb->head += tlen;
		if(rb->eof == 2 && rb->head == rb->tail) { rbdebug("output buffer starved"); rb->eof = 3; }
		rbdebug("len(%zu), rem(%zu), head(%zu), tail(%zu), eof(%u)", len, rem, rb->head, rb->tail, rb->eof);
	}
	rbdebug("len(%zu), rem(%zu), head(%zu), tail(%zu), eof(%u)", len, rem, rb->head, rb->tail, rb->eof);
	return(len - rem);
}

/**
 * @fn rbreadline
 * @brief readline; contents are kept inside rbread_t, returns NULL for end of file.
 */
static inline
uint64_t rb_scan_newline(rbread_t *rb, size_t head)
{
	v64i8_t const lv = _set_v64i8('\n');
	for(size_t i = head; i < rb->tail; i += 64) {
		v64i8_t const v = _loadu_v64i8(&rb->buf[i]);
		v64i8_t const x = _eq_v64i8(v, lv);

		uint64_t const mask = ((v64i8_masku_t){ .mask = _mask_v64i8(x) }).all;
		if(_likely(mask == 0)) { continue; }

		/* '\n' found */
		ZCNT_RESULT size_t const idx = _tzcnt_u64(mask);
		rb->buf[i + idx] = '\0';
		return(i + idx);
	}
	return(rb->tail);
}

typedef struct {
	uint8_t const *ptr;
	size_t len;
} rb_readline_t;

static inline
rb_readline_t rbreadline_intl(rbread_t *rb)
{
	if(rb->eof > 2) {
		return((rb_readline_t){
			.ptr = NULL,
			.len = 0
		});
	}

	/* load context variables */
	size_t head = rb->head;
	while((head = rb_scan_newline(rb, head)) == rb->tail) {
		if(rb->eof > 1) { rb->eof++; break; }

		/* not found */
		if(rb->head != 0) {
			memcpy(rb->buf, &rb->buf[rb->head], rb->tail - rb->head);
			rb->tail -= rb->head;
			rb->head  = 0;
		}

		/* expand array if need */
		if(rb->buf_size - rb->tail < rb->bulk_size) {
			rb->buf_size *= 2;
			rb->buf = realloc(rb->buf, rb->buf_size);		/* aligned to ARCH_HUGEPAGE_SIZE */
		}

		/* read more */
		rb_reader_res_t const r = rb->read(rb, &rb->buf[rb->tail], rb->bulk_size);
		rb->tail += r.obtained;
	}

	/* found */
	uint8_t const *ptr = &rb->buf[rb->head];
	size_t const len = head - rb->head;

	rb->head = head + 1;
	rb->buf[head] = '\0';
	return((rb_readline_t){
		.ptr = ptr,
		.len = len
	});
}

static inline
uint8_t const *rbreadline(rbread_t *rb)
{
	rb_readline_t const r = rbreadline_intl(rb);
	return(r.ptr);
}

/**
 * @fn rbopen, rbclose, rbeof
 * @brief open decompressor (reader stream); always single-threaded (seek, flush, ... are not implemented)
 */
static inline
void rb_ptr_destroy(rbread_t *rb)
{
	if(rb->fp) { fclose(rb->fp); }
	free(rb->ibuf);
	free(rb->buf);
	return;
}
static inline
int rbopen_bulk_static(rbread_t *rb, char const *fn, size_t bulk_size)
{
	/* open file */
	rb->fp = (fn == NULL || strcmp(fn, "-") == 0) ? fdopen(fileno(stdin), "r") : fopen(fn, "rb");
	if(rb->fp == NULL) { goto _rbopen_fail; }

	/* calculate buffer size */
	size_t const buffer_size = _roundup(2 * bulk_size, RB_BUF_ALIGN_SIZE);
	rb->buf_size  = buffer_size;
	rb->bulk_size = bulk_size;

	/* read the first chunk */
	uint8_t *buf = aligned_malloc(buffer_size);
	size_t const rlen = fread(buf, sizeof(char), rb->bulk_size, rb->fp);
	if(feof(rb->fp)) { rb->eof = 1; }

	static struct { uint64_t magic, mask; } const x[4] = {
		{ 0x000000088b1fULL, 0x000000ffffffULL },		/* gzip */
		{ 0x000000685a42ULL, 0x000000ffffffULL },		/* bz2 */
		{ 0x005a587a37fdULL, 0xffffffffffffULL },		/* lzma */
		{ 0, 0 }										/* sentinel */
	};

	/* determine file type; FIXME: refactor redundant code */
	#define _rb_init(_sfx, _ctx, _type_t, _body) { \
		rb->read = rb_read_##_sfx; \
		rb->ibuf = buf; \
		rb->buf  = aligned_malloc(buffer_size); \
		rb->s._ctx.next_in   = (_type_t *)buf; \
		rb->s._ctx.avail_in  = rlen; \
		rb->s._ctx.next_out  = (_type_t *)rb->buf; \
		rb->s._ctx.avail_out = 2 * rb->bulk_size; \
		if(rb->buf == NULL || !(_body)) { goto _rbopen_fail; } \
	}
	if(rlen >= sizeof(uint64_t)) {
		uint64_t h = _loadu_u64(buf), i = 0;
		while((h & x[i].mask) != x[i].magic) { i++; }
		switch(i) {
			#ifdef RB_USE_ZLIB
				case 0: _rb_init(gzip, z, uint8_t, ({ inflateInit2(&rb->s.z, 15 + 16) == Z_OK; })); break;
			#else
				case 0: goto _rbopen_fail;
			#endif

			#ifdef RB_USE_BZLIB
				case 1: _rb_init(bz2,  b, char,    ({ BZ2_bzDecompressInit(&rb->s.b, 0, 0) == BZ_OK; })); break;
			#else
				case 1: goto _rbopen_fail;
			#endif

			#ifdef RB_USE_LZMA
				case 2: _rb_init(xz,   x, uint8_t, ({ lzma_stream_decoder(&rb->s.x, UINT64_MAX, LZMA_CONCATENATED) == LZMA_OK; })); break;
			#else
				case 2: goto _rbopen_fail;
			#endif
			default: break;
		}
	}
	#undef _rb_init

	if(rb->read == NULL) {
		rb->read = rb_read_transparent;
		rb->buf  = buf;
		rb->tail = rlen;
	}
	return(0);

_rbopen_fail:;
	rb_ptr_destroy(rb);
	return(1);
}
static inline
rbread_t *rbopen(char const *fn)
{
	rbread_t *rb = calloc(1, sizeof(rbread_t));
	if(rb == NULL) { return(NULL); }
	if(rbopen_bulk_static(rb, fn, RB_BUF_SIZE) != 0) {
		free(rb);
		return(NULL);
	}
	return(rb);
}
static inline
void rbclose_static(rbread_t *rb)
{
	#ifdef RB_USE_ZLIB
		if(rb->read == rb_read_gzip) {
			deflateEnd(&rb->s.z);
			return;
		}
	#endif

	#ifdef RB_USE_BZLIB
		if(rb->read == rb_read_bz2) {
			BZ2_bzDecompressEnd(&rb->s.b);
		}
	#endif

	#ifdef RB_USE_LZMA
		if(rb->read == rb_read_xz) {
			lzma_end(&rb->s.x);
		}
	#endif

	/* transparent */
	rb_ptr_destroy(rb);
	memset(rb, 0, sizeof(rbread_t));
	return;
}
static inline
void rbclose(rbread_t *rb)
{
	rbclose_static(rb);
	free(rb);
	return;
}


#endif
/**
 * end of rbread.h
 */
