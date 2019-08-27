
/**
 * @rbread.h
 * @brief gzread-compatible stream reader
 */

#ifndef _RBREAD_H_INCLUDED
#define _RBREAD_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"

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


#define RB_BUF_SIZE				( 2ULL * 1024 * 1024 )

/**
 * @struct rbread_s
 * @brief context (gzFile)
 */
struct rb_reader_result_s { size_t obtained, remaining; };		/* (%rax, %rdx) pair */
typedef struct rbread_s {
	FILE *fp;
	uint8_t *ibuf, *buf;	/* input buffer and output buffer; ibuf is unused for transparent mode */
	size_t head, len;
	size_t bulk_size;
	uint8_t eof, _pad[7];
	struct rb_reader_result_s (*read)(struct rbread_s *, void *, size_t);
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

/**
 * @fn rb_read_gzip, rb_read_bz2, rb_read_xz, rb_read_transparent
 * @brief block fetcher
 */
#define _rb_reader(_sfx, _ctx, _type_t, _body) \
	static inline struct rb_reader_result_s rb_read_##_sfx(rbread_t *rb, void *dst, size_t len) \
	{ \
		if(rb->s._ctx.avail_in < sizeof(uint64_t) && rb->eof == 0) { \
			/* copy to the head then fill */ \
			uint64_t const _remaining = _loadu_u64(rb->s._ctx.next_in); \
			_storeu_u64(rb->ibuf, _remaining); \
			rb->s._ctx.next_in = (_type_t *)&rb->ibuf[rb->s._ctx.avail_in]; \
			rb->s._ctx.avail_in += fread((void *)rb->s._ctx.next_in, sizeof(char), rb->bulk_size, rb->fp); \
			if(feof(rb->fp)) { rb->eof = 1; } \
		} \
		rb->s._ctx.next_out = (_type_t *)dst; \
		rb->s._ctx.avail_out = len; \
		_body; \
		return((struct rb_reader_result_s){ \
			.obtained = rb->eof == 2 ? 0 : len - rb->s._ctx.avail_out, \
			.remaining = rb->s._ctx.avail_in \
		}); \
	}

#ifdef RB_USE_ZLIB
_rb_reader(gzip, z, uint8_t, {
	switch(inflate(&rb->s.z, Z_NO_FLUSH)) {
		default: rb->eof = 2;
		case Z_STREAM_END: inflateReset(&rb->s.z);
		case Z_OK: break;
	}
});
#endif

#ifdef RB_USE_BZLIB
_rb_reader(bz2, b, char, {
	switch(BZ2_bzDecompress(&rb->s.b)) {
		default: rb->eof = 2;
		case BZ_STREAM_END: case BZ_OK: break;
	}
});
#endif

#ifdef RB_USE_LZMA
_rb_reader(xz, x, uint8_t, {
	switch(lzma_code(&rb->s.x, rb->eof == 1 ? LZMA_FINISH : LZMA_RUN)) {
		default: rb->eof = 2;
		case LZMA_STREAM_END: case LZMA_OK: break;
	}
});
#endif
#undef _rb_reader

static inline
struct rb_reader_result_s rb_read_transparent(rbread_t *rb, void *dst, size_t len)
{
	size_t const bytes = fread(dst, sizeof(char), len, rb->fp);
	rb->eof = feof(rb->fp);
	return((struct rb_reader_result_s){
		.obtained = bytes,
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
	if(rb->eof > 1) { return(0); }
	size_t rem = len;
	uint8_t *p = dst;

	/* flush remaining content in buffer */
	size_t const hlen = MIN2(rem, rb->len - rb->head);
	if(hlen > 0) { memcpy(&p[len - rem], &rb->buf[rb->head], hlen); }
	rem -= hlen; rb->head += hlen;

	/*
	 * EOF state fixup before continue to bulk loop;
	 * we only have input buffer in transparent mode so starvation of the buffer is the end of file.
	 */
	if(rb->ibuf == NULL && rb->head == rb->len) {
		rb->eof = rb->eof ? 2 : 0;
	}

	/*
	 * rb->eof == 1 indicates there is only small chunk remaining in the decoding buffer,
	 * so try get everything in this case.
	 * otherwise the stream might continue longer than the output buffer,
	 * so check remaining buffer size first and postpone reading if it is smaller than bulk_size / 8.
	 */

	/* if not EOF, try bulk read */
	while(rb->eof < 2 && rem > rb->bulk_size / 8) {
		struct rb_reader_result_s const r = rb->read(rb, &p[len - rem], rem);
		rem -= r.obtained;
		if(rb->eof == 1 && r.remaining == 0) { rb->eof = 2; }
	}
	return(len - rem);
}
static inline
size_t rbread(rbread_t *rb, void *dst, size_t len)
{
	if(rb->eof > 1) { return(0); }
	size_t rem = len - rbread_bulk(rb, dst, len);
	uint8_t *p = dst;
	while(rb->eof < 2 && rem > 0) {
		struct rb_reader_result_s const r = rb->read(rb, rb->buf, rb->bulk_size);
		rb->len = r.obtained;
		rb->head = 0;
		size_t tlen = MIN2(rem, rb->len);
		memcpy(&p[len - rem], &rb->buf[rb->head], tlen);
		rem -= tlen; rb->head += tlen;
		if(rb->eof == 1 && rb->head == rb->len) { rb->eof = 2; }
	}
	return(len - rem);
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
	rb->bulk_size = bulk_size;

	/* open file */
	rb->fp = (fn == NULL || strcmp(fn, "-") == 0) ? fdopen(fileno(stdin), "r") : fopen(fn, "rb");
	if(rb->fp == NULL) { goto _rbopen_fail; }

	/* calculate buffer size */
	size_t const buffer_size = _roundup(2 * rb->bulk_size, ARCH_HUGEPAGE_SIZE);

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
		rb->buf = aligned_malloc(buffer_size); \
		rb->s._ctx.next_in = (_type_t *)buf; \
		rb->s._ctx.avail_in = rlen; \
		rb->s._ctx.next_out = (_type_t *)rb->buf; \
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
		rb->buf = buf;
		rb->len = rlen;
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
	return;
}
static inline
void rbclose(rbread_t *rb)
{
	rbclose_static(rb);
	free(rb);
	return;
}

#define rbeof(_rb)				( (_rb)->eof >= 2 )


#endif
/**
 * end of rbread.h
 */
