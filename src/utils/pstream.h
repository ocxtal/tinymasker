
/**
 * @file pstream.h
 * multithreading (queue and worker thread)
 */

#ifndef _PSTREAM_H_INCLUDED
#define _PSTREAM_H_INCLUDED


/* use pthread_setaffinity_np to pin threads to cpu cores */
// #define PT_SET_AFFINITY


/* transparent mode for debugging */
#define PG_TRANSPARENT


/* use pithy for faster compression / decompression than zlib */
// #define PG_USE_PITHY

/* include global header *before* we include individual dependencies */
#include "common.h"
#include <pthread.h>
#include <unistd.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/syscall.h>

#ifdef PG_TRANSPARENT
#  define PG_MAGIC_MINOR			( 0 )
#else
#  ifdef PG_USE_PITHY
#    include "pithy.h"
#    define PG_COMP_LEVEL			( 5 )
#    define PG_MAGIC_MINOR			( 2 )
#  else
#    include <zlib.h>
#    define PG_COMP_LEVEL			( 1 )
#    define PG_MAGIC_MINOR			( 1 )
#  endif
#endif

/**
 * @type pt_source_t, pt_worker_t, pt_drain_t
 * @brief callback functions types for multithreaded stream
 */
typedef void *(*pt_source_t)(uint32_t tid, void *opaque);
typedef void *(*pt_worker_t)(uint32_t tid, void *opaque, void *item);
typedef void (*pt_drain_t)(uint32_t tid, void *opaque, void *item);

/**
 * @struct pt_q_s
 * @brief lock-based queue context
 */
typedef struct pt_q_s {
	size_t lock, head, tail, size;
	void **elems;
	size_t wait_cnt;
	size_t _pad1[2];
} pt_q_t;
_static_assert(sizeof(pt_q_t) == 64);		/* to occupy a single cache line */

/**
 * @struct pt_thread_s
 * @brief thread-local worker context
 */
typedef struct pt_thread_s {
	pthread_t th;
	size_t tid, wait_cnt;
	pt_q_t *in, *out;
	pt_worker_t wfp;
	void *wopaque;
} pt_thread_t;

/**
 * @struct pt_s
 * @brief parallel task processor context
 */
typedef struct pt_s {
	pt_q_t in, out;
	size_t nth;
	pt_thread_t c[];	/* [0] is reserved for master */
} pt_t;
#define pt_nth(_pt)				( (_pt)->nth )
#define PT_EMPTY				( (void *)(UINT64_MAX) )
#define PT_EXIT					( (void *)(UINT64_MAX - 1) )
#define PT_DEFAULT_INTERVAL		( 512 * 1024 )

/**
 * @fn pt_enq
 * @brief enqueue; concurrent queue is better while currently lock-based for simplicity
 */
static _force_inline
size_t pt_enq(pt_q_t *q, size_t tid, void *elem)
{
	size_t z, ret = UINT64_MAX;

	/* lock by thread id */
	do { z = UINT32_MAX; } while(!cas(&q->lock, &z, tid));

	/* push element to queue */
	size_t head = q->head, tail = q->tail, size = q->size;
	if(((head + 1) % size) != tail) {
		q->elems[head] = elem;
		q->head = (head + 1) % size;
		ret = 0;
	}

	/* release */
	do { z = tid; } while(!cas(&q->lock, &z, UINT32_MAX));
	return(ret);
}
static _force_inline
size_t pt_enq_retry(pt_q_t *q, size_t tid, void *elem, uint64_t nsec)
{
	struct timespec tv = { .tv_nsec = nsec };
	while(pt_enq(q, tid, elem) != 0) { q->wait_cnt++; nanosleep(&tv, NULL); }
	return(0);
}

/**
 * @fn pt_deq
 */
static _force_inline
void *pt_deq(pt_q_t *q, size_t tid)
{
	void *elem = PT_EMPTY;
	size_t z;

	/* lock by thread id */
	do { z = UINT32_MAX; } while(!cas(&q->lock, &z, tid));

	/* get element from queue */
	size_t head = q->head, tail = q->tail, size = q->size;
	if(head != tail) {
		elem = q->elems[tail]; q->tail = (tail + 1) % size;
	}

	/* release */
	do { z = tid; } while(!cas(&q->lock, &z, UINT32_MAX));
	return(elem);
}


static _force_inline
void pt_setaffinity(uint32_t id)
{
	#if defined(__linux__) && defined(PT_SET_AFFINITY)
		cpu_set_t cpu;
		CPU_ZERO(&cpu);
		CPU_SET(id, &cpu);

		pid_t tid = syscall(SYS_gettid);
		sched_setaffinity(tid, sizeof(cpu_set_t), &cpu);
	#else
		_unused(id);
	#endif
	return;
}

/**
 * @fn pt_dispatch
 * @brief per-thread function dispatcher, with ping-pong prefetching
 */
static _force_inline
void *pt_dispatch(void *s)
{
	pt_thread_t *c = (pt_thread_t *)s;
	pt_setaffinity(c->tid - 1);

	uint64_t const intv = PT_DEFAULT_INTERVAL;
	struct timespec tv = { .tv_nsec = intv };

	void *ping = PT_EMPTY, *pong = PT_EMPTY;
	while(1) {
		ping = pt_deq(c->in, c->tid);		/* prefetch */
		if(ping == PT_EMPTY && pong == PT_EMPTY) {
			c->wait_cnt++; nanosleep(&tv, NULL);			/* no task is available, sleep for a while */
		}
		if(pong != PT_EMPTY) {
			pt_enq_retry(c->out, c->tid, c->wfp(c->tid, (void *)c->wopaque, pong), intv);
		}
		if(ping == PT_EXIT) { break; }		/* terminate thread */

		pong = pt_deq(c->in, c->tid);		/* prefetch */
		if(ping == PT_EMPTY && pong == PT_EMPTY) {
			c->wait_cnt++; nanosleep(&tv, NULL);			/* no task is available, sleep for a while */
		}
		if(ping != PT_EMPTY) {
			pt_enq_retry(c->out, c->tid, c->wfp(c->tid, (void *)c->wopaque, ping), intv);
		}
		if(pong == PT_EXIT) { break; }		/* terminate thread */
	}
	return(NULL);
}

/**
 * @fn pt_destroy_static
 */
static _force_inline
void pt_destroy_static(pt_t *pt)
{
	void *status;
	if(pt == NULL) { return; }

	/* send termination signal */
	for(size_t i = 1; i < pt->nth; i++) {
		pt_enq_retry(pt->c->in, pt->c->tid, PT_EXIT, PT_DEFAULT_INTERVAL);
	}

	/* wait for the threads terminate */
	for(size_t i = 1; i < pt->nth; i++) { pthread_join(pt->c[i].th, &status); }

	/* clear queues */
	while(pt_deq(pt->c->in, pt->c->tid) != PT_EMPTY) {}
	while(pt_deq(pt->c->out, pt->c->tid) != PT_EMPTY) {}

	/* clear queues and objects */
	free(pt->in.elems);
	free(pt->out.elems);
	return;
}
#define pt_destroy(_pt)				{ pt_destroy_static(_pt); free(_pt); }

/**
 * @fn pt_init_static
 */
static _force_inline
void pt_init_static(pt_t *pt, uint32_t nth)
{
	/* init object */
	pt->nth = (nth == 0) ? 1 : nth;

	/* init queues (note: #elems can be larger for better performance?) */
	size_t const size = 16 * nth;
	pt->in = (pt_q_t){
		.lock = UINT32_MAX,
		.elems = calloc(size, sizeof(void *)),
		.size = size
	};
	pt->out = (pt_q_t){
		.lock = UINT32_MAX,
		.elems = calloc(size, sizeof(void *)),
		.size = size
	};

	pt_setaffinity(nth - 1);

	/* init children info, create children */
	for(size_t i = 0; i < nth; i++) {
		pt->c[i].tid = i;
		pt->c[i].in = &pt->in;
		pt->c[i].out = &pt->out;

		if(i == 0) { continue; }			/* do not create child thread for tid 0 */
		pthread_create(&pt->c[i].th, NULL, pt_dispatch, (void *)&pt->c[i]);
	}
	return;
}
#define pt_init_size(_nth)		( sizeof(pt_t) + (_nth) * sizeof(pt_thread_t) )
#define pt_init(_nth) ({ \
	size_t nth = ((_nth) == 0) ? 1 : (_nth); \
	pt_t *pt = calloc(1, pt_init_size(_nth)); \
	pt_init_static(pt, nth); \
	pt; \
})

/**
 * @fn pt_set_worker
 * @brief update worker function and zument pointers (FIXME: might cause race)
 */
static _force_inline
int pt_set_worker(pt_t *pt, void *opaque, pt_worker_t wfp)
{
	void *item = NULL;

	/* fails when unprocessed object exists in the queue */
	if((item = pt_deq(&pt->in, 0)) != PT_EMPTY) {
		pt_enq_retry(&pt->in, 0, item, PT_DEFAULT_INTERVAL);
		return(-1);
	}

	/* update pointers */
	for(size_t i = 0; i < pt->nth; i++) {
		pt->c[i].wfp = wfp;
		pt->c[i].wopaque = opaque;
	}

	/* syncronize */
	fence();
	return(0);
}

/**
 * @fn pt_stream
 * @brief multithreaded stream, source and drain are always called in the parent thread
 */
static _force_inline
int pt_stream(pt_t *pt, void *opaque, pt_source_t sfp, pt_worker_t wfp, pt_drain_t dfp)
{
	if(pt_set_worker(pt, opaque, wfp)) { return(-1); }

	/* keep balancer between [lb, ub) */
	size_t const lb = 3 * pt->nth, ub = 4 * pt->nth;
	size_t bal = 0;
	void *it;

	while((it = sfp(0, opaque)) != NULL) {
		// if(pt_enq(&pt->in, 0, it) != 0) { dfp(0, opaque, wfp(0, opaque, it)); continue; }	/* queue full, process locally */
		pt_enq(&pt->in, 0, it);
		if(++bal < ub) { continue; }				/* queue not full */

		/* bal <= lb indicates input queue might be starving */
		while(bal > lb) {
			/* flush to drain (note: while loop is better?) */
			if((it = pt_deq(&pt->out, 0)) != PT_EMPTY) {
				bal--; dfp(0, opaque, it); continue;
			}

			/* output queue empty; process one in the master (parent) thread */
			if((it = pt_deq(&pt->in, 0)) != PT_EMPTY) { bal--; dfp(0, opaque, wfp(0, opaque, it)); }
		}
	}

	/* source depleted, process remainings */
	while((it = pt_deq(&pt->in, 0)) != PT_EMPTY) { pt_enq_retry(&pt->out, 0, wfp(0, opaque, it), PT_DEFAULT_INTERVAL); }

	/* flush results */
	while(bal > 0) {
		if((it = pt_deq(&pt->out, 0)) != PT_EMPTY) {
			bal--; dfp(0, opaque, it);
		} else {
			struct timespec tv = { .tv_nsec = 2 * 1024 * 1024 };
			nanosleep(&tv, NULL);
		}
	}
	return(0);
}

/**
 * @fn pt_parallel
 */
typedef struct {
	void *opaque;
	pt_worker_t wfp;
	void *items;
	size_t item_cnt, item_size, curr_idx;
} pt_parallel_t;

static void *pt_parallel_source(uint32_t tid, pt_parallel_t *w)
{
	_unused(tid);
	if(w->curr_idx >= w->item_cnt) { return(NULL); }
	size_t idx = w->curr_idx++;
	return(_add_offset(w->items, w->item_size * idx));
}
static void *pt_parallel_worker(uint32_t tid, pt_parallel_t *w, void *item)
{
	w->wfp(tid, w->opaque, item);
	return(NULL);
}
static void pt_parallel_drain(uint32_t tid, pt_parallel_t *w, void *item)
{
	_unused(tid);
	_unused(w);
	_unused(item);
	return;				/* nothing to do */
}

static _force_inline
int pt_parallel(pt_t *pt, void *opaque, pt_worker_t wfp, void *items, size_t item_cnt, size_t item_size)
{
	pt_parallel_t w = {
		.opaque = opaque,
		.wfp = wfp,
		.items = items,
		.item_cnt = item_cnt,
		.item_size = item_size,
		.curr_idx = 0
	};
	return(pt_stream(pt, &w,
		(pt_source_t)pt_parallel_source,
		(pt_worker_t)pt_parallel_worker,
		(pt_drain_t)pt_parallel_drain
	));
}


#if 0
static _force_inline
int pt_parallel(pt_t *pt, void *opaque, pt_worker_t wfp, void *items, size_t item_cnt, size_t item_size)
{
	/* fails when unprocessed element exists */
	if(pt_set_worker(pt, opaque, wfp)) { return(-1); }

	/* push items */
	for(size_t i = 1; i < pt->nth; i++) {
		pt_enq_retry(&pt->in, 0, (void *)i, PT_DEFAULT_INTERVAL);
	}
	debug("pushed items");

	/* process the first one in the master (parent) thread */
	wfp(0, opaque, (void *)0);
	debug("finished master");

	/* wait for the children done */
	for(size_t i = 1; i < pt->nth; i++) {
		while(pt_deq(&pt->out, 0) == PT_EMPTY) {
			struct timespec tv = { .tv_nsec = 512 * 1024 };
			nanosleep(&tv, NULL);	/* wait for a while */
			/* sched_yield(); */
		}
		debug("joined i(%lu)", i);
	}
	return 0;
}
#endif

#if 0
/* unittests */
static void *pt_unittest_source(uint32_t tid, void *arg)
{
	uint64_t *s = (uint64_t*)arg;
	if(*s >= 1024) return NULL;
	uint64_t *p = malloc(sizeof(uint64_t));
	*p = *s; *s = *s + 1;
	return p;
}
static void *pt_unittest_worker(uint32_t tid, void *arg, void *item)
{
	uint64_t *p = (uint64_t *)item, *a = (uint64_t *)arg, i = *a;
	*p += i;
	return p;
}
static void pt_unittest_drain(uint32_t tid, void *arg, void *item)
{
	uint64_t *d = (uint64_t*)arg, *p = (uint64_t*)item, i = *p;
	*d += i;
	free(item);
}

unittest( .name = "pt.single" ) {
	pt_t *pt = pt_init(1);
	ut_assert(pt != NULL);

	uint64_t icnt = 0, ocnt = 0, inc = 1, *arr[1] = { &inc };
	pt_stream(pt, arr, pt_unittest_source, pt_unittest_worker, pt_unittest_drain);
	ut_assert(icnt == 1024, "icnt(%lu)", icnt);
	ut_assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0 };	//, *sp[1] = { &s[0] };
	pt_parallel(pt, arr, pt_unittest_worker);
	ut_assert(s[0] == 1, "d[0](%lu)", s[0]);
	pt_destroy(pt);
}

unittest( .name = "pt.multi" ) {
	pt_t *pt = pt_init(4);
	ut_assert(pt != NULL);

	uint64_t icnt = 0, ocnt = 0, inc = 1, *arr[4] = { &inc, &inc, &inc, &inc };
	pt_parallel(pt, arr, pt_unittest_worker);
	ut_assert(icnt == 1024, "icnt(%lu)", icnt);
	ut_assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0,1,2,3 };	//, *sp[4] = { &s[0],&s[1],&s[2],&s[3] };
	pt_parallel(pt, arr, pt_unittest_worker);
	for(uint64_t i = 0; i < 4; i++) {
		ut_assert(s[i] == (i+1), "i(%lu), d[i](%lu)", i, s[i]);
	}
	pt_destroy(pt);
}

#endif

/* stdio stream with multithreaded compression / decompression */

#define PG_BLOCK_SIZE				( 1024 * 1024 )				/* 1M; half of hugepage */
#define PG_MAGIC					( 0x30304750 + (PG_MAGIC_MINOR<<24) )			/* "PG00" */
#define PG_MAGIC_SIZE				( 4 )
#define pg_eof(_pg)					( (_pg)->eof == 2 ? 1 : 0 )

/**
 * @struct pg_header_t
 */
typedef struct {
	uint32_t magic, len;			/* each block is no longer than gigabyte */
} pg_header_t;
_static_assert(sizeof(pg_header_t) == 8);

/**
 * @struct pg_block_t
 * @brief compression / decompression unit block
 */
typedef struct {
	uint32_t head, len;				/* head pointer (index) and block length */
	uint32_t id;					/* block id */
	uint8_t raw, flush, _pad[2];	/* raw: 1 if compressed, flush: 1 if needed to dump */
	uint8_t buf[];
} pg_block_t;

/**
 * @struct pg_t
 * @brief context
 */
typedef struct {
	size_t in, out;
} pg_stat_t;

typedef struct {
	size_t id;
	pg_block_t *block;
} pg_hqbin_t;

typedef struct {
	FILE *fp;						/* must be transparent */
	pt_t *pt;

	/* working buffers */
	pg_block_t *s;					/* current */
	uint32_t ub, lb, bal, eof, is_out, nth;		/* eof: 0 for input (FILE *) remaining, 1 for FILE * depleted but remaining in buffer, 2 for completed, >= 3 for error */

	/* bucket counter and sorter */
	struct {
		size_t in, out;
	} cnt;
	kvechq_t(pg_hqbin_t) hq;

	/* constant */
	size_t block_size;

	/* stat */
	pg_stat_t bytes;
} pg_t;
#define incq_comp(a, b)		( (int64_t)(a).id - (int64_t)(b).id )

/**
 * @fn pg_deflate
 * @brief compress block
 */
#ifdef PG_TRANSPARENT
static _force_inline
size_t pg_deflate_calc_size(uint8_t const *in, size_t in_size, size_t block_size)
{
	_unused(in);
	_unused(block_size);
	return(in_size);
}
static
size_t pg_deflate_core(uint8_t const *in, size_t in_size, uint8_t *out, size_t out_size)
{
	_unused(out_size);
	memcpy(out, in, in_size);
	return(in_size);
}

#elif defined(PG_USE_PITHY)
static _force_inline
size_t pg_deflate_calc_size(uint8_t const *in, size_t in_size, size_t block_size)
{
	_unused(in);
	_unused(in_size);
	return(1.2 * block_size);
}
static _force_inline
size_t pg_deflate_core(uint8_t const *in, size_t in_size, uint8_t *out, size_t out_size)
{
	/* compress (pithy) */
	return(pithy_Compress(
		(char const *)in, in_size,
		(char *)out, out_size,
		PG_COMP_LEVEL		/* compression level; 0 <= x <= 9 */
	));
}

#else
static _force_inline
size_t pg_deflate_calc_size(uint8_t const *in, size_t in_size, size_t block_size)
{
	_unused(in);
	_unused(in_size);
	return(1.2 * block_size);
}
static _force_inline
size_t pg_deflate_core(uint8_t const *in, size_t in_size, uint8_t *out, size_t out_size)
{
	/* compress (deflate) */
	z_stream zs = {
		.next_in = (uint8_t *)in, .avail_in = in_size,
		.next_out = out,          .avail_out = out_size
	};
	if(deflateInit2(&zs, PG_COMP_LEVEL, Z_DEFLATED, 15, 8, Z_DEFAULT_STRATEGY) != Z_OK
	|| deflate(&zs, Z_FINISH) != Z_STREAM_END
	|| deflateEnd(&zs) != Z_OK) {
		return(0);
	}

	return(out_size - zs.avail_out);
}

#endif

static
pg_block_t *pg_deflate(pg_block_t *in, size_t block_size)
{
	size_t buf_size = pg_deflate_calc_size(in->buf, in->len, block_size);
	size_t malloc_size = _roundup(sizeof(pg_block_t) + buf_size, ARCH_HUGEPAGE_SIZE);
	pg_block_t *out = aligned_malloc(malloc_size);

	/* dispatch */
	size_t out_size = pg_deflate_core(in->buf, in->len, out->buf, buf_size);
	if(out_size == 0) {
		free(out);
		return(NULL);
	}

	/* set metadata */
	out->head = 0;
	out->len = out_size;
	out->id = in->id;
	out->raw = 0;
	out->flush = 1;

	/* cleanup input block */
	free(in);
	return(out);
}

/**
 * @fn pg_inflate
 * @brief decompress block
 */
#ifdef PG_TRANSPARENT
static _force_inline
size_t pg_inflate_calc_size(uint8_t const *in, size_t in_size, size_t block_size)
{
	_unused(in);
	_unused(block_size);
	return(in_size);
}
static _force_inline
size_t pg_inflate_core(uint8_t const *in, size_t in_size, uint8_t *out, size_t out_size)
{
	_unused(out_size);
	memcpy(out, in, in_size);
	return(in_size);
}

#elif defined(PG_USE_PITHY)
static _force_inline
size_t pg_inflate_calc_size(uint8_t const *in, size_t in_size, size_t block_size)
{
	_unused(block_size);

	size_t buf_size = 0;
	pithy_GetDecompressedLength((char const *)in, in_size, &buf_size);
	return(buf_size + 64);
}
static _force_inline
size_t pg_inflate_core(uint8_t const *in, size_t in_size, uint8_t *out, size_t out_size)
{
	int ret = pithy_Decompress((char const *)in, in_size, (char *)out, out_size);
	if(!ret) { message(stderr, "failed decompress, in(%p, %zu), out(%p, %zu)", in, in_size, out, out_size); }
	return(ret ? out_size : 0);
}

#else
static _force_inline
size_t pg_inflate_calc_size(uint8_t const *in, size_t in_size, size_t block_size)
{
	_unused(in);
	_unused(in_size);

	return(2 * block_size);
}
static _force_inline
size_t pg_inflate_core(uint8_t const *in, size_t in_size, uint8_t *out, size_t out_size)
{
	/* inflate */
	z_stream zs = {
		.next_in = (uint8_t *)in, .avail_in = in_size,
		.next_out = out,          .avail_out = out_size
	};
	if(inflateInit2(&zs, 15) != Z_OK
	|| inflate(&zs, Z_FINISH) != Z_STREAM_END
	|| inflateEnd(&zs) != Z_OK) {
		return(0);
	}
	return(out_size - zs.avail_out);
}

#endif

static
pg_block_t *pg_inflate(pg_block_t *in, size_t block_size)
{
	size_t buf_size = pg_inflate_calc_size(in->buf, in->len, block_size);
	size_t malloc_size = _roundup(sizeof(pg_block_t) + buf_size, ARCH_HUGEPAGE_SIZE);
	pg_block_t *out = aligned_malloc(malloc_size);

	/* dispatch */
	size_t out_size = pg_inflate_core(in->buf, in->len, out->buf, buf_size);
	if(out_size == 0) {
		free(out);
		return(NULL);
	}

	/* save metadata */
	out->head = 0;
	out->len = out_size;
	out->id = in->id;
	out->raw = 1;
	out->flush = 0;

	/* cleanup input block */
	free(in);
	return(out);
}

/**
 * @fn pg_worker
 * @brief thread-local worker
 */
static
void *pg_worker(uint32_t tid, void *opaque, void *item)
{
	_unused(tid);
	pg_t *pg = (pg_t *)opaque;
	pg_block_t *s = (pg_block_t *)item;
	if(s == NULL || s->len == 0) { return(s); }
	return((s->raw ? pg_deflate : pg_inflate)(s, pg->block_size));
}

/**
 * @fn pg_read_block
 * @brief read compressed block from input stream
 */
static _force_inline
pg_block_t *pg_read_block(pg_t *pg)
{
	size_t malloc_size = _roundup(sizeof(pg_block_t) + pg->block_size, ARCH_HUGEPAGE_SIZE);
	pg_block_t *s = aligned_malloc(malloc_size);

	/* read block; we expect header is always valid; check magic */
	pg_header_t header = { 0, 0 };
	if(fread(&header, sizeof(pg_header_t), 1, pg->fp) != 1) { goto _fail; }
	if(header.magic != PG_MAGIC || header.len == 0) { goto _fail; }

	/* 0xffffffff as terminator */
	if(header.len == 0xffffffff) {
		pg->eof = MAX2(pg->eof, 1);
		free(s);
		return(NULL);
	}

	/* batch read block */
	if(fread(s->buf, sizeof(uint8_t), header.len, pg->fp) != header.len) { goto _fail; }

	/* set metadata */
	s->head = 0;
	s->len = header.len;
	s->id = pg->cnt.in++;
	s->raw = 0;
	s->flush = 0;

	/* update accumulator */
	pg->bytes.in += header.len;
	return(s);
_fail:
	free(s);
	pg->eof = MAX2(pg->eof, 3);
	return(NULL);
}
static _force_inline
pg_block_t *pg_record_read_block(pg_t *pg, pg_block_t *s)
{
	pg->s = s;
	if(s != NULL) { pg->bytes.out += s->len; }
	return(s);
}

/**
 * @fn pg_write_block
 * @brief write compressed block to output stream
 */
static _force_inline
pg_block_t *pg_record_write_block(pg_t *pg, pg_block_t *s)
{
	if(s != NULL) { pg->bytes.in += s->len; }
	return(s);
}
static _force_inline
void pg_write_block(pg_t *pg, pg_block_t *s)
{
	if(s->len == 0) return;
	pg->bytes.out += s->len;
	// message(stderr, "id(%u), len(%u)", s->id, s->len);

	/* dump header then block */
	pg_header_t header = {
		.magic = PG_MAGIC,
		.len = s->len
	};

	fwrite(&header, sizeof(pg_header_t), 1, pg->fp);
	fwrite(s->buf, sizeof(uint8_t), s->len, pg->fp);
	free(s);
	return;
}

/**
 * @fn pg_init_static
 * @brief initialize pg stream with fp
 */
static _force_inline
void pg_init_static(pg_t *pg, FILE *fp, pt_t *pt)
{
	/* create context */
	*pg = (pg_t){
		.fp = fp, .pt = pt,

		/* constants */
		.block_size = PG_BLOCK_SIZE,
		.nth = pt_nth(pt),

		/* balancer and states */
		.lb = pt_nth(pt),
		.ub = 3 * pt_nth(pt),
		.bal = 0,
		.cnt = { .in = 0, .out = 0 },
		.eof = 0,

		/* stats */
		.bytes = { .in = 0, .out = 0 }
	};
	kv_hq_init(pg->hq);

	/* init worker args */
	pt_set_worker(pg->pt, pg, pg_worker);
	return;
}
#define pg_init(_fp, _pt)		({ pg_t * pg = calloc(1, sizeof(pg_t)); if(pg) { pg_init_static(pg, _fp, _pt); } pg; })

static _force_inline
uint64_t pg_is_outbound(pg_block_t const *s)
{
	return(s != NULL && s->flush != 0);
}

/**
 * @fn pg_freeze
 * @brief clear ptask queues
 */
static _force_inline
void pg_freeze(pg_t *pg)
{
	pg_block_t *s = pg->s, *t;

	/* save direction for pg_destroy */
	if(s != NULL) {
		pg->is_out = pg_is_outbound(s);
	}

	/* process current working block */
	if(pg_is_outbound(s) && s->head != 0) {
		s->len = s->head;

		/* record last block before push */
		s = pg_record_write_block(pg, s);

		/* push to the queue */
		if(pg->nth == 1) {
			pg_write_block(pg, pg_deflate(s, pg->block_size));
		} else {
			pg->bal++;
			pt_enq_retry(&pg->pt->in, 0, s, PT_DEFAULT_INTERVAL);
		}
		pg->s = NULL;
	}

	/* process remainings */
	while(pg->bal > 0) {
		if((t = pt_deq(&pg->pt->out, 0)) == PT_EMPTY) {
			/* wait for a while(note: nanosleep is better?) */
			sched_yield(); continue;
		}
		pg->bal--;
		kv_hq_push(pg_hqbin_t, incq_comp, pg->hq, ((pg_hqbin_t){ .id = t->id, .block = t }));
	}
	return;
}

/**
 * @fn pg_flush
 */
static _force_inline
void pg_flush(pg_t *pg)
{
	pg_freeze(pg);		/* flush input buffer */
	free(pg->s);		/* free if inbound */

	/* flush heapqueue */
	while(kv_hq_cnt(pg->hq) > 0) {
		pg->cnt.out++;
		pg_block_t *t = kv_hq_pop(pg_hqbin_t, incq_comp, pg->hq).block;
		if(pg_is_outbound(t)) { pg_write_block(pg, t); } else { free(t); }
	}
	return;
}

/**
 * @fn pg_destroy_static
 */
static _force_inline
pg_stat_t pg_destroy_static(pg_t *pg)
{
	if(pg == NULL) { return((pg_stat_t){ 0, 0 }); }
	pg_flush(pg);

	/* write terminator if (the last block is) outbound */
	if(pg->is_out) {
		pg_header_t header = {
			.magic = PG_MAGIC,
			.len = 0xffffffff
		};
		fwrite(&header, sizeof(pg_header_t), 1, pg->fp);
	}

	/* cleanup contexts */
	kv_hq_destroy(pg->hq);
	return(pg->bytes);
}
#define pg_destroy(_pg)				({ pg_stat_t stat = pg_destroy_static(_pg); free(_pg); stat; })

/**
 * @fn pgread
 */
static _force_inline
pg_block_t *pg_read_single(pg_t *pg)
{
	/* single-threaded */
	pg_block_t *t;
	if((t = pg_read_block(pg)) == NULL) {
		pg->eof = MAX2(pg->eof, 2);		/* not an error if eof == 1, otherwise something was wrong in read_block */
		return(NULL);
	}
	return(pg_inflate(t, pg->block_size));
}
static _force_inline
pg_block_t *pg_read_multi(pg_t *pg)
{
	/* multithreaded; read compressed blocks and push them to queue */
	pg_block_t *t;
	while(kv_hq_cnt(pg->hq) < pg->ub && !pg->eof && pg->bal < pg->ub) {
		if((t = pg_read_block(pg)) == NULL) { break; }
		pg->bal++;
		pt_enq_retry(&pg->pt->in, 0, t, PT_DEFAULT_INTERVAL);
	}

	/* check if input depleted */
	if(pg->cnt.out >= pg->cnt.in) { pg->eof = MAX2(pg->eof, 2); return(NULL); }

	/* fetch inflated blocks and push heapqueue to sort */
	while(kv_hq_cnt(pg->hq) == 0 || kv_hq_head(pg->hq).id > pg->cnt.out) {
		while((t = pt_deq(&pg->pt->out, 0)) != PT_EMPTY) {
			pg->bal--;
			kv_hq_push(pg_hqbin_t, incq_comp, pg->hq, ((pg_hqbin_t){ .id = t->id, .block = t }));
		}
		if(kv_hq_cnt(pg->hq) >= 1 && kv_hq_head(pg->hq).id == pg->cnt.out) { break; }
		sched_yield();
	}

	/* block is available! */
	pg->cnt.out++;
	return(kv_hq_pop(pg_hqbin_t, incq_comp, pg->hq).block);
}
static _force_inline
size_t pgread(pg_t *pg, void *dst, size_t len)
{
	size_t rem = len;
	pg_block_t *s = pg->s;
	if(pg->eof > 1) return(0);
	if(pg->pt->c->wfp != pg_worker && pt_set_worker(pg->pt, pg, pg_worker)) {
		return(0);
	}

	static pg_block_t *(*const read_block[2])(pg_t *) = { pg_read_single, pg_read_multi };
	while(rem > 0) {
		/* check and prepare a valid inflated block */
		while(s == NULL || s->head == s->len) {
			free(s);
			s = pg_record_read_block(pg, read_block[pg->nth > 1](pg));
			if(pg->eof <= 1 && s == NULL) { pg->eof = MAX2(pg->eof, 3); }
			if(pg->eof > 1) {
				debugblock({ if(pg->eof > 2) { trap(); } })
				return(len - rem);
			}
		}

		/* copy to dst buffer */
		uint64_t adv = MIN2(rem, s->len - s->head);
		memcpy(_add_offset(dst, len - rem), s->buf + s->head, adv);
		rem -= adv; s->head += adv;
	}
	return(len);
}

/**
 * @fn pgwrite
 */
static _force_inline
void pg_write_single(pg_t *pg, pg_block_t *s)
{
	/* single-threaded */
	pg_write_block(pg, pg_deflate(s, pg->block_size));
	return;
}
static _force_inline
void pg_write_multi(pg_t *pg, pg_block_t *s)
{
	/* push the current block to deflate queue */
	pg->bal++;
	pt_enq_retry(&pg->pt->in, 0, s, PT_DEFAULT_INTERVAL);

	/* fetch copressed block and push to heapqueue to sort */
	pg_block_t *t;
	while(pg->bal > pg->lb) {
		if((t = pt_deq(&pg->pt->out, 0)) == PT_EMPTY) {
			if(pg->bal < pg->ub) { break; }	/* unexpected error */
			sched_yield(); continue;		/* queue full, wait for a while */
		}
		pg->bal--;
		kv_hq_push(pg_hqbin_t, incq_comp, pg->hq, ((pg_hqbin_t){ .id = t->id, .block = t }));
	}

	/* flush heapqueue */
	while(kv_hq_cnt(pg->hq) > 0 && kv_hq_head(pg->hq).id <= pg->cnt.out) {
		pg->cnt.out++;
		t = kv_hq_pop(pg_hqbin_t, incq_comp, pg->hq).block;
		pg_write_block(pg, t);
	}
	return;
}
static _force_inline
size_t pgwrite(pg_t *pg, const void *src, size_t len)
{
	size_t rem = len;
	pg_block_t *s = pg->s;
	if(pg->pt->c->wfp != pg_worker && pt_set_worker(pg->pt, pg, pg_worker)) {
		return(0);
	}

	void (*const write_block[2])(pg_t *, pg_block_t *) = { pg_write_single, pg_write_multi };
	while(rem > 0) {
		/* push the current block to queue and prepare an empty one if the current one is full */
		if(s == NULL || s->head == s->len) {
			if(s != NULL) { write_block[pg->nth > 1](pg, pg_record_write_block(pg, s)); }

			/* create new block */
			size_t malloc_size = _roundup(sizeof(pg_block_t) + pg->block_size, ARCH_HUGEPAGE_SIZE);
			s = aligned_malloc(malloc_size);
			s->head = 0;
			s->len = pg->block_size;
			s->id = pg->cnt.in++;
			s->raw = 1;
			s->flush = 1;
		}

		/* copy the content */
		size_t adv = MIN2(rem, s->len-s->head);
		memcpy(s->buf + s->head, _add_offset(src, len - rem), adv);
		rem -= adv; s->head += adv;
	}
	pg->s = s;
	return len;
}

#if 0
unittest( .name = "pg.single" ) {
	uint64_t const size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size), *q = malloc(size);
	for(uint64_t i = 0; i < size; i++) p[i] = i % 253;

	char const *filename = "./minialign.unittest.pg.tmp";
	FILE *fp = fopen(filename, "w");
	ut_assert(fp != NULL);
	pg_t *pg = pg_init(fp, 1);
	ut_assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgwrite(pg, p, size);
		ut_assert(l == size, "l(%lu), size(%lu)", l, size);
	}
	pg_destroy(pg);
	fclose(fp);


	fp = fopen(filename, "r");
	ut_assert(fp != NULL);
	pg = pg_init(fp, 1);
	ut_assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgread(pg, q, size);
		ut_assert(l == size, "l(%lu), size(%lu)", l, size);
		ut_assert(memcmp(p, q, size) == 0);
	}
	pg_destroy(pg);
	fclose(fp);
	free(p); free(q); remove(filename);
}

unittest( .name = "pg.multi" ) {
	uint64_t const size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size), *q = malloc(size);
	for(uint64_t i = 0; i < size; i++) p[i] = i % 253;

	char const *filename = "./minialign.unittest.pg.tmp";
	FILE *fp = fopen(filename, "w");
	ut_assert(fp != NULL);
	pg_t *pg = pg_init(fp, 4);
	ut_assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgwrite(pg, p, size);
		ut_assert(l == size, "l(%lu), size(%lu)", l, size);
	}
	pg_destroy(pg);
	fclose(fp);


	fp = fopen(filename, "r");
	ut_assert(fp != NULL);
	pg = pg_init(fp, 4);
	ut_assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgread(pg, q, size);
		ut_assert(l == size, "l(%lu), size(%lu)", l, size);
		ut_assert(memcmp(p, q, size) == 0);
	}
	pg_destroy(pg);
	fclose(fp);
	free(p); free(q); remove(filename);
}
#endif

#endif
/**
 * end of pstream.h
 */
