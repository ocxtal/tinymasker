
/**
 * @file dbg.h
 * @brief de Bruijn hash structure for reference sequences
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#ifndef _DBG_H_INCLUDED
#define _DBG_H_INCLUDED


/*
 * index data structure; we expect each repeat sequence is shorter than 32kbp
 */
typedef struct {
	uint32_t kbits;
	uint32_t max_amb_col;			/* maximum # of ambiguous columns */
	struct {
		size_t head, tail;
	} margin;
} tm_ref_conf_t;

typedef struct {
	uint32_t size, head_margin;		/* object size */
	uint32_t kbits, slen;			/* kmer size and sequence length */
} tm_ref_sketch_t;


/* forward decl for APIs */
typedef struct tm_ref_tbuf_s tm_ref_tbuf_t;
tm_ref_tbuf_t *tm_ref_tbuf_init(void);
void tm_ref_tbuf_destroy(tm_ref_tbuf_t *tbuf);

tm_ref_sketch_t *tm_ref_sketch(tm_ref_tbuf_t *self, tm_ref_conf_t const *conf, uint8_t const *seq, size_t slen);
void tm_ref_destroy(tm_ref_sketch_t *sk);



/*
 * 2nd stage hash array. put here so that matcher functions below are inlined.
 */
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

#define TM_REF_ALIGN_SIZE		( 2 * sizeof(uint64_t) )
#define TM_REF_BASE_OFS			( sizeof(tm_ref_sketch_t) )



/*
 * matcher implemetation comes here for better performance
 * we expect the following lines are inlined
 */
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

static _force_inline
tm_ref_next_t tm_ref_match_next(tm_ref_state_t s, uint64_t keep, uint8_t next)
{
	int64_t const pitch = TM_REF_ALIGN_SIZE;

	/* get link to the next bin */
	tm_ref_link_t const *lk = _add_offset(s.bin->link, next);
	tm_ref_bin_t const *bin = _add_offset(s.bin, pitch * lk->next);	/* signed */

	/* update matching state */
	uint64_t const prefix     = MAX2(s.prefix + s.unmatching, lk->plen);
	uint64_t const unmatching = 0ULL - (prefix != 0);

	/* calc squashable subbin of the previous node */
	tm_ref_squash_t const sq = tm_ref_load_squash(s.bin,
		unmatching - keep,		/* do not squash the previous if the current bin is unmatching; also when keep == 1 */
		next
	);

	/* done */
	tm_ref_next_t const r = {
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

#endif	/* _DBG_H_INCLUDED */
/**
 * end of dbg.h
 */
