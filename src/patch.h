
/**
 * @file patch.h
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#ifndef _PATCH_H_INCLUDED
#define _PATCH_H_INCLUDED


/* configurations: make user-provided value different from default one */
#define tm_patch_wrap(x)				( (x) + 1 )
#define tm_patch_unwrap(x)				( (x) - 1 )
#define tm_patch_is_default(x)			( (x) == 0 )

typedef struct {
	size_t max_gap_len;
	size_t margin_len;
} tm_patch_conf_t;



#endif	/* _PATCH_H_INCLUDED */
/**
 * end of patch.h
 */
