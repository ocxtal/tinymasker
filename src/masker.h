
/**
 * @file masker.h
 *
 * @author Hajime Suzuki
 * @license MIT
 */

#ifndef _MASKER_H_INCLUDED
#define _MASKER_H_INCLUDED

#ifndef _ALIGN_H_INCLUDED
#  error "#include \"align.h\" must be before #include \"masker.h\""
#endif

#if 0
/* printer */
typedef struct {
	uint32_t flip;				/* flip reference and query */
} tm_print_conf_t;

typedef struct {
	uint32_t flip;
} tm_print_t;
#endif

// void tm_print_destory_static(tm_print_t *self);
// void tm_print_init_static(tm_print_t *self, tm_print_conf_t const *conf, char const *args);



/* multithreading */
typedef struct tm_mtscan_s tm_mtscan_t;

tm_mtscan_t *tm_mtscan_init(tm_idx_t const *mi/*, tm_print_t *printer*/, pt_t *pt);
void tm_mtscan_destroy(tm_mtscan_t *self);
int tm_mtscan_file(tm_mtscan_t *self, char const *fn);


#endif		/* _MASKER_H_INCLUDED */

/**
 * end of masker.h
 */
