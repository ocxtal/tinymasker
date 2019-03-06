#include <stdbool.h>

/* opaque type */
struct re_nfa_s;
typedef struct re_nfa_s re_pattern_t;

typedef struct re_result_s {
        char const *start;
        char const *end;
} re_result_t;

re_pattern_t *  re_compile (char const *);
bool      re_match   (re_pattern_t const *, char const *, re_result_t *);
void      re_free    (re_pattern_t *);

