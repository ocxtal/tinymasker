#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>

#include "re.h"

#define L(c) ((c) >> 8)
#define H(c) ((c) &  0xFF)

typedef enum {
        NFA_CHAR,
        NFA_EPSILON,
        NFA_ANYCHAR,
        NFA_CLASS,
        NFA_NCLASS,
        NFA_BEGIN,
        NFA_END,
} re_nfa_t;
typedef enum {
        RE_INVALID,
        RE_CHAR,
        RE_ALT,
        RE_STAR,
        RE_PLUS,
        RE_OPTION,
        RE_DOT,
        RE_CONCAT,
        RE_CLASS,
        RE_NCLASS,
        RE_BEGIN,
        RE_END
} re_type_t;

struct re_s {
        re_type_t type;
        union {
                uint8_t c;

                struct {
                        uint8_t n;
                        uint16_t *class;
                };

                struct {
                        struct re_s *left;
                        struct re_s *right;
                };

                struct re_s *re;
        };
};

struct transition_s {
        uint8_t t;
        uint8_t c;
        uint16_t *class;
        union { struct st_s *s; size_t idx; };
};

struct st_s {
        struct transition_s one, two;
};


struct re_nfa_s {
        struct st_s *states;
        size_t count;
        size_t alloc;
};

struct frame_s {
        char const *s;
        struct st_s const *state;
};

/*
 * Attempt to allocate a new `struct re`.
 * Return a null-pointer if the attempt fails.
 */
#define mkre(name) \
        name = malloc(sizeof *name); \
        if (name == NULL) return NULL;


size_t
addstate(struct re_nfa_s *nfa)
{
        if (nfa->count == nfa->alloc) {
                nfa->alloc = nfa->alloc ? nfa->alloc * 2 : 4;
                struct st_s *tmp = realloc(nfa->states, nfa->alloc * sizeof *tmp);
                if (tmp == NULL) {
                        assert(false);
                }
                nfa->states = tmp;
        }

        nfa->states[nfa->count].one.idx = -1;
        nfa->states[nfa->count].two.idx = -1;

        nfa->count += 1;

        return nfa->count - 1;
}

static struct transition_s *
transition(struct re_nfa_s *nfa, size_t from, size_t to, uint16_t type)
{
        if (nfa->states[from].one.idx == SIZE_MAX) {
                nfa->states[from].one.idx = to;
                nfa->states[from].one.t = type;
                return &nfa->states[from].one;
        } else if (nfa->states[from].two.idx == SIZE_MAX) {
                nfa->states[from].two.idx = to;
                nfa->states[from].two.t = type;
                return &nfa->states[from].two;
        } else {
                assert(false);
        }
        return(NULL);
}

static void
complete(struct re_nfa_s *nfa)
{
        for (size_t i = 0; i < nfa->count; ++i) {
                if (nfa->states[i].one.idx != SIZE_MAX) {
                        nfa->states[i].one.s = &nfa->states[nfa->states[i].one.idx];
                } else {
                        nfa->states[i].one.s = NULL;
                }
                if (nfa->states[i].two.idx != SIZE_MAX) {
                        nfa->states[i].two.s = &nfa->states[nfa->states[i].two.idx];
                } else {
                        nfa->states[i].two.s = NULL;
                }
        }
}

static size_t
tonfa(struct re_nfa_s *nfa, size_t start, struct re_s *re)
{
        struct transition_s *tr;
        size_t a, b, c;
        size_t t, v;

        switch (re->type) {
        case RE_CHAR:
                a = addstate(nfa);
                transition(nfa, start, a, RE_CHAR)->c = re->c;
                return a;
        case RE_BEGIN:
                a = addstate(nfa);
                transition(nfa, start, a, NFA_BEGIN);
                return a;
        case RE_END:
                a = addstate(nfa);
                transition(nfa, start, a, NFA_END);
                return a;
        case RE_CLASS:
                a = addstate(nfa);
                tr = transition(nfa, start, a, NFA_CLASS);
                tr->class = re->class;
                tr->c = re->n;
                return a;
        case RE_NCLASS:
                a = addstate(nfa);
                tr = transition(nfa, start, a, NFA_NCLASS);
                tr->class = re->class;
                tr->c = re->n;
                return a;
        case RE_ALT:
                /* End state */
                c = addstate(nfa);
                
                /* Start of left */
                a = addstate(nfa);
                /* Start of right */
                b = addstate(nfa);

                /* End of left */
                t = tonfa(nfa, a, re->left);
                /* End of right */
                v = tonfa(nfa, b, re->right);

                /* Link left to end */
                transition(nfa, t, c, NFA_EPSILON);
                /* Link right to end */
                transition(nfa, v, c, NFA_EPSILON);

                /* Link start to left */
                transition(nfa, start, a, NFA_EPSILON);
                /* Link start to right */
                transition(nfa, start, b, NFA_EPSILON);

                return c;
        case RE_STAR:
                /* End state */
                b = addstate(nfa);

                /* Star's operand */
                a = addstate(nfa);

                /* End of star's operand */
                t = tonfa(nfa, a, re->re);

                /* Make the loop (connect end to start) */
                transition(nfa, t, a, NFA_EPSILON);

                /* The other way out: the end state */
                transition(nfa, t, b, NFA_EPSILON);

                /* Link start to star's operand (match 1 or more) */
                transition(nfa, start, a, NFA_EPSILON);

                /* Link start directly to end (match 0 times) */
                transition(nfa, start, b, NFA_EPSILON);

                return b;
        case RE_PLUS:
                /* End state */
                b = addstate(nfa);

                /* Plus's operand */
                a = addstate(nfa);

                /* End of plus's operand */
                t = tonfa(nfa, a, re->re);

                /* Make the loop (connect end to start) */
                transition(nfa, t, a, NFA_EPSILON);

                /* The other way out: the end state */
                transition(nfa, t, b, NFA_EPSILON);

                /* Link start to plus's operand (match 1 or more) */
                transition(nfa, start, a, NFA_EPSILON);

                return b;
        case RE_OPTION:
                /* End state */
                b = addstate(nfa);

                /* ?'s operand */
                a = addstate(nfa);

                /* End of ?'s operand */
                t = tonfa(nfa, a, re->re);

                /* Link end of ?'s operand to the end state */
                transition(nfa, t, b, NFA_EPSILON);

                /* Link start directly to end (no match) */
                transition(nfa, start, b, NFA_EPSILON);

                /* Link start to ?'s operand (match) */
                transition(nfa, start, a, NFA_EPSILON);

                return b;
        case RE_DOT:
                /* End state */
                b = addstate(nfa);

                /* Link start to end */
                transition(nfa, start, b, NFA_ANYCHAR);

                return b;
        case RE_CONCAT:
                t = tonfa(nfa, start, re->left);
                return tonfa(nfa, t, re->right);
        }

        /* never reach here */
        return 0;
}


static struct re_s *regexp(char const **, bool allow_trailing);
static struct re_s *subexp(char const **);
static struct re_s *atom(char const **);
static struct re_s *charclass(char const **);

static void
freere(struct re_s *re)
{
        if (re->type == RE_ALT || re->type == RE_CONCAT) {
                freere(re->left);
                freere(re->right);
        }

        free(re);
}

static struct re_s *
or(struct re_s *left, struct re_s *right)
{
        if (left == NULL) {
                return right;
        }
        if (right == NULL) {
                return left;
        }

        struct re_s *e;
        mkre(e);
        e->type  = RE_ALT;
        e->left  = left;
        e->right = right;
        return e;
}

static struct re_s *
and(struct re_s *left, struct re_s *right)
{
        if (left == NULL) {
                return right;
        }
        if (right == NULL) {
                return left;
        }

        struct re_s *e;
        mkre(e);
        e->type  = RE_CONCAT;
        e->left  = left;
        e->right = right;
        return e;
}

static struct re_s *
regexp(char const **s, bool allow_trailing)
{
        struct re_s *re, *e;

        e = NULL;
        while (**s && **s != '|') {
                struct re_s *sub = subexp(s);
                if (sub == NULL) {
                        if (allow_trailing) {
                                return e;
                        } else {
                                return NULL;
                        }
                }

                e = and(e, sub);
                if (e == NULL) {
                        return NULL;
                }
        }

        if (!**s) {
                return e;
        }

        *s += 1;

        mkre(re);
        re->type  = RE_ALT;
        re->left  = e;
        re->right = regexp(s, allow_trailing);

        if (re->right == NULL) {
                return NULL;
        }

        return re;
}

static struct re_s *
subexp(char const **s)
{
        struct re_s *re, *e = atom(s);
        if (e == NULL) {
                return NULL;
        }
        
        re_type_t type = RE_INVALID;

        switch (**s) {
        case '*': ++*s; type = RE_STAR;   break;
        case '+': ++*s; type = RE_PLUS;   break;
        case '?': ++*s; type = RE_OPTION; break;
        default:                          break;
        }

        if (type == RE_INVALID) {
                return e;
        }

        mkre(re);
        re->type = type;
        re->re   = e;

        return re;
}

static struct re_s *
atom(char const **s)
{
        if (**s == '\0' || **s == ')') {
                return NULL;
        }

        /*
         * If the next character is a left
         * parenthesis, then we must match
         * a `regexp`; otherwise, we match
         * one character.
         */
        struct re_s *e;
        if (**s == '(') {
                *s += 1;
                e = regexp(s, true);
                if (**s != ')') {
                        return NULL;
                }
                *s += 1;
                return e;
        } else if (**s == '.') {
                mkre(e);
                e->type = RE_DOT;
                *s += 1;
                return e;
        } else if (**s == '^') {
                mkre(e);
                e->type = RE_BEGIN;
                *s += 1;
                return e;
        } else if (**s == '$') {
                mkre(e);
                e->type = RE_END;
                *s += 1;
                return e;
        } else if (**s == '[') {
                *s += 1;
                return charclass(s);
        } else if (**s == '\\') {
                *s += 1;
                if (**s == '\0') {
                        /* The regular expression cannot end with a backslash */
                        return NULL;
                }
                mkre(e);
                e->type = RE_CHAR;
                e->c    = **s;
                *s += 1;
                return e;
        } else {
                mkre(e);
                e->type = RE_CHAR;
                e->c    = **s;
                *s += 1;
                return e;
        }
}

static bool
searchclass(uint16_t const *class, size_t n, uint8_t c)
{
        int64_t lo = 0;
        int64_t hi = n - 1;

        while (lo <= hi) {
                int64_t m = (lo + hi) / 2;
                if      (c < L(class[m])) hi = m - 1;
                else if (c > H(class[m])) lo = m + 1;
                else                      return true;
        }

        return false;
}

static int
classcmp(void const *ap, void const *bp)
{
        uint16_t const *a = ap;
        uint16_t const *b = bp;
        return (*a << 8) - (*b << 8);
}

static struct re_s *
charclass(char const **s)
{
        if (**s == '\0')
                return NULL;

        bool negate = **s == '^';
        if (negate)
                *s += 1;

        uint8_t n = 0;
        for (size_t i = 0; (*s)[i] != '\0'; ++n) {
                if ((*s)[i] == ']' && i != 0)
                        break;
                if ((*s)[i + 1] == '-' && (*s)[i + 2] != '\0' && (*s)[i + 2] != ']')
                        i += 3;
                else
                        i += 1;
        }

        uint16_t *class = malloc(n * sizeof *class);
        if (class == NULL)
                return NULL;

        for (size_t i = 0; i < n; ++i, ++*s) {
                if ((*s)[1] == '-' && (*s)[2] != '\0' && (*s)[2] != ']') {
                        class[i] = ((*s)[0] << 8) + (*s)[2];
                        *s += 2;
                } else {
                        class[i] = (**s << 8) + **s;
                }
        }

        if (**s == '\0') {
                free(class);
                return NULL;
        }

        *s += 1;


        qsort(class, n, sizeof *class, classcmp);
        
        for (size_t i = 0; i < n; ++i) {
                uint8_t high = H(class[i]);
                size_t j = i + 1;
                while (j < n && L(class[j]) <= high + 1) {
                        if (H(class[j]) > high)
                                high = H(class[j]);
                        ++j;
                }
                class[i] = (L(class[i]) << 8) + high;
                memmove(class + i + 1, class + j, (n - j) * sizeof *class);
                n -= (j - i - 1);
        }

        struct re_s *e;
        mkre(e);
        e->type = negate? RE_NCLASS: RE_CLASS;
        e->class = class;
        e->n = n;

        return e;
}

static struct re_s *
parse(char const *s)
{
        return regexp(&s, false);
}

inline static bool
charmatch(char const **s, struct transition_s const *tr, char const *begin)
{
        switch (tr->t) {
        case NFA_EPSILON: return true;
        case NFA_ANYCHAR: return (**s != '\0') && ++*s;
        case NFA_BEGIN:   return *s == begin;
        case NFA_END:     return **s == '\0';
        case NFA_CLASS:   return **s != '\0' && searchclass(tr->class, tr->c, **s) && ++*s;
        case NFA_NCLASS:  return **s != '\0' && !searchclass(tr->class, tr->c, **s) && ++*s;
        case NFA_CHAR:    return **s == tr->c && ++*s;
        default:          assert(false);
        }
        return false;
}

static char const *
domatch(struct st_s const *state, char const *string, char const *begin)
{
        struct frame_s *stack = malloc(8 * sizeof *stack);
        size_t capacity = 8;
        size_t i = 0;

        if (stack == NULL)
                goto end;

        stack[i++] = (struct frame_s){ .s = string, .state = state };

        while (i != 0) {
                struct frame_s f =  stack[--i];
                char const *s = f.s;

                if (i + 2 >= capacity) {
                        capacity *= 2;
                        struct frame_s *tmp = realloc(stack, capacity * sizeof *stack);
                        if (tmp == NULL)
                                goto end;
                        else
                                stack = tmp;
                }

                if (f.state->one.s == NULL)
                        return f.s;
                if (f.state->two.s != NULL && charmatch(&s, &f.state->two, begin))
                        stack[i++] = (struct frame_s){ .s = s, .state = f.state->two.s };
                if (charmatch(&f.s, &f.state->one, begin))
                        stack[i++] = (struct frame_s){ .s = f.s, .state = f.state->one.s };
        }

end:
        free(stack);
        return NULL;
}

bool
re_match(struct re_nfa_s const *nfa, char const *s, struct re_result_s *result)
{
        char const *end;
        if (end = domatch(nfa->states, s, s), end != NULL) {
                if (result != NULL) {
                        result->start = s;
                        result->end   = end;
                }
                return true;
        }

        while (*++s) {
                if (end = domatch(nfa->states, s, NULL), end != NULL) {
                        if (result != NULL) {
                                result->start = s;
                                result->end   = end;
                        }
                        return true;
                }
        }

        return false;
}

re_pattern_t *
re_compile(char const *s)
{
        assert(s);

        struct re_nfa_s *nfa = malloc(sizeof *nfa);
        if (nfa == NULL) {
                return NULL;
        }

        nfa->states = NULL;
        nfa->count = 0;
        nfa->alloc = 0;

        struct re_s *re = parse(s);
        if (re == NULL) {
                free(nfa);
                return NULL;
        }

        addstate(nfa);
        tonfa(nfa, 0, re);
        complete(nfa);

        freere(re);

        return nfa;
}

void
re_free(struct re_nfa_s *nfa)
{
        free(nfa->states);
        free(nfa);
}
