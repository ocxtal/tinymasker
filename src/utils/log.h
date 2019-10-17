
/**
 * @file log.h
 */

#ifndef _LOG_H_INCLUDED
#define _LOG_H_INCLUDED

/* include global header *before* we include individual dependencies */
#include "common.h"
#include <sys/time.h>
#include <sys/resource.h>		/* getrusage */

#if defined(__linux__)
#  define LOGGER_MAXRSS_COEF	( 1000.0 )
#elif defined(__APPLE__)
#  define LOGGER_MAXRSS_COEF	( 1000.0 * 1000.0 )
#else
#  define LOGGER_MAXRSS_COEF	( 1000.0 )
#endif


/* time functions */
static _force_inline
double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return(
		  r.ru_utime.tv_sec
		+ r.ru_stime.tv_sec
		+ 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec)
	);
}

static _force_inline
double realtime(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return(tp.tv_sec + tp.tv_usec * 1e-6);
}

#if defined(LOGGER_MAIN)
double inittime = 0.0;
#else
extern double inittime;
#endif


#define logger_init()			{ inittime = realtime(); }
#define logger_destroy() { \
	struct rusage r; \
	getrusage(RUSAGE_SELF, &r); \
	message(stderr, "Real time: %.3f sec; CPU: %.3f sec; maxrss: %.1f MB", realtime() - inittime, cputime(), (double)r.ru_maxrss / LOGGER_MAXRSS_COEF); \
}

#define error(...)				{ logger_impl(stderr, 'E', __VA_ARGS__, ""); }
#define warn(...)				{ logger_impl(stderr, 'W', __VA_ARGS__, ""); }
#define message(_fp, ...)		{ if((_fp) != NULL) { logger_time_impl(_fp, 'M', __VA_ARGS__, ""); } }

#define logger_time_impl(fp, level, fmt, ...) { \
	logger_printf(fp, "[%c::%s::%.3f*%.2f] " fmt "%s\n", level, __func__, realtime() - inittime, cputime() / (realtime() - inittime), __VA_ARGS__); \
}
#define logger_impl(fp, level, fmt, ...) { \
	logger_printf(fp, "[%c::%s] " fmt "%s\n", level, __func__, __VA_ARGS__); \
}
#define logger_printf(fp, fmt, ...) { \
	fprintf(fp, fmt, __VA_ARGS__); \
}

#endif
/**
 * end of log.h
 */
