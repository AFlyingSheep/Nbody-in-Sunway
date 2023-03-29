#ifndef _SPH_UTILS_H_
#define _SPH_UTILS_H_

#include <stdio.h>
#include <string.h>
#include <time.h> // for clock_gettime

#ifdef __sw_64__

#include <crts.h>
#define get_cycles CRTS_time_cycle

#elif defined(x86)

#include <x86intrin.h>
static unsigned long get_cycles(void) { return __rdtsc(); }

#else

#endif /* __sw_64__ */

static double dtime(void) {
  struct timespec tp;
  double ret;
  int err;
  err = clock_gettime(CLOCK_REALTIME, &tp);
  if (err == 0) {
    ret = ((double)tp.tv_sec + tp.tv_nsec * 1.e-9);
  } else {
    ret = 0.0;
  }
  return ret;
}

#endif /* ifndef _SPH_UTILS_H_ */
