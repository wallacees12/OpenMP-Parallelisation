//
// Created by Juan Rodriguez on 2025-03-02.
//

#ifndef HPP_TIMER_H
#define HPP_TIMER_H
#include <sys/time.h>

/* The argument now should be a double (not a pointer to a double) */
#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}

#endif //HPP_TIMER_H
