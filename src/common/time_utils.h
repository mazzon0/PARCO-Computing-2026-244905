#ifndef TIME_UTILS
#define TIME_UTILS

#include <time.h>
#include <stdlib.h>

double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

#endif
