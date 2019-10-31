#ifndef GLOBALS_H
#define GLOBALS_H
#include<stddef.h>

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define HUGE 1.e30
#define TINY 1.e-30


#define ETA 0.02
#define ETA_INI 0.01
#define NMAX 1024

#endif
