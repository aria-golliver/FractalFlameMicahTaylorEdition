#ifndef __HISTOGRAM_H__
#define __HISTOGRAM_H__ __HISTOGRAM_H__

#include "datatypes.h"
#include <immintrin.h>
void histoinit();

f128tuple histohit(f128tuple xyvec, const f128 rvec, const f128 gvec, const f128 bvec, const i32 th_id);

void saveimage();

#define swid 1920
#define shei 1080
#define ss 1

#define hwid (swid * ss)
#define hhei (shei * ss)
#define xshrink 1.0
#define yshrink 1.0

#define xoffset 0.0
#define yoffset 0.5

#endif