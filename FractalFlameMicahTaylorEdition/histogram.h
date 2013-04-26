#ifndef __HISTOGRAM_H__
#define __HISTOGRAM_H__ __HISTOGRAM_H__

#include "datatypes.h"
#include <immintrin.h>
void histoinit();

f128tuple histohit(f128tuple xyvec, const f128 rvec, const f128 gvec, const f128 bvec, const i32 th_id);

void saveimage();

#endif