#ifndef __VMATH_H__
#define __VMATH_H__ __VMATH_H__

#include <immintrin.h>

#define FLOATS_PER_VECTOR_REGISTER 8

static __forceinline __m256 vadd (const __m256 a, const __m256 b) { return _mm256_add_ps(a, b); }
static __forceinline __m256 vsub (const __m256 a, const __m256 b) { return _mm256_sub_ps(a, b); }
static __forceinline __m256 vmul (const __m256 a, const __m256 b) { return _mm256_mul_ps(a, b); }
static __forceinline __m256 vdiv (const __m256 a, const __m256 b) { return _mm256_div_ps(a, b); }

static __forceinline __m128 vadd128 (const __m128 a, const __m128 b) { return _mm_add_ps(a, b); }
static __forceinline __m128 vsub128 (const __m128 a, const __m128 b) { return _mm_sub_ps(a, b); }
static __forceinline __m128 vmul128 (const __m128 a, const __m128 b) { return _mm_mul_ps(a, b); }
static __forceinline __m128 vdiv128 (const __m128 a, const __m128 b) { return _mm_div_ps(a, b); }


static __forceinline __m256 vsqrt    (const __m256 a)                 { return _mm256_sqrt_ps(a);    }
static __forceinline __m256 vexp     (const __m256 a)                 { return _mm256_exp_ps(a);     }
static __forceinline __m256 vinvsqrt (const __m256 a)                 { return _mm256_invsqrt_ps(a); }
static __forceinline __m256 vrsqrt   (const __m256 a)                 { return _mm256_rsqrt_ps(a);   }
static __forceinline __m256 vpow     (const __m256 a, const __m256 b) { return _mm256_pow_ps(a, b);  }


static __forceinline __m256 vsin   (const __m256 a)                 { return _mm256_sin_ps(a);      }
static __forceinline __m256 vcos   (const __m256 a)                 { return _mm256_cos_ps(a);      }
static __forceinline __m256 vsinh  (const __m256 a)                 { return _mm256_sinh_ps(a);     }
static __forceinline __m256 vcosh  (const __m256 a)                 { return _mm256_cosh_ps(a);     }
static __forceinline __m256 vtan   (const __m256 a)                 { return _mm256_tan_ps(a);      }
static __forceinline __m256 vatan2 (const __m256 a, const __m256 b) { return _mm256_atan2_ps(a, b); }

static __forceinline __m256 vtrunc(__m256 a){ return _mm256_round_ps(a, _MM_FROUND_TO_ZERO); }

static __forceinline __m256 vmod(__m256 a, __m256 b){
    // divide a by b
    // round to integer
    // multiply by b
    // return difference between that and a
    __m256 divide     = _mm256_div_ps(a, b);
    __m256 rounded    = _mm256_round_ps(divide, _MM_FROUND_TO_ZERO);
    __m256 multiplied = _mm256_mul_ps(rounded, b);
    __m256 difference = _mm256_sub_ps(a, multiplied);
    return difference;
}


// checks all values for NaN or INFINITY
// adapted from here http://www.gamedev.net/page/resources/_/technical/game-programming/practical-cross-platform-simd-math-part-2-r3101?recache=true
static int vvalid( __m256 v )
{
    __m256 test  = _mm256_mul_ps(v, _mm256_setzero_ps());
    test         = _mm256_cmp_ps(test, _mm256_setzero_ps(), _CMP_EQ_OQ);
    return 0x0ff == _mm256_movemask_ps(test);
}

static __m256 vload(const float *p){
    return _mm256_load_ps(p);
}

static void vstore(float *p, __m256 v){
    _mm256_store_ps(p, v);
}

static __m128 vload128(const float *p){
    return _mm_load_ps(p);
}

static void vstore128(float *p, __m128 v){
    _mm_store_ps(p, v);
}
#endif