#ifndef __VMATH_H__
#define __VMATH_H__ __VMATH_H__

#include <immintrin.h>

static __forceinline __m128 vadd (const __m128 a, const __m128 b) { return _mm_add_ps(a, b); }
static __forceinline __m128 vsub (const __m128 a, const __m128 b) { return _mm_sub_ps(a, b); }
static __forceinline __m128 vmul (const __m128 a, const __m128 b) { return _mm_mul_ps(a, b); }
static __forceinline __m128 vdiv (const __m128 a, const __m128 b) { return _mm_div_ps(a, b); }


static __forceinline __m128 vsqrt    (const __m128 a)                 { return _mm_sqrt_ps(a);    }
static __forceinline __m128 vexp     (const __m128 a)                 { return _mm_exp_ps(a);     }
static __forceinline __m128 vinvsqrt (const __m128 a)                 { return _mm_invsqrt_ps(a); }
static __forceinline __m128 vrsqrt   (const __m128 a)                 { return _mm_rsqrt_ps(a);   }
static __forceinline __m128 vpow     (const __m128 a, const __m128 b) { return _mm_pow_ps(a, b);  }


static __forceinline __m128 vsin   (const __m128 a)                 { return _mm_sin_ps(a);      }
static __forceinline __m128 vcos   (const __m128 a)                 { return _mm_cos_ps(a);      }
static __forceinline __m128 vsinh  (const __m128 a)                 { return _mm_sinh_ps(a);     }
static __forceinline __m128 vcosh  (const __m128 a)                 { return _mm_cosh_ps(a);     }
static __forceinline __m128 vtan   (const __m128 a)                 { return _mm_tan_ps(a);      }
static __forceinline __m128 vatan2 (const __m128 a, const __m128 b) { return _mm_atan2_ps(a, b); }

static __forceinline __m128 vtrunc(__m128 a){ return _mm_round_ps(a, _MM_FROUND_TO_ZERO); }

static __forceinline __m128 vmod(__m128 a, __m128 b){
    // divide a by b
    // round to integer
    // multiply by b
    // return difference between that and a
    __m128 divide     = _mm_div_ps(a, b);
    __m128 rounded    = _mm_round_ps(divide, _MM_FROUND_TO_ZERO);
    __m128 multiplied = _mm_mul_ps(rounded, b);
    __m128 difference = _mm_sub_ps(a, multiplied);
    return difference;
}


// checks all values for NaN or INFINITY
// adapted from here http://www.gamedev.net/page/resources/_/technical/game-programming/practical-cross-platform-simd-math-part-2-r3101?recache=true
static int vvalid( __m128 v )
{
    __m128 test  = _mm_mul_ps(v, _mm_setzero_ps());
    test         = _mm_cmpeq_ps(test, _mm_setzero_ps());
    return 0x0f == _mm_movemask_ps(test);
}
#endif