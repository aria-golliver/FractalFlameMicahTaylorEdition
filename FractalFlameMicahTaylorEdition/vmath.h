#ifndef __VMATH_H__
#define __VMATH_H__ __VMATH_H__

#include <immintrin.h>
/*
static __forceinline __m128 vadd(const f128 a, const f128 b) { return _mm_add_ps(a.v, b.v); }
static __forceinline __m128 vsub(const f128 a, const f128 b) { return _mm_sub_ps(a.v, b.v); }
static __forceinline __m128 vmul(const f128 a, const f128 b) { return _mm_mul_ps(a.v, b.v); }
static __forceinline __m128 vdiv(const f128 a, const f128 b) { return _mm_div_ps(a.v, b.v); }


static __forceinline __m128 vsqrt(const f128 a) { return _mm_sqrt_ps(a.v); }
static __forceinline __m128 vinvsqrt(const f128 a) { return _mm_invsqrt_ps(a.v); }
static __forceinline __m128 vrsqrt(const f128 a) { return _mm_rsqrt_ps(a.v); }


static __forceinline __m128 vsin(const f128 a) { return _mm_sin_ps(a.v); }
static __forceinline __m128 vcos(const f128 a) { return _mm_cos_ps(a.v); }
static __forceinline __m128 vatan2(const f128 a, const f128 b) { return _mm_atan2_ps(a.v, b.v); }

*/
static __forceinline __m128 vadd(const __m128 a, const __m128 b) { return _mm_add_ps(a, b); }
static __forceinline __m128 vsub(const __m128 a, const __m128 b) { return _mm_sub_ps(a, b); }
static __forceinline __m128 vmul(const __m128 a, const __m128 b) { return _mm_mul_ps(a, b); }
static __forceinline __m128 vdiv(const __m128 a, const __m128 b) { return _mm_div_ps(a, b); }


static __forceinline __m128 vsqrt(const __m128 a) { return _mm_sqrt_ps(a); }
static __forceinline __m128 vexp(const __m128 a) { return _mm_exp_ps(a); }
static __forceinline __m128 vinvsqrt(const __m128 a) { return _mm_invsqrt_ps(a); }
static __forceinline __m128 vrsqrt(const __m128 a) { return _mm_rsqrt_ps(a); }
static __forceinline __m128 vpow(const __m128 a, const __m128 b) { return _mm_pow_ps(a, b); }


static __forceinline __m128 vsin(const __m128 a) { return _mm_sin_ps(a); }
static __forceinline __m128 vcos(const __m128 a) { return _mm_cos_ps(a); }
static __forceinline __m128 vsinh(const __m128 a) { return _mm_sinh_ps(a); }
static __forceinline __m128 vcosh(const __m128 a) { return _mm_cosh_ps(a); }
static __forceinline __m128 vtan(const __m128 a) { return _mm_tan_ps(a); }
static __forceinline __m128 vatan2(const __m128 a, const __m128 b) { return _mm_atan2_ps(a, b); }

static __forceinline __m128 vmod(__m128 a, __m128 b){
    // divide a by b
    // round to integer
    // multiply by b
    // return difference between that and a
    __m128 divide = _mm_div_ps(a, b);
    __m128 rounded = _mm_round_ps(divide, _MM_FROUND_TO_ZERO);
    __m128 multiplied = _mm_mul_ps(rounded, b);
    __m128 difference = _mm_sub_ps(a, multiplied);
    return difference;
}

#endif