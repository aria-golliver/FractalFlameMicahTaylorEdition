#ifndef __DATATYPES_H__
#define __DATATYPES_H__ __DATATYPES_H__

#include <stdint.h>
#include <immintrin.h>

typedef uint8_t u8  ;
typedef uint16_t u16 ;
typedef uint32_t u32 ;
typedef uint64_t u64 ;

typedef int8_t i8  ;
typedef int16_t i16 ;
typedef int32_t i32 ;
typedef int64_t i64 ;

typedef float f32 ;
typedef double f64 ;

typedef union {
    u32 data;
    struct {
        u8 r;
        u8 g;
        u8 b;
        u8 a;
    };
} rgba8;

typedef struct GLRGB8_t{
    u8 c[3];
} GLRGB8;

union f256_t {
    _declspec(align(16)) f32 f[8];
    __m256 v;
};

typedef union f256_t f256;

typedef union {
    _declspec(align(16)) struct {
        f32 r, g, b, a;
    };
    __m128 vec;
} histocell;

typedef union {
    struct {
        f32 r, g, b, a;
    };
    _declspec(align(16)) __m128 vec;
} colorset;

typedef struct {
    f32 x, y;
} point;

typedef struct {
    _declspec(align(32)) f256 x;
    _declspec(align(32)) f256 y;
} f256tuple;

typedef struct {
    f32 a, b, c, d, e, f;
    _declspec(align(16)) colorset color;
} affinematrix;

#endif