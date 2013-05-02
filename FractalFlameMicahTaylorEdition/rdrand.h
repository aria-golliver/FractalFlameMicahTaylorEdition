#ifndef __RDRAND_H__
#define __RDRAND_H__ __RDRAND_H__

#include "datatypes.h"
#include <stdio.h>
#include <stdlib.h>

#define RDRAND_SUCCESS 1

__forceinline u16 rdrand_u16(u16 *p){
    while( _rdrand16_step(p) != RDRAND_SUCCESS);
    return *p;
}

__forceinline u32 rdrand_u32(u32 *p){
    while( _rdrand32_step(p) != RDRAND_SUCCESS);
    return *p;
}

__forceinline u64 rdrand_u64(u64 *p){
    while( _rdrand64_step(p) != RDRAND_SUCCESS);
    return *p;
}

__forceinline i16 rdrand_i16(i16 *p){
    return rdrand_u16((u16 *)p);
}

__forceinline i32 rdrand_i32(i32 *p){
    return rdrand_u32((u32 *)p);
}

__forceinline i64 rdrand_i64(i64 *p){
    return rdrand_u64((u64 *)p);
}

// floating point random functions return (-1, 1)
__forceinline f32 rdrand_f32(f32 *p){
    i32 i = rdrand_i32((i32 *)p);
    return (f32)((f64)i / (f64)INT32_MAX);
}

// floating point random functions return (-1, 1)
__forceinline f64 rdrand_f64(f64 *p){
    i64 i = rdrand_i64((i64 *) p);
    return ((f64)i / (f64)INT64_MAX);
}

#endif