#ifndef __RDRAND_H__
#define __RDRAND_H__ __RDRAND_H__

#include "datatypes.h"
#include <stdio.h>
#include <stdlib.h>

#define RDRAND_SUCCESS 1

__forceinline void rdrand_u16(u16 *p){
    while( _rdrand16_step(p) != RDRAND_SUCCESS);
}

__forceinline void rdrand_u32(u32 *p){
    while( _rdrand32_step(p) != RDRAND_SUCCESS);
}

__forceinline void rdrand_u64(u64 *p){
    while( _rdrand64_step(p) != RDRAND_SUCCESS);
}

__forceinline void rdrand_i16(i16 *p){
    rdrand_u16((u16 *)p);
}

__forceinline void rdrand_i32(i32 *p){
    rdrand_u32((u32 *)p);
}

__forceinline void rdrand_i64(i64 *p){
    rdrand_u64((u64 *)p);
}

// floating point random functions return (-1, 1)
__forceinline void rdrand_f32(f32 *p){
    i32 tmp;
    rdrand_i32(&tmp);
    *p = tmp / (f32)INT32_MAX;
}

// floating point random functions return (-1, 1)
__forceinline void rdrand_f64(f64 *p){
    i64 tmp;
    rdrand_i64(&tmp);
    *p = tmp / (f64)INT64_MAX;
}

#endif