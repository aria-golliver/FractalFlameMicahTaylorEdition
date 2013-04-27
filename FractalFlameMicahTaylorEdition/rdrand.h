#ifndef __RDRAND_H__
#define __RDRAND_H__ __RDRAND_H__

#include "datatypes.h"
#include "systemconfig.h"
#include <stdio.h>
#if __RDRAND_AVAILIABLE__ != 0
#define RDRAND_SUCCESS 1

inline u16 rdrand_u16(u16 *p){
    while( _rdrand16_step(p) != RDRAND_SUCCESS);
    return *p;
}

inline u32 rdrand_u32(u32 *p){
    while( _rdrand32_step(p) != RDRAND_SUCCESS);
    return *p;
}

inline u64 rdrand_u64(u64 *p){
    while( _rdrand64_step(p) != RDRAND_SUCCESS);
    return *p;
}

inline i16 rdrand_i16(i16 *p){
    return rdrand_u16((u16 *)p);
}

inline i32 rdrand_i32(i32 *p){
    return rdrand_u32((u32 *)p);
}

inline i64 rdrand_i64(i64 *p){
    return rdrand_u64((u64 *)p);
}

inline f32 rdrand_f32(f32 *p){
    i32 i = rdrand_i32((i32 *)p);
    return (f32)((f64)i / (f64)INT32_MAX);
}

inline f64 rdrand_f64(f64 *p){
    i64 i = rdrand_i64((i64 *) p);
    return ((f64)i / (f64)INT64_MAX);
}
#else
#include <stdlib.h>

inline u16 rdrand_u16(){
    return rand();
}

inline u32 rdrand_u32(){
    return rand();
}

inline u64 rdrand_u64(){
    return rand();
}

inline i16 rdrand_i16(){
    return rand();
}

inline i32 rdrand_i32(){
    return rand();
}

inline i64 rdrand_i64(){
    return rand();
}

inline f32 rdrand_f32(){
    return (float)((double)rand()/RAND_MAX);
}

inline f64 rdrand_f64(){
    return (double)rand()/RAND_MAX;
}
#endif

#endif