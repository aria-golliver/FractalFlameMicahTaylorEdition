#ifndef __RDRAND_H__
#define __RDRAND_H__ __RDRAND_H__

#include "datatypes.h"
#include "systemconfig.h"

#if __RDRAND_AVAILIABLE__ != 0
	inline u16 rdrand_u16(){
		__asm rdrand ax;
	}

	inline u32 rdrand_u32(){
		__asm rdrand eax;
	}

	inline u64 rdrand_u64(){
		__asm rdrand rax;
	}

	inline i16 rdrand_i16(){
		__asm rdrand ax;
	}

	inline i32 rdrand_i32(){
		__asm rdrand eax;
	}

	inline i64 rdrand_i64(){
		__asm rdrand rax;
	}

	inline f32 rdrand_f32(){
		return (f32)(((f64)rdrand_i32()) / INT32_MAX);
	}

	inline f64 rdrand_f64(){
		return ((f64)rdrand_i64()) / INT64_MAX;
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