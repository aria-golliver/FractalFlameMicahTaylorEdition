#ifndef __VARIATIONS_H__
#define __VARIATIONS_H__ __VARIATIONS_H__

#include "vmath.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#endif


#define v1											\
	sumvecx = vadd(sumvecx, vsin(affinedx));		\
	sumvecy = vadd(sumvecy, vcos(affinedy));	

#define v2											\
	sumvecx = vadd(sumvecx, vdiv(affinedx, rsq));	\
	sumvecy = vadd(sumvecy, vdiv(affinedy, rsq));

#define v3	{																		\
	sumvecx = vadd(sumvecx, vsub(vmul(affinedx, sinrsq), vmul(affinedy, cosrsq)));	\
	sumvecy = vadd(sumvecy, vsub(vmul(affinedx, cosrsq), vmul(affinedy, sinrsq)));	\
}

#define v4																\
	__m128 newx = vdiv(													\
					vmul(												\
						vsub(affinedx, affinedy),						\
						vadd(affinedx, affinedy)),						\
					r);													\
																		\
	__m128 newy = vdiv(													\
					vmul(												\
						vmul(affinedx, affinedy),						\
						twovec),										\
					r);													\
																		\
	sumvecx = vadd(sumvecx, newx);										\
	sumvecy = vadd(sumvecx, newy);


#define v5																\
	sumvecx = vadd(sumvecx, vdiv(theta,pivec));							\
	sumvecy = vadd(sumvecy, vsub(r,onevec));

#define v6																\
	sumvecx = vadd(sumvecx, vmul(vsin(thetaaddr), r));					\
	sumvecy = vadd(sumvecy, vmul(vcos(thetasubr), r));

#define v7																\
	sumvecx = vadd(sumvecx, vmul(vsin(thetamulpi), r));					\
	sumvecy = vsub(sumvecy, vmul(vcos(thetamulpi), r));

#define v8																\
	sumvecx = vadd(sumvecx, vmul(thetadivpi, vsin(pimulr)));			\
	sumvecy = vadd(sumvecy, vmul(thetadivpi, vcos(pimulr)));

#define v9																\
	sumvecx = vadd(sumvecx, vmul(invr, vadd(vcos(theta), vsin(r))));	\
	sumvecy = vadd(sumvecy, vmul(invr, vsub(vsin(theta), vcos(r))));

#define v10																\
	sumvecx = vadd(sumvecx, vdiv(vsin(theta), r));						\
	sumvecy = vadd(sumvecy, vsub(vsin(theta), vcos(r)));


#endif