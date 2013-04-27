#ifndef __VARIATIONS_H__
#define __VARIATIONS_H__ __VARIATIONS_H__

#include "vmath.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#endif



#define v1                                                                          \
    sumvecx = vadd(sumvecx, vsin(affinedx));                                        \
    sumvecy = vadd(sumvecy, vcos(affinedy));    

#define v2                                                                          \
    sumvecx = vadd(sumvecx, vdiv(affinedx, rsq));                                   \
    sumvecy = vadd(sumvecy, vdiv(affinedy, rsq));

#define v3                                                                          \
    sumvecx = vadd(sumvecx, vsub(vmul(affinedx, sinrsq), vmul(affinedy, cosrsq)));  \
    sumvecy = vadd(sumvecy, vsub(vmul(affinedx, cosrsq), vmul(affinedy, sinrsq)));


#define v4                                                                          \
    __m128 newx = vdiv(                                                             \
                    vmul(                                                           \
                        vsub(affinedx, affinedy),                                   \
                        vadd(affinedx, affinedy)),                                  \
                    r);                                                             \
                                                                                    \
    __m128 newy = vdiv(                                                             \
                    vmul(                                                           \
                        vmul(affinedx, affinedy),                                   \
                        twovec),                                                    \
                    r);                                                             \
                                                                                    \
    sumvecx = vadd(sumvecx, newx);                                                  \
    sumvecy = vadd(sumvecx, newy);


#define v5                                                                          \
    sumvecx = vadd(sumvecx, vdiv(theta,pivec));                                     \
    sumvecy = vadd(sumvecy, vsub(r,onevec));

#define v6                                                                          \
    sumvecx = vadd(sumvecx, vmul(vsin(thetaaddr), r));                              \
    sumvecy = vadd(sumvecy, vmul(vcos(thetasubr), r));

#define v7                                                                          \
    sumvecx = vadd(sumvecx, vmul(vsin(thetamulpi), r));                             \
    sumvecy = vsub(sumvecy, vmul(vcos(thetamulpi), r));

#define v8                                                                          \
    sumvecx = vadd(sumvecx, vmul(thetadivpi, vsin(pimulr)));                        \
    sumvecy = vadd(sumvecy, vmul(thetadivpi, vcos(pimulr)));

#define v9                                                                          \
    sumvecx = vadd(sumvecx, vmul(invr, vadd(vcos(theta), vsin(r))));                \
    sumvecy = vadd(sumvecy, vmul(invr, vsub(vsin(theta), vcos(r))));

#define v10                                                                         \
    sumvecx = vadd(sumvecx, vdiv(vsin(theta), r));                                  \
    sumvecy = vadd(sumvecy, vsub(vsin(theta), vcos(r)));

#define v11                                                                         \
    sumvecx = vadd(sumvecx, vmul(sintheta, cosr));                                  \
    sumvecy = vadd(sumvecy, vmul(costheta, sinr));

#define v12                                                                         \
    const __m128 p12_0 = vpow(vsin( vadd(theta, r)), threevec);                     \
    const __m128 p12_1 = vpow(vsin( vsub(theta, r)), threevec);                     \
    sumvecx = vadd(sumvecx, vadd(p12_0, p12_1));                                    \
    sumvecy = vadd(sumvecy, vsub(p12_0, p12_1));

#define v13                                                                         \
    __m128 o13_omega;                                                               \
    f32 rand = rdrand_f32();                                                        \
    if(rand < 0.5) {                                                                \
        o13_omega = pivec;                                                          \
    } else {                                                                        \
        o13_omega = zerovec;                                                        \
    }                                                                               \
    sumvecx = vadd(sumvecx, vmul(sqrtr, vcos(vadd(halftheta, o13_omega))));         \
    sumvecy = vadd(sumvecy, vmul(sqrtr, vsin(vadd(halftheta, o13_omega))));

#define v14                                                                         \
    __m128 o14_x;                                                                   \
    __m128 o14_y;                                                                   \
    for(int i14 = 0; i14 < 4; i14++){                                               \
        f32 fx = affinedx.m128_f32[i14];                                            \
        f32 fy = affinedy.m128_f32[i14];                                            \
        if(fx >= 0 && fy >= 0){                                                     \
            o14_x.m128_f32[i14] = fx;                                               \
            o14_y.m128_f32[i14] = fy;                                               \
        } else if (fx < 0 && fy >= 0){                                              \
            o14_x.m128_f32[i14] = 2.0 * fx;                                         \
            o14_y.m128_f32[i14] = fy;                                               \
        } else if(fx >= 0 && fy < 0){                                               \
            o14_x.m128_f32[i14] = fx;                                               \
            o14_y.m128_f32[i14] = fy / 2.0;                                         \
        } else {                                                                    \
            o14_x.m128_f32[i14] = 2.0 * fx;                                         \
            o14_y.m128_f32[i14] = fy / 2.0;                                         \
        }                                                                           \
    }                                                                               \
    sumvecx = vadd(sumvecx, o14_x);                                                 \
    sumvecy = vadd(sumvecy, o14_y);


#define v15                                                                         \
    sumvecx = vadd(sumvecx,                                                         \
                vadd(affinedx,                                                      \
                    vmul(                                                           \
                        affineb,                                                    \
                        vsin(                                                       \
                            vdiv(affinedy, vpow(affinec, twovec))))));              \
    sumvecy = vadd(sumvecy,                                                         \
                vadd(affinedy,                                                      \
                    vmul(                                                           \
                        affinee,                                                    \
                        vsin(                                                       \
                            vdiv(affinedx, vpow(affinef, twovec))))));

#define v16                                                                         \
    sumvecx = vadd(sumvecx, vmul(vdiv(twovec, vadd(r, onevec)), affinedy));         \
    sumvecy = vadd(sumvecy, vmul(vdiv(twovec, vadd(r, onevec)), affinedx))


#define v17                                                                         \
    sumvecx = vadd(sumvecx,                                                         \
                vadd(                                                               \
                    affinedx,                                                       \
                    vmul(                                                           \
                        affinec,                                                    \
                        vsin(                                                       \
                            vtan(                                                   \
                                vmul(                                               \
                                    threevec,                                       \
                                    affinedy))))));                                 \
    sumvecy = vadd(sumvecy,                                                         \
                vadd(                                                               \
                    affinedy,                                                       \
                    vmul(                                                           \
                        affinef,                                                    \
                        vsin(                                                       \
                            vtan(                                                   \
                                vmul(                                               \
                                    threevec,                                       \
                                    affinedx)))))); 

#define v18                                                                         \
    sumvecx = vadd( sumvecx,                                                        \
                vmul(                                                               \
                    vsub(                                                           \
                        affinedx,                                                   \
                        onevec),                                                    \
                    vcos(                                                           \
                        vmul(                                                       \
                            pivec,                                                  \
                            affinedy))));                                           \
    sumvecy = vadd( sumvecy,                                                        \
                vmul(                                                               \
                    vsub(                                                           \
                        affinedx,                                                   \
                        onevec),                                                    \
                    vsin(                                                           \
                        vmul(                                                       \
                            pivec,                                                  \
                            affinedy))));


#endif