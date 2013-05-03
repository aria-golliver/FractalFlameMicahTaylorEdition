#ifndef __VARIATIONS_H__
#define __VARIATIONS_H__ __VARIATIONS_H__

#include "vmath.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#endif



#define v1                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vsin(affinedx), variation_weights[1].v));                                          \
    sumvecy = vadd(sumvecy, vmul(vcos(affinedy), variation_weights[1].v));    

#define v2                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vdiv(affinedx, rsq), variation_weights[2].v));                                     \
    sumvecy = vadd(sumvecy, vmul(vdiv(affinedy, rsq), variation_weights[2].v));
                            
#define v3                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vsub(vmul(affinedx, sinrsq), vmul(affinedy, cosrsq)), variation_weights[3].v));    \
    sumvecy = vadd(sumvecy, vmul(vsub(vmul(affinedx, cosrsq), vmul(affinedy, sinrsq)), variation_weights[3].v));


#define v4                                                                                                          \
    __m128 newx = vmul(vdiv(                                                                                        \
                    vmul(                                                                                           \
                        vsub(affinedx, affinedy),                                                                   \
                        vadd(affinedx, affinedy)),                                                                  \
                    r), variation_weights[4].v);                                                                    \
                                                                                                                    \
    __m128 newy = vmul(vdiv(                                                                                        \
                    vmul(                                                                                           \
                        vmul(affinedx, affinedy),                                                                   \
                        twovec),                                                                                    \
                    r), variation_weights[4].v);                                                                    \
                                                                                                                    \
    sumvecx = vadd(sumvecx, newx);                                                                                  \
    sumvecy = vadd(sumvecx, newy);


#define v5                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vdiv(theta,pivec), variation_weights[5].v));                                       \
    sumvecy = vadd(sumvecy, vmul(vsub(r,onevec), variation_weights[5].v));

#define v6                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vmul(vsin(thetaaddr), r), variation_weights[6].v));                                \
    sumvecy = vadd(sumvecy, vmul(vmul(vcos(thetasubr), r), variation_weights[6].v));

#define v7                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vmul(vsin(thetamulpi), r), variation_weights[7].v));                               \
    sumvecy = vsub(sumvecy, vmul(vmul(vcos(thetamulpi), r), variation_weights[7].v));

#define v8                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vmul(thetadivpi, vsin(pimulr)), variation_weights[8].v));                          \
    sumvecy = vadd(sumvecy, vmul(vmul(thetadivpi, vcos(pimulr)), variation_weights[8].v));

#define v9                                                                                                          \
    sumvecx = vadd(sumvecx, vmul(vmul(invr, vadd(vcos(theta), vsin(r))), variation_weights[9].v));                  \
    sumvecy = vadd(sumvecy, vmul(vmul(invr, vsub(vsin(theta), vcos(r))), variation_weights[9].v));

#define v10                                                                                                         \
    sumvecx = vadd(sumvecx, vmul(vdiv(vsin(theta), r), variation_weights[10].v));                                   \
    sumvecy = vadd(sumvecy, vmul(vsub(vsin(theta), vcos(r)), variation_weights[10].v));

#define v11                                                                                                         \
    sumvecx = vadd(sumvecx, vmul(vmul(sintheta, cosr), variation_weights[11].v));                                   \
    sumvecy = vadd(sumvecy, vmul(vmul(costheta, sinr), variation_weights[11].v));

#define v12                                                                                                         \
    const __m128 p12_0 = vpow(vsin( vadd(theta, r)), threevec);                                                     \
    const __m128 p12_1 = vpow(vsin( vsub(theta, r)), threevec);                                                     \
    sumvecx = vadd(sumvecx, vmul(vadd(p12_0, p12_1), variation_weights[12].v));                                     \
    sumvecy = vadd(sumvecy, vmul(vsub(p12_0, p12_1), variation_weights[12].v));

#define v13                                                                                                         \
    __m128 o13_omega;                                                                                               \
    f32 rand = rdrand_f32(&rand);                                                                                   \
    if(rand < 0.5) {                                                                                                \
        o13_omega = pivec;                                                                                          \
    } else {                                                                                                        \
        o13_omega = zerovec;                                                                                        \
    }                                                                                                               \
    sumvecx = vadd(sumvecx, vmul(vmul(sqrtr, vcos(vadd(halftheta, o13_omega))), variation_weights[13].v));          \
    sumvecy = vadd(sumvecy, vmul(vmul(sqrtr, vsin(vadd(halftheta, o13_omega))), variation_weights[13].v));

#define v14                                                                                                         \
    __m128 o14_x = { 0, 0, 0, 0 };                                                                                  \
    __m128 o14_y = { 0, 0, 0, 0 };                                                                                  \
    for(int i14 = 0; i14 < 4; i14++){                                                                               \
        f32 fx = affinedx.m128_f32[i14];                                                                            \
        f32 fy = affinedy.m128_f32[i14];                                                                            \
        if(fx >= 0 && fy >= 0){                                                                                     \
            o14_x.m128_f32[i14] = fx;                                                                               \
            o14_y.m128_f32[i14] = fy;                                                                               \
        } else if (fx < 0 && fy >= 0){                                                                              \
            o14_x.m128_f32[i14] = 2.0 * fx;                                                                         \
            o14_y.m128_f32[i14] = fy;                                                                               \
        } else if(fx >= 0 && fy < 0){                                                                               \
            o14_x.m128_f32[i14] = fx;                                                                               \
            o14_y.m128_f32[i14] = fy / 2.0;                                                                         \
        } else {                                                                                                    \
            o14_x.m128_f32[i14] = 2.0 * fx;                                                                         \
            o14_y.m128_f32[i14] = fy / 2.0;                                                                         \
        }                                                                                                           \
    }                                                                                                               \
    sumvecx = vadd(sumvecx, vmul(o14_x, variation_weights[14].v));                                                  \
    sumvecy = vadd(sumvecy, vmul(o14_y, variation_weights[14].v));


#define v15                                                                                                         \
    sumvecx = vadd(sumvecx,                                                                                         \
                vmul(                                                                                               \
                vadd(affinedx,                                                                                      \
                    vmul(                                                                                           \
                        affineb,                                                                                    \
                        vsin(                                                                                       \
                            vdiv(affinedy, vpow(affinec, twovec))))), variation_weights[15].v));                    \
    sumvecy = vadd(sumvecy,                                                                                         \
                vmul(                                                                                               \
                vadd(affinedy,                                                                                      \
                    vmul(                                                                                           \
                        affinee,                                                                                    \
                        vsin(                                                                                       \
                            vdiv(affinedx, vpow(affinef, twovec))))), variation_weights[15].v));

#define v16                                                                                                         \
    sumvecx = vadd(sumvecx, vmul(vmul(vdiv(twovec, vadd(r, onevec)), affinedy), variation_weights[16].v));          \
    sumvecy = vadd(sumvecy, vmul(vmul(vdiv(twovec, vadd(r, onevec)), affinedx), variation_weights[16].v));


#define v17                                                                                                         \
    sumvecx = vadd(sumvecx,                                                                                         \
                vmul(                                                                                               \
                vadd(                                                                                               \
                    affinedx,                                                                                       \
                    vmul(                                                                                           \
                        affinec,                                                                                    \
                        vsin(                                                                                       \
                            vtan(                                                                                   \
                                vmul(                                                                               \
                                    threevec,                                                                       \
                                    affinedy))))), variation_weights[17].v));                                       \
    sumvecy = vadd(sumvecy,                                                                                         \
                vmul(                                                                                               \
                vadd(                                                                                               \
                    affinedy,                                                                                       \
                    vmul(                                                                                           \
                        affinef,                                                                                    \
                        vsin(                                                                                       \
                            vtan(                                                                                   \
                                vmul(                                                                               \
                                    threevec,                                                                       \
                                    affinedx))))), variation_weights[17].v)); 

#define v18                                                                                                         \
    sumvecx = vadd( sumvecx,                                                                                        \
                vmul(                                                                                               \
                vmul(                                                                                               \
                    vsub(                                                                                           \
                        affinedx,                                                                                   \
                        onevec),                                                                                    \
                    vcos(                                                                                           \
                        vmul(                                                                                       \
                            pivec,                                                                                  \
                            affinedy))), variation_weights[18].v));                                                 \
    sumvecy = vadd( sumvecy,                                                                                        \
                vmul(                                                                                               \
                vmul(                                                                                               \
                    vsub(                                                                                           \
                        affinedx,                                                                                   \
                        onevec),                                                                                    \
                    vsin(                                                                                           \
                        vmul(                                                                                       \
                            pivec,                                                                                  \
                            affinedy))), variation_weights[18].v));

#define v19                                                                                                         \
    sumvecx = vadd(sumvecx,                                                                                         \
                vmul(                                                                                               \
                vmul(                                                                                               \
                    vpow(                                                                                           \
                        r,                                                                                          \
                        vsin(theta)),                                                                               \
                    vcos(theta)), variation_weights[19].v));                                                        \
    sumvecy = vadd(sumvecy,                                                                                         \
            vmul(                                                                                                   \
            vmul(                                                                                                   \
                vpow(                                                                                               \
                    r,                                                                                              \
                    vsin(theta)),                                                                                   \
                vsin(theta)), variation_weights[19].v));

// bugged
#define v20                                                                                                         \
    sumvecx = vadd(sumvecx,                                                                                         \
                vmul(                                                                                               \
                vmul(                                                                                               \
                    vcos(                                                                                           \
                        vmul(                                                                                       \
                            pivec,                                                                                  \
                            affinedx)),                                                                             \
                    vcosh(affinedy)), variation_weights[20].v));                                                    \
                                                                                                                    \
    sumvecy = vadd(sumvecy,                                                                                         \
                vmul(                                                                                               \
                vmul( negonevec,                                                                                    \
                    vmul(                                                                                           \
                        vsin(                                                                                       \
                            vmul(                                                                                   \
                                pivec,                                                                              \
                                affinedx)),                                                                         \
                        vsinh(affinedy))), variation_weights[20].v));

#define v21                                                                                                         \
    const __m128 csquared = vmul(affinec, affinec);                                                                 \
    const __m128 rpluscsquared = vadd(r, csquared);                                                                 \
    const __m128 multiplier =                                                                                       \
                    vadd(                                                                                           \
                        vsub(                                                                                       \
                            vmod(rpluscsquared, vmul(twovec, csquared)),                                            \
                            csquared),                                                                              \
                        vmul(                                                                                       \
                            r,                                                                                      \
                            vsub(onevec, csquared)));                                                               \
    sumvecx = vadd(sumvecx, vmul(vmul(multiplier, vcos(theta)), variation_weights[21].v));                          \
    sumvecy = vadd(sumvecy, vmul(vmul(multiplier, vsin(theta)), variation_weights[21].v));

#define v22                                                                                                         \
    f128 t22;                                                                                                       \
    t22.v = vmul(pivec, vmul(affinee, affinee));                                                                    \
    const __m128 halft22 = vdiv(t22.v, twovec);                                                                     \
    f128 switchvec;                                                                                                 \
    switchvec.v = vmod(vadd(theta, affinef), t22.v);                                                                \
    f128 tvec;                                                                                                      \
    for(int temp22 = 0; temp22 < 4; temp22++){                                                                      \
        if(switchvec.f[temp22] > halft22.f[temp22]){                                                                \
            tvec.f[temp22] = theta.f[temp22] - halft22.f[temp22];                                                   \
        } else {                                                                                                    \
            tvec.f[temp22] = theta.f[temp22] + halft22.f[temp22];                                                   \
        }                                                                                                           \
    }                                                                                                               \
    sumvecx = vadd(sumvecx,                                                                                         \
                vmul(                                                                                               \
                    vmul(                                                                                           \
                        r,                                                                                          \
                        vcos(tvec.v)),                                                                              \
                    variation_weights[22].v));                                                                      \
    sumvecy = vadd(sumvecy,                                                                                         \
            vmul(                                                                                                   \
                vmul(                                                                                               \
                    r,                                                                                              \
                    vsin(tvec.v)),                                                                                  \
                variation_weights[22].v));
    
#endif