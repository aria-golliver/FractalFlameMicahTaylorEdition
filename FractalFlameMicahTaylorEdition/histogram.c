#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "vmath.h"

#include <immintrin.h>
#include <omp.h>
#include <cilk\cilk.h>
#include <cilk\cilk_api_windows.h>
#include <cilk\cilk_stub.h>

#include "datatypes.h"
#include "rdrand.h"
#include "bmpfile.h"
#include "histogram.h"


histocell *h;


__m128 xshrinkvec = { xshrink, xshrink, xshrink, xshrink };
__m128 yshrinkvec = { yshrink, yshrink, yshrink, yshrink };

const __m128 xoffsetvec = { xoffset, xoffset, xoffset, xoffset };
const __m128 yoffsetvec = { yoffset, yoffset, yoffset, yoffset };

const __m128 halfhwidvec = { hwid/2.0, hwid/2.0, hwid/2.0, hwid/2.0 };
const __m128 halfhheivec = { hhei/2.0, hhei/2.0, hhei/2.0, hhei/2.0 };
const __m128 hwidvecvec = { hwid, hwid, hwid, hwid };
const __m128 hheivecvec = { hhei, hhei, hhei, hhei };

const __m128 hwidShrunkvec = { hwid/xshrink, hwid/xshrink, hwid/xshrink, hwid/xshrink };
const __m128 hheiShrunkvec = { hhei/yshrink, hhei/yshrink, hhei/yshrink, hhei/yshrink };

void histoinit(){
    int numcells = hwid * hhei;
    size_t histogramSize = numcells * sizeof(histocell);

    if(h == NULL)
        h = (histocell *) calloc(numcells, sizeof(histocell));
    else
        memset(h, 0, histogramSize);

    if(!h){
        printf("Could not allocate %zu bytes\nPress enter to exit.", hwid * hhei * sizeof(histocell));
        getchar();
        exit(1);
    }

}

u64 goodHits = 0;
u64 missHits = 0;
u64 badHits = 0;

u64 threadHits[12];

typedef union {
    f32 f;
    u32 u;
} f32u32;

f128 zerovec = { 0, 0, 0, 0 };

f128tuple histohit(f128tuple xyvec, const f128 rvec, const f128 gvec, const f128 bvec, const i32 th_id){
    if(threadHits[th_id]++ > 20){
        f128 xarr = xyvec.x;
        f128 yarr = xyvec.y;
        const f128 rarr = rvec;
        const f128 garr = gvec;
        const f128 barr = bvec;

        if(vvalid(xarr.v) && vvalid(yarr.v)){

            f128 scaledX;
            scaledX.v = vadd(
                           vmul(
                               vadd(xarr.v, xoffsetvec),
                               hwidShrunkvec),
                           halfhwidvec);
            f128 scaledY;
            scaledY.v = vadd(
                            vmul(
                                vadd(yarr.v, yoffsetvec),
                                hheiShrunkvec),
                            halfhheivec);

            for (i32 i = 0; i < 4; i++){
                u64 ix = scaledX.f[i];
                u64 iy = scaledY.f[i];

                f32 r = rarr.f[i];
                f32 g = garr.f[i];
                f32 b = barr.f[i];

                u64 cell = ix + (iy * hwid);
                if(ix < hwid && iy < hhei){
                    f32 cellr, cellg, cellb;

                    cellr = h[cell].r;
                    cellg = h[cell].g;
                    cellb = h[cell].b;
                    ++h[cell].a;

                    cellr += r;
                    cellg += g;
                    cellb += b;

                    cellr /= 2;
                    cellg /= 2;
                    cellb /= 2;

                    h[cell].r = cellr;
                    h[cell].g = cellb;
                    h[cell].b = cellg;
                    ++goodHits;
                } else {
                    ++missHits;
                }
            }
        } else {
            ++badHits;
            threadHits[th_id] = 0;
            xarr = zerovec;
            yarr = zerovec;
            xyvec.x = xarr;
            xyvec.y = yarr;
        }
    }
    return xyvec;
}

void saveimage(){
    printf("Good hits: %llu\t Miss hits: %llu\t Bad hits: %llu\t %f\n", goodHits, missHits, badHits, (f32)goodHits/(badHits > 0 ? badHits : 1));

    bmpfile_t *bmp;

    u64 amax = 1;

    for(int i = 0; i < hwid * hhei; i++){
        amax = amax > h[i].a ? amax : h[i].a;
    }

    if((bmp = bmp_create(hwid, hhei, 24)) == NULL){
        printf("Invalid depth value: %d.\n", 24);
        exit(1);
    }

    printf("generating image");

    cilk_for (i32 i = 0; i < hwid * hhei; i++){
        f64 a = (log((f64)h[i].a) / log((f64)amax));

        u8 r = h[i].r * 0xFF * a;
        u8 g = h[i].g * 0xFF * a;
        u8 b = h[i].b * 0xFF * a;

        rgb_pixel_t pixel = {r, g, b, 0xFF};
        u32 x = i % hwid;
        u32 y = i / hwid;
        bmp_set_pixel(bmp, x, y, pixel);
        if((i % ((hwid * hhei) / 20)) == 0)
            printf(".");
    }

    printf(" done\n");

    printf("saving image... ");
    bmp_save(bmp, "fractal.bmp");
    bmp_destroy(bmp);
    printf("done\n");
    printf("\n");
    printf("press enter to quit\n");
}