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

const __m128 halfhwid = { hwid/2.0, hwid/2.0, hwid/2.0, hwid/2.0 };
const __m128 halfhhei = { hhei/2.0, hhei/2.0, hhei/2.0, hhei/2.0 };
const __m128 hwidvec = { hwid, hwid, hwid, hwid };
const __m128 hheivec = { hhei, hhei, hhei, hhei };

const __m128 hwidShrunk = { hwid/xshrink, hwid/xshrink, hwid/xshrink, hwid/xshrink };
const __m128 hheiShrunk = { hhei/yshrink, hhei/yshrink, hhei/yshrink, hhei/yshrink };

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

f128tuple histohit(f128tuple xyvec, const f128 rvec, const f128 gvec, const f128 bvec, const i32 th_id){
    if(threadHits[th_id]++ > 20){
        f128 xarr = xyvec.x;
        f128 yarr = xyvec.y;
        const f128 rarr = rvec;
        const f128 garr = gvec;
        const f128 barr = bvec;

        i32 wasNaN = 0;

        for (i32 i = 0; i < 4; i++){
            f32u32 x;
            x.f = xarr.f[i];
            f32u32 y;
            y.f = yarr.f[i];

            f32 r = rarr.f[i];
            f32 g = garr.f[i];
            f32 b = barr.f[i];
            const u32 NaN = 0x7f800000;
            if(((x.u & NaN) != NaN)  & ((y.u & NaN) != NaN) ){
                u64 ix = (xarr.f[i] + xoffset) * (hwid / xshrink) + (hwid / 2);
                u64 iy = (yarr.f[i] + yoffset) * (hhei / yshrink) + (hhei / 2);
                u64 cell = ix + (iy * hwid);
                if(ix < hwid && iy < hhei){
                    f32 cellr, cellg, cellb;
                    u64 cella;

                    cellr = h[cell].r;
                    cellg = h[cell].g;
                    cellb = h[cell].b;
                    cella = h[cell].a;

                    cellr += r;
                    cellg += g;
                    cellb += b;
                    ++cella;

                    cellr /= 2;
                    cellg /= 2;
                    cellb /= 2;

                    h[cell].r = cellr;
                    h[cell].g = cellb;
                    h[cell].b = cellg;
                    h[cell].a = cella;
                    goodHits++;
                }
            } else {
                missHits++;
                //rand_sse((u32 *) xarr.f, th_id);
                //rand_sse((u32 *) yarr.f, th_id);
                xarr.f[0] = 0;
                xarr.f[1] = 0;
                xarr.f[2] = 0;
                xarr.f[3] = 0;
                yarr.f[0] = 0;
                yarr.f[1] = 0;
                yarr.f[2] = 0;
                yarr.f[3] = 0;
                threadHits[th_id] = 0;
                wasNaN = 1;
            }
        }
        if(wasNaN){
            badHits++;
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