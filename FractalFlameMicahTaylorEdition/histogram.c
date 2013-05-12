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
#include <intrin.h>

#include "datatypes.h"
#include "rdrand.h"
#include "bmpfile.h"
#include "histogram.h"

// make sure everything is aligned to a 64 bit cache line
_declspec(align(64)) static histocell h[hwid * hhei];
_declspec(align(64)) static omp_lock_t locks[hhei];

_declspec(align(64)) static const __m128 xoffsetvec    = { xoffset, xoffset, xoffset, xoffset };
_declspec(align(64)) static const __m128 yoffsetvec    = { yoffset, yoffset, yoffset, yoffset };

_declspec(align(64)) static const __m128 halfhwidvec   = { hwid/2.0, hwid/2.0, hwid/2.0, hwid/2.0 };
_declspec(align(64)) static const __m128 halfhheivec   = { hhei/2.0, hhei/2.0, hhei/2.0, hhei/2.0 };

_declspec(align(64)) static const __m128 hwidShrunkvec = { hwid/xshrink, hwid/xshrink, hwid/xshrink, hwid/xshrink };
_declspec(align(64)) static const __m128 hheiShrunkvec = { hhei/yshrink, hhei/yshrink, hhei/yshrink, hhei/yshrink };
_declspec(align(64)) static const __m128 halfRGB       =  { 0.5, 0.5, 0.5, 1.0};

void histoinit(){
    int numcells = hwid * hhei;
    size_t histogramSize = numcells * sizeof(histocell);

    memset(h, 0, histogramSize);

    memset(locks, 0, hwid * sizeof(volatile long));
    for(int i = 0; i < hhei; i++)
        omp_init_lock(&(locks[i]));

    if(!h){
        printf("Could not allocate %zu bytes\nPress enter to exit.", hwid * hhei * sizeof(histocell));
        getchar();
        exit(1);
    }

}

static u64 goodHits = 0;
static u64 missHits = 0;
static u64 badHits = 0;

static u64 threadHits[12];

static f128 zerovec = { 0, 0, 0, 0 };

f128tuple histohit(f128tuple xyvec, const colorset pointcolors[4], const i32 th_id){
    if(threadHits[th_id]++ > 20){
        f128 xarr = xyvec.x;
        f128 yarr = xyvec.y;

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
                u32 ix = scaledX.f[i];
                u32 iy = scaledY.f[i];

                if(ix < hwid && iy < hhei){
                    u64 cell = ix + (iy * hwid);
                    // lock the cell
                    //omp_set_lock(&(locks[iy]));

                    __m128 histocolor = vload((float *)&(h[cell]));
                    

                    // add half the new color with the exising color
                    histocolor = vadd(
                                    histocolor,
                                    pointcolors[i].vec);

                    // increment alpha channel
                    //histocolor = vadd(histocolor, incrementAlpha);

                    // write back
                    vstore((float *)&(h[cell]), histocolor);

                    ++goodHits;

                    // unlock the cell
                    //omp_unset_lock(&(locks[iy]));
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

#define max(a,b) (a > b ? a : b) 

void saveimage(){
    printf("Good hits: %llu\t Miss hits: %llu\t Bad hits: %llu\t %f\n", goodHits, missHits, badHits, (f32)goodHits/(badHits > 0 ? badHits : 1));

    bmpfile_t *bmp;

    f32 amax = 1;

    for(u32 i = 0; i < hwid * hhei; i++){
        amax = amax > h[i].a ? amax : h[i].a;
    }

    if((bmp = bmp_create(hwid, hhei, 24)) == NULL){
        printf("Invalid depth value: %d.\n", 24);
        exit(1);
    }

    printf("generating image");

    cilk_for (i32 i = 0; i < hwid * hhei; i++){
        f32 a = log(h[i].a) / log(amax);

        f32 maxColor = max(h[i].r, max(h[i].g, h[i].b));

        if(maxColor <= 0)
            maxColor = 1;

        u8 r = (h[i].r / maxColor) * 0xFF * a;
        u8 g = (h[i].g / maxColor) * 0xFF * a;
        u8 b = (h[i].b / maxColor) * 0xFF * a;

        rgb_pixel_t pixel = {b, g, r, 0xFF};
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
}