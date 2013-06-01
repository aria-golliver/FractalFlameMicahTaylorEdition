#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "vmath.h"

#include <immintrin.h>
#include <omp.h>
#include <cilk\cilk.h>
#include <cilk\cilk_api.h>
#include <cilk\cilk_api_windows.h>
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

/*
 * initializes the histogram, allocates memory and zeros it
 * also initializes the locks, each lock covers one horizontal row of the histogram
 */
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
        exit(EXIT_FAILURE);
    }

}

static u64 goodHits = 0;
static u64 missHits = 0;
static u64 badHits = 0;

_declspec(align(64)) static u64 threadHits[12];

_declspec(align(64)) static f128 zerovec = { 0, 0, 0, 0 };

f128tuple histohit(f128tuple xyvec, const colorset pointcolors[4], const i32 th_id){
    // don't plot the first 20 iterations
    if(threadHits[th_id]++ > 20){
        f128 xarr = xyvec.x;
        f128 yarr = xyvec.y;

        // check to see if any points have escaped to infinity or are NaN
        // if they have, reset them and the threadHits counter to 0
        if(vvalid(xarr.v) && vvalid(yarr.v)){

            // scale the points from fractal-space to histogram-space
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
            
            // extract each point and plot
            for (i32 i = 0; i < FLOATS_PER_VECTOR_REGISTER; i++){
                const u32 ix = scaledX.f[i];
                const u32 iy = scaledY.f[i];

                if(ix < hwid && iy < hhei){
                    const u64 cell = ix + (iy * hwid);
                    // lock the row
                    omp_set_lock(&(locks[iy]));
                    
                    // load the existing data
                    // the cache miss here takes up maybe 2/3s of the program's execution time
                    __m128 histocolor = vload((float *)&(h[cell]));
                    

                    // add the new color to the old
                    histocolor = vadd(
                                    histocolor,
                                    pointcolors[i].vec);

                    // write back
                    vstore((float *)&(h[cell]), histocolor);

                    ++goodHits;

                    // unlock the cell
                    omp_unset_lock(&(locks[iy]));
                } else {
                    ++missHits;
                }
            }
        } else {
            ++badHits;
            threadHits[th_id] = 0;
            xyvec.x = zerovec;
            xyvec.y = zerovec;
        }
    }
    return xyvec;
}

#define MAX(a,b) (a > b ? a : b) 
#define MAX3(a,b,c) MAX(a, MAX(b, c))

void saveimage(){
    printf("Good hits: %llu\t Miss hits: %llu\t Bad hits: %llu\t %f\n", goodHits, missHits, badHits, (f32)goodHits/(badHits > 0 ? badHits : 1));

    bmpfile_t *bmp;

    f32 amax = 1;

    for(u32 i = 0; i < hwid * hhei; i++){
        amax = MAX(amax, h[i].a);
    }

    bmp = bmp_create(hwid, hhei, 24);

    printf("generating image");

    cilk_for(i32 i = 0; i < hwid * hhei; i++){
        const f32 a = log(h[i].a) / log(amax);

        f32 maxColor = MAX3(h[i].r, h[i].g, h[i].b);

        if(maxColor <= 0)
            maxColor = 1;

        const u8 r = (h[i].r / h[i].a) * 0xFF * a;
        const u8 g = (h[i].g / h[i].a) * 0xFF * a;
        const u8 b = (h[i].b / h[i].a) * 0xFF * a;

        const rgb_pixel_t pixel = {b, g, r, 0xFF};
        const u32 x = i % hwid;
        const u32 y = i / hwid;
        bmp_set_pixel(bmp, x, y, pixel);

        // progress bar
        if((i % ((hwid * hhei) / 20)) == 0)
            printf(".");
    }

    printf(" done\n");

    printf("saving image... ");
    bmp_save(bmp, "fractal.bmp");
    bmp_destroy(bmp);
    printf("done\n");
}

histocell histoget(u64 cell){
    return h[cell];
};
