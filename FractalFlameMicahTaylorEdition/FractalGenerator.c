#include "datatypes.h"
#include "histogram.h"
#include "rdrand.h"
#include "fastrandsse.h"
#include "variations.h"
#include "vmath.h"

#include <stdio.h>
#include <time.h>
#include <immintrin.h>
#include <Windows.h>

#include <omp.h>

#define n_affine_matrix 6
#define jump_table_size 1024

#define abs(x) (x >= 0 ? x : - x)

typedef struct {
    f32 a, b, c, d, e, f;
    f32 red, green, blue;
} affinematrix;

void affineinit();

static affinematrix am[n_affine_matrix];
static affinematrix * affine_jump_table[jump_table_size];

int main(i32 argc, i8 **argv){
    SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
    srand(time(NULL));

    while(1){
        printf("allocating memory... ");
        histoinit();
        affineinit();
        printf("done\n");

        printf("plotting points\n");
        clock_t start;
        start = clock();

        i32 th_id;

#pragma omp parallel private(th_id)
        {
            th_id = omp_get_thread_num();

            f128 pointvecx;
            f128 pointvecy;
            f128tuple xyvec;

            srand_sse(rdrand_u32(), th_id);

            rand_sse((unsigned int *)pointvecx.f, th_id);
            rand_sse((unsigned int *)pointvecy.f, th_id);

            printf("thread id: %d\n", th_id);
            _sleep(1000);
            for(u32 j = 0; j < 3 * 1000000; j++){
                // seed the random number generator every so often
                srand_sse(rdrand_u32(), th_id);

                if(j % 1000000 == 0){
                    printf("...%u", j/1000000);
                }

                for(u64 i = 0; i < 10; i++){
                    affinematrix *am_itt[4];
                    u32 jumpTable[4];

                    rand_sse(jumpTable, th_id);

                    for(u32 j = 0; j < 4; j++){
                        i32 n = jumpTable[j] % jump_table_size;
                        am_itt[j] = affine_jump_table[n];
                    }

                    // this is slow I think
                    // will be fast with gather instruction on Haswell
                    const __m128 affinea = { am_itt[0]->a, am_itt[1]->a, am_itt[2]->a, am_itt[3]->a };
                    const __m128 affineb = { am_itt[0]->b, am_itt[1]->b, am_itt[2]->b, am_itt[3]->b };
                    const __m128 affinec = { am_itt[0]->c, am_itt[1]->c, am_itt[2]->c, am_itt[3]->c };
                    const __m128 affined = { am_itt[0]->d, am_itt[1]->d, am_itt[2]->d, am_itt[3]->d };
                    const __m128 affinee = { am_itt[0]->e, am_itt[1]->e, am_itt[2]->e, am_itt[3]->e };
                    const __m128 affinef = { am_itt[0]->f, am_itt[1]->f, am_itt[2]->f, am_itt[3]->f };

                    const f128 colorsetr = { am_itt[0]->red,   am_itt[1]->red,   am_itt[2]->red,   am_itt[3]->red   };
                    const f128 colorsetg = { am_itt[0]->green, am_itt[1]->green, am_itt[2]->green, am_itt[3]->green };
                    const f128 colorsetb = { am_itt[0]->blue,  am_itt[1]->blue,  am_itt[2]->blue,  am_itt[3]->blue  };

                    const __m128 twovec = {  2,  2,  2,  2 };
                    const __m128 onevec = {  1,  1,  1,  1 };
                    const __m128 pivec =  { PI, PI, PI, PI };

                    const __m128 affinedx = vadd(
                                                vadd(
                                                    vmul(affinea, pointvecx.v),
                                                    vmul(affineb, pointvecy.v)),
                                                affinec);

                    const __m128 affinedy = vadd(
                                                vadd(
                                                    vmul(affined, pointvecx.v),
                                                    vmul(affinee, pointvecy.v)),
                                                affinef);

                    const __m128 rsq = vadd(
                                        vmul(affinedx, affinedx),
                                        vmul(affinedy, affinedy));

                    const __m128 r = vsqrt(rsq);

                    const __m128 theta = vatan2(affinedx, affinedy);

                    const __m128 thetaaddr = vadd(theta, r);
                    const __m128 thetasubr = vsub(theta, r);
                    const __m128 sinrsq = vsin(rsq);
                    const __m128 cosrsq = vcos(rsq);
                    const __m128 thetamulpi = vmul(theta, r);
                    const __m128 thetadivpi = vdiv(theta, pivec);
                    const __m128 pimulr = vmul(pivec, r);
                    const __m128 invr = vdiv(onevec, r);

                    __m128 sumvecx = {0, 0, 0, 0};
                    __m128 sumvecy = {0, 0, 0, 0};


                    //v1;
                    //v2;
                    //v3;
                    //v4;
                    //v5;
                    //v6;
                    //v7;
                    //v8;
                    //v9;
                    //v10;
                    v11;

                    pointvecx.v = sumvecx;
                    pointvecy.v = sumvecy;

                    xyvec.x = pointvecx;
                    xyvec.y = pointvecy;

                    xyvec = histohit(xyvec, colorsetr, colorsetg, colorsetb, th_id);
                    pointvecx = xyvec.x;
                    pointvecy = xyvec.y;
                }
            }
        }

        clock_t end = clock();

        printf(" done took %f seconds\n", (f64)(end - start)/(f64)CLOCKS_PER_SEC);
        saveimage();

        char *resize_command = "mogrify -format png -path images -resize 1920x1080 fractal.bmp";
        system(resize_command);

    }
    //getchar();
}

void affineinit(){
    // init matrix values
    for(i32 i = 0; i < n_affine_matrix; i++){
        am[i].a = rdrand_f32();
        am[i].b = rdrand_f32();
        am[i].c = rdrand_f32();
        am[i].d = rdrand_f32();
        am[i].e = rdrand_f32();
        am[i].f = rdrand_f32();

        f32 r = rdrand_f32();
        f32 g = rdrand_f32();
        f32 b = rdrand_f32();
        am[i].red = abs(r);
        am[i].green = abs(g);
        am[i].blue = abs(b);
    }

    // init jump table
    for(i32 i = 0; i < jump_table_size; i++){
        affine_jump_table[i] = &(am[rdrand_u32() % n_affine_matrix]);
    }
}