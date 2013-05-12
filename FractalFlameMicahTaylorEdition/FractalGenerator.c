#include "datatypes.h"
#include "histogram.h"
#include "rdrand.h"
#include "variations.h"
#include "vmath.h"
#include "FractalGenerator.h"

#include <stdio.h>
#include <time.h>
#include <immintrin.h>
#include <Windows.h>

#include <omp.h>

#define abs(x) (x >= 0 ? x : - x)

char *fractal_name;

_declspec(align(64)) static affinematrix am[n_affine_matrix];
_declspec(align(64)) static affinematrix * affine_jump_table[jump_table_size];
_declspec(align(64)) static f128 variation_weights[MAX_VARIATIONS];
_declspec(align(64)) static f128 parametric_paramaters[MAX_VARIATIONS][4];

int main(i32 argc, i8 **argv){
    SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
    u64 dl;
    rdrand_u64(&dl);
    i32 dd;
    rdrand_i32(&dd);
    f32 df;
    rdrand_f32(&df);
    printf("%d, %ull, %f\n", dd, dl, df);

    do {
        u64 fractal_code;
        if(fractal_name){
            free(fractal_name);
        }
        fractal_name = (char *) calloc(1024, sizeof(char));
        rdrand_u64(&fractal_code);
        sprintf(fractal_name, "%llu", fractal_code);

        printf("Creating fractal named: %s\n", fractal_name);
        printf("allocating memory... ");
        histoinit();
        affineinit();
        variationinit();


        printf("done\n");

        printf("plotting points\n");
        clock_t start;
        start = clock();

        i32 th_id;

#pragma omp parallel private(th_id)
        {
            th_id = omp_get_thread_num();

            f128 pointvecx = { 0, 0, 0, 0 };
            f128 pointvecy = { 0, 0, 0, 0 };

            f128tuple xyvec;

            u32 tmp;

            printf("thread id: %d\n", th_id);
            Sleep(1000);

            colorset pointcolors[4];
            memset(pointcolors, 0, sizeof(pointcolors));

            for(u64 j = 0; j < FLAME_ITTS * 10000000ull; j++){

                if(j % 10000000ull == 0){
                    printf("...%u%", j/10000000ull);
                }

                affinematrix *am_itt[4];
                u32 jumpTable[4];

                rdrand_u32(&jumpTable[0]);
                rdrand_u32(&jumpTable[1]);
                rdrand_u32(&jumpTable[2]);
                rdrand_u32(&jumpTable[3]);

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
                
                __m128 colorset;
                for(int c = 0; c < 4; c++){
                    __m128 newcolor = am_itt[c]->color.vec;
                    pointcolors[c].vec = vmul(vadd(pointcolors[c].vec, newcolor), halfvec);
                }
                
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
                const __m128 sintheta = vsin(theta);
                const __m128 costheta = vcos(theta);
                const __m128 sinr = vsin(r);
                const __m128 cosr = vcos(r);
                const __m128 thetamulpi = vmul(theta, r);
                const __m128 thetadivpi = vdiv(theta, pivec);
                const __m128 pimulr = vmul(pivec, r);
                const __m128 invr = vdiv(onevec, r);
                const __m128 sqrtr = vsqrt(r);
                const __m128 halftheta = vdiv(theta, twovec);
                
                __m128 sumvecx = { 0, 0, 0, 0 };
                __m128 sumvecy = { 0, 0, 0, 0 };

                //v1;
                //v2;
                //v3;
                v4;
                //v5;
                //v6;
                //v7;
                //v8;
                //v9;
                //v10;
                //v11;
                //v12;
                //v13;
                //v14;
                //v15;
                //v16;
                //v17;
                //v18;
                //v19;
                //v20;
                //v21;
                //v22;
                //v25;
                //v28;
                //v30;
                
                xyvec.x.v = sumvecx;
                xyvec.y.v = sumvecy;

                xyvec = histohit(xyvec, pointcolors, th_id);
                pointvecx = xyvec.x;
                pointvecy = xyvec.y;
            }
        }

        clock_t end = clock();

        printf("\n done took %f seconds\n", (f64)(end - start)/(f64)CLOCKS_PER_SEC);
        saveimage();

        printf("compressing image... ");
        compressimage();
        printf("done\n");
        printf("saving genome... ");
        savegenome();
        printf("done\n");
        printf("-------------------------------\n");
    } while(RUN_FOREVER);
}

void compressimage(){
        char *rand_filename = (char *) calloc(1024, sizeof(char));

        char *base_command = "mogrify -format png -path images -resize 1920x1080 -write images/";
        char *extension = ".png ";
        u64 tmp;
        strcat(rand_filename, base_command);
        strcat(rand_filename, fractal_name);
        strcat(rand_filename, extension);
        strcat(rand_filename, " fractal.bmp");

        system(rand_filename);
        free(rand_filename);
}

#include "fractalgenome.h"
#ifndef USING_GENOME

void affineinit(){
    // init matrix values
    for(u32 i = 0; i < n_affine_matrix; i++){
        rdrand_f32(&(am[i].a));
        rdrand_f32(&(am[i].b));
        rdrand_f32(&(am[i].c));
        rdrand_f32(&(am[i].d));
        rdrand_f32(&(am[i].e));
        rdrand_f32(&(am[i].f));

        f32 r;
        rdrand_f32(&r);
        f32 g;
        rdrand_f32(&g);
        f32 b;
        rdrand_f32(&b);
        
        r = abs(r);
        g = abs(g);
        b = abs(b);

        f32 maxColor = max(r,max(g,b));

        am[i].color.r  = r / maxColor;
        am[i].color.g  = g / maxColor;
        am[i].color.b  = b / maxColor;
        am[i].color.a  = 1;
    }

    // init jump table
    for(i32 i = 0; i < jump_table_size; i++){
        u32 tmp;
        rdrand_u32(&tmp);
        affine_jump_table[i] = &(am[tmp % n_affine_matrix]);
    }
}

void variationinit(){
    float total = 0;
    for(u32 i = 0; i < MAX_VARIATIONS; i++){
        float weight;
        rdrand_f32(&weight);
        weight = abs(weight);
        variation_weights[i].f[0] = weight;
        variation_weights[i].f[1] = weight;
        variation_weights[i].f[2] = weight;
        variation_weights[i].f[3] = weight;
        total += weight;
    }

    if(total < 1.0){
        for(u32 i = 0; i < MAX_VARIATIONS; i++){
            variation_weights[i].f[0] /= total;
            variation_weights[i].f[1] /= total;
            variation_weights[i].f[2] /= total;
            variation_weights[i].f[3] /= total;
        }
    }

    for(u32 i = 0; i < MAX_VARIATIONS; i++){
        for(u32 j = 0; j < 4; j++){
            f32 r;
            rdrand_f32(&r);
            for(u32 k = 0; k < 4; k++){
                parametric_paramaters[i][j].f[k] = r;
            }
        }
    }
}

#endif

void savegenome(){
    
    FILE *file;
    char *genome_filename = (char *) calloc(1024, sizeof(char));
    strcat(genome_filename, "images/");
    strcat(genome_filename, fractal_name);
    strcat(genome_filename, ".fractalgenome");
    file = fopen(genome_filename, "w");
    free(genome_filename);

    fprintf(file, "#ifndef __FRACTALGENOME_H__\n\
#define __FRACTALGENOME_H__ __FRACTALGENOME_H__\n\
\
// comment this line out to use a random genome\n\
#define USING_GENOME 1\n\
\
#ifdef USING_GENOME\n\n");

    fprintf(file, "void variationinit(){\n");
    for(u32 i = 0; i < n_affine_matrix; i++){
        affinematrix matrix = am[i];
        fprintf(file, "\tam[%d].a = %.20ff;\n", i, am[i].a);
        fprintf(file, "\tam[%d].b = %.20ff;\n", i, am[i].b);
        fprintf(file, "\tam[%d].c = %.20ff;\n", i, am[i].c);
        fprintf(file, "\tam[%d].d = %.20ff;\n", i, am[i].d);
        fprintf(file, "\tam[%d].e = %.20ff;\n", i, am[i].e);
        fprintf(file, "\tam[%d].f = %.20ff;\n", i, am[i].f);
        fprintf(file, "\tam[%d].color.r = %.20ff;\n", i, am[i].color.r);
        fprintf(file, "\tam[%d].color.g = %.20ff;\n", i, am[i].color.g);
        fprintf(file, "\tam[%d].color.b = %.20ff;\n", i, am[i].color.b);
        fprintf(file, "\tam[%d].color.a = %.20ff;\n", i, am[i].color.a);
    }

    for(u32 i = 0; i < jump_table_size; i++){
        fprintf(file, "\taffine_jump_table[%d] = &(am[%d]);\n", i, affine_jump_table[i] - am);
    }
    fprintf(file, "}\n");

    fprintf(file, "void affineinit(){\n");
    for(u32 i = 0; i < MAX_VARIATIONS; i++){
        fprintf(file, "\tvariation_weights[%d].f[0] = %.20ff;\n", i, variation_weights[i].f[0]);
        fprintf(file, "\tvariation_weights[%d].f[1] = %.20ff;\n", i, variation_weights[i].f[1]);
        fprintf(file, "\tvariation_weights[%d].f[2] = %.20ff;\n", i, variation_weights[i].f[2]);
        fprintf(file, "\tvariation_weights[%d].f[3] = %.20ff;\n", i, variation_weights[i].f[3]);
    }

    for(u32 i = 0; i < MAX_VARIATIONS; i++){
        for(u32 j = 0; j < 4; j++){
            //f32 r = rdrand_f32(&r);
            for(u32 k = 0; k < 4; k++){
                fprintf(file, "\tparametric_paramaters[%d][%d].f[%d] = %.20f;\n", i, j, k, parametric_paramaters[i][j].f[k]);
            }
        }
    }
    
    fprintf(file, "}\n");
    fprintf(file, "#endif\n");
    fprintf(file, "#endif\n");
    fclose(file);
    
}
