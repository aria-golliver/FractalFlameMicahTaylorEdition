#include "datatypes.h"
#include "histogram.h"
#include "rdrand.h"
#include "variations.h"
#include "vmath.h"
#include "FractalGenerator.h"
#include "fractaldisplay.h"

#include <stdio.h>
#include <time.h>
#include <immintrin.h>
#include <Windows.h>

#include <omp.h>

#define abs(x) (x >= 0 ? x : - x)

char *fractal_name;

_declspec(align(64)) static affinematrix am[n_affine_matrix];
_declspec(align(64)) static affinematrix * affine_jump_table[jump_table_size];
_declspec(align(64)) static f256 variation_weights[MAX_VARIATIONS];
_declspec(align(64)) static f256 parametric_paramaters[MAX_VARIATIONS][4];

int main(i32 argc, i8 **argv){
    SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
    
    #ifdef DISPLAY
        displayinit();
    #endif

    do {
        #ifdef DISPLAY
            displayreset();
        #endif

        u64 fractal_code;
        if(fractal_name){
            free(fractal_name);
            fractal_name = 0;
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
        i32 active_work_threads = 1;

        #pragma omp parallel private(th_id)
        {
            th_id = omp_get_thread_num();

            #pragma omp critical
            {
                active_work_threads = 1;
            }

            #ifdef DISPLAY
                if(th_id == 0){
                    while(active_work_threads){
                        updateDisplay();
                    }
                    goto end;
                }
            #endif

            _declspec(align(64)) f256 pointvecx = { 0, 0, 0, 0, 0, 0, 0, 0 };
            _declspec(align(64)) f256 pointvecy = { 0, 0, 0, 0, 0, 0, 0, 0 };

            _declspec(align(64)) f256tuple xyvec;

            printf("thread id: %d\n", th_id);

            _declspec(align(64)) colorset pointcolors[8];
            memset(pointcolors, 0, sizeof(pointcolors));


            const u64 FLAME_ITTS_multiplier = 1000000;

            for(u64 j = 0; j < FLAME_ITTS * FLAME_ITTS_multiplier; j++){
                if(j % FLAME_ITTS_multiplier == 0){
                    printf("...%u%", j/FLAME_ITTS_multiplier);

                }

                _declspec(align(64)) affinematrix *am_itt[FLOATS_PER_VECTOR_REGISTER];

                typedef union {
                    u32 u32[FLOATS_PER_VECTOR_REGISTER];
                    u64 u64[FLOATS_PER_VECTOR_REGISTER / 2];
                } jumpTable_t;

                _declspec(align(64)) jumpTable_t jumpTable;

                // fill the jump table with random bits
                for(u32 j = 0; j < FLOATS_PER_VECTOR_REGISTER / 2; j++){
                    rdrand_u64(&jumpTable.u64[j]);
                }


                // use the random bits to choose which affine transformation to apply
                for(u32 j = 0; j < FLOATS_PER_VECTOR_REGISTER; j++){
                    i32 n = jumpTable.u32[j] % jump_table_size;
                    am_itt[j] = affine_jump_table[n];
                }

                // load each affine transformation into the vector registers
                _declspec(align(64)) const __m256 affinea = { am_itt[0]->a, am_itt[1]->a, am_itt[2]->a, am_itt[3]->a, am_itt[4]->a, am_itt[5]->a, am_itt[6]->a, am_itt[7]->a };
                _declspec(align(64)) const __m256 affineb = { am_itt[0]->b, am_itt[1]->b, am_itt[2]->b, am_itt[3]->b, am_itt[4]->b, am_itt[5]->b, am_itt[6]->b, am_itt[7]->b };
                _declspec(align(64)) const __m256 affinec = { am_itt[0]->c, am_itt[1]->c, am_itt[2]->c, am_itt[3]->c, am_itt[4]->c, am_itt[5]->c, am_itt[6]->c, am_itt[7]->c };
                _declspec(align(64)) const __m256 affined = { am_itt[0]->d, am_itt[1]->d, am_itt[2]->d, am_itt[3]->d, am_itt[4]->d, am_itt[5]->d, am_itt[6]->d, am_itt[7]->d };
                _declspec(align(64)) const __m256 affinee = { am_itt[0]->e, am_itt[1]->e, am_itt[2]->e, am_itt[3]->e, am_itt[4]->e, am_itt[5]->e, am_itt[6]->e, am_itt[7]->e };
                _declspec(align(64)) const __m256 affinef = { am_itt[0]->f, am_itt[1]->f, am_itt[2]->f, am_itt[3]->f, am_itt[4]->f, am_itt[5]->f, am_itt[6]->f, am_itt[7]->f };
                
                // apply the affine transformation to the existing points
                _declspec(align(64)) const __m256 affinedx = vadd(
                                                                vadd(
                                                                    vmul(affinea, pointvecx.v),
                                                                    vmul(affineb, pointvecy.v)),
                                                                affinec);

                _declspec(align(64)) const __m256 affinedy = vadd(
                                                                vadd(
                                                                    vmul(affined, pointvecx.v),
                                                                    vmul(affinee, pointvecy.v)),
                                                                affinef);
                // each affine transformation has a color associated with it, update the current color with the new one (old + new)/2
                __m128 colorset;
                for(int c = 0; c < FLOATS_PER_VECTOR_REGISTER; c++){
                    __m128 newcolor = am_itt[c]->color.vec;
                    pointcolors[c].vec = vmul128(vadd128(pointcolors[c].vec, newcolor), halfvec128);
                }

                _declspec(align(64)) const __m256 rsq = vadd(
                                                    vmul(affinedx, affinedx),
                                                    vmul(affinedy, affinedy));
                
                _declspec(align(64)) const __m256 r = vsqrt(rsq);
                
                _declspec(align(64)) const __m256 theta = vatan2(affinedx, affinedy);
                
                _declspec(align(64)) const __m256 thetaaddr = vadd(theta, r);
                _declspec(align(64)) const __m256 thetasubr = vsub(theta, r);
                _declspec(align(64)) const __m256 sinrsq = vsin(rsq);
                _declspec(align(64)) const __m256 cosrsq = vcos(rsq);
                _declspec(align(64)) const __m256 sintheta = vsin(theta);
                _declspec(align(64)) const __m256 costheta = vcos(theta);
                _declspec(align(64)) const __m256 sinr = vsin(r);
                _declspec(align(64)) const __m256 cosr = vcos(r);
                _declspec(align(64)) const __m256 thetamulpi = vmul(theta, r);
                _declspec(align(64)) const __m256 thetadivpi = vdiv(theta, pivec);
                _declspec(align(64)) const __m256 pimulr = vmul(pivec, r);
                _declspec(align(64)) const __m256 invr = vdiv(onevec, r);
                _declspec(align(64)) const __m256 sqrtr = vsqrt(r);
                _declspec(align(64)) const __m256 halftheta = vdiv(theta, twovec);
                
                _declspec(align(64)) __m256 sumvecx = { 0, 0, 0, 0, 0, 0, 0, 0 };
                _declspec(align(64)) __m256 sumvecy = { 0, 0, 0, 0, 0, 0, 0, 0 };
                //theta.m256_f32
                //v0;
                //v1;
                //v2;
                //v3;
                //v4;
                //v5;
                //v6;
                v7;
                v8;
                //v9;
                //v10;
                //v11;
                //v12;
                //v13;
                //v14;
                //v15;
                //v16;
                v17;
                //v18;
                //v19;
                //v20;
                //v21;
                v22;
                //v25;
                //v28;
                //v30;
                
                xyvec.x.v = sumvecx;
                xyvec.y.v = sumvecy;

                xyvec = histohit(xyvec, pointcolors, th_id);
                pointvecx = xyvec.x;
                pointvecy = xyvec.y;
            }
        
#pragma omp critical
        {
            active_work_threads = 0;
        }
end:;
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


#define MAX(a,b) (a > b ? a : b) 
#define MAX3(a,b,c) MAX(a, MAX(b, c))

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
        f32 g;
        f32 b;
        rdrand_f32(&r);
        rdrand_f32(&g);
        rdrand_f32(&b);
        
        r = abs(r);
        g = abs(g);
        b = abs(b);

        f32 maxColor = 1;MAX3(r,g,b);

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

/*
 * this should be updated so it only accounts for the variations which are selected
 * instead of accounting for all of them
 */
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
    /*
    if(total < 1.0){
        for(u32 i = 0; i < MAX_VARIATIONS; i++){
            variation_weights[i].f[0] /= total;
            variation_weights[i].f[1] /= total;
            variation_weights[i].f[2] /= total;
            variation_weights[i].f[3] /= total;
        }
    }
    */
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
        fprintf(file, "\tam[%u].a = %.20ff;\n", i, am[i].a);
        fprintf(file, "\tam[%u].b = %.20ff;\n", i, am[i].b);
        fprintf(file, "\tam[%u].c = %.20ff;\n", i, am[i].c);
        fprintf(file, "\tam[%u].d = %.20ff;\n", i, am[i].d);
        fprintf(file, "\tam[%u].e = %.20ff;\n", i, am[i].e);
        fprintf(file, "\tam[%u].f = %.20ff;\n", i, am[i].f);
        fprintf(file, "\tam[%u].color.r = %.20ff;\n", i, am[i].color.r);
        fprintf(file, "\tam[%u].color.g = %.20ff;\n", i, am[i].color.g);
        fprintf(file, "\tam[%u].color.b = %.20ff;\n", i, am[i].color.b);
        fprintf(file, "\tam[%u].color.a = %.20ff;\n", i, am[i].color.a);
    }

    for(u32 i = 0; i < jump_table_size; i++){
        fprintf(file, "\taffine_jump_table[%u] = &(am[%u]);\n", i, (u32)( affine_jump_table[i] - am));
    }
    fprintf(file, "}\n");

    fprintf(file, "void affineinit(){\n");
    for(u32 i = 0; i < MAX_VARIATIONS; i++){
        fprintf(file, "\tvariation_weights[%u].f[0] = %.20ff;\n", i, variation_weights[i].f[0]);
        fprintf(file, "\tvariation_weights[%u].f[1] = %.20ff;\n", i, variation_weights[i].f[1]);
        fprintf(file, "\tvariation_weights[%u].f[2] = %.20ff;\n", i, variation_weights[i].f[2]);
        fprintf(file, "\tvariation_weights[%u].f[3] = %.20ff;\n", i, variation_weights[i].f[3]);
    }

    for(u32 i = 0; i < MAX_VARIATIONS; i++){
        for(u32 j = 0; j < 4; j++){
            //f32 r = rdrand_f32(&r);
            for(u32 k = 0; k < 4; k++){
                fprintf(file, "\tparametric_paramaters[%u][%u].f[%u] = %.20f;\n", i, j, k, parametric_paramaters[i][j].f[k]);
            }
        }
    }
    
    fprintf(file, "}\n");
    fprintf(file, "#endif\n");
    fprintf(file, "#endif\n");
    fclose(file);
    
}
