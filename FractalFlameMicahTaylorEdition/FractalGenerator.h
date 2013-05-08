#ifndef __FRACTAL_GENERATOR_H
#define __FRACTAL_GENERATOR_H __FRACTAL_GENERATOR_H

#include "datatypes.h"

typedef struct {
    f32 a, b, c, d, e, f;
    f32 red, green, blue;
} affinematrix;

void affineinit();
void variationinit();
void compressimage();
void savegenome();

#define n_affine_matrix (3)
#define jump_table_size (1000)
#define MAX_VARIATIONS (50)

#define FLAME_ITTS (2)
#define RUN_FOREVER (1)

#endif