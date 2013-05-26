#ifndef __FRACTAL_GENERATOR_H
#define __FRACTAL_GENERATOR_H __FRACTAL_GENERATOR_H

#include "datatypes.h"

void affineinit();
void variationinit();
void compressimage();
void savegenome();

#define n_affine_matrix (6)
#define jump_table_size (1000)
#define MAX_VARIATIONS (50)

#define FLAME_ITTS (1)
#define RUN_FOREVER (1)


static const __m128 zerovec   = { 0, 0, 0, 0 };
static const __m128 onevec    = { 1, 1, 1, 1 };
static const __m128 twovec    = { 2, 2, 2, 2 };
static const __m128 threevec  = { 3, 3, 3, 3 };
static const __m128 fourvec   = { 4, 4, 4, 4 };
static const __m128 pivec     = { PI, PI, PI, PI };
static const __m128 negonevec = { -1, -1, -1, -1 };
static const __m128 halfvec   = { 0.5, 0.5, 0.5, 0.5 };

#endif