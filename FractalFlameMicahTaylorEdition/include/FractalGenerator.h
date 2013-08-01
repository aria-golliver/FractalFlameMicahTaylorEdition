#ifndef __FRACTAL_GENERATOR_H
#define __FRACTAL_GENERATOR_H __FRACTAL_GENERATOR_H

#include "datatypes.h"

void affineinit();
void variationinit();
void compressimage();
void savegenome();

#define n_affine_matrix (4)
#define jump_table_size (1000)
#define MAX_VARIATIONS (50)

#define FLAME_ITTS (4)
#define RUN_FOREVER (1)

// remove this define to disable fractal visualization and free up another core for generation
#define DISPLAY enabled


#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#endif

static const __m256 zerovec   = { 0, 0, 0, 0, 0, 0, 0, 0 };
static const __m256 onevec    = { 1, 1, 1, 1, 1, 1, 1, 1 };
static const __m256 twovec    = { 2, 2, 2, 2, 2, 2, 2, 2 };
static const __m256 threevec  = { 3, 3, 3, 3, 3, 3, 3, 3 };
static const __m256 fourvec   = { 4, 4, 4, 4, 4, 4, 4, 4 };
static const __m256 pivec     = { PI, PI, PI, PI, PI, PI, PI, PI };
static const __m256 negonevec = { -1, -1, -1, -1, -1, -1, -1, -1 };
static const __m256 halfvec   = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
static const __m128 halfvec128   = { 0.5, 0.5, 0.5, 0.5 };

#endif