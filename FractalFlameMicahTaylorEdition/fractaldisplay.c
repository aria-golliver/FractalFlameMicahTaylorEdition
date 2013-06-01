#include "fractaldisplay.h"
#include "GL\glfw.h"
#include "histogram.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static GLRGB8 *texture;
static histocell *tmphistogram;

static u32 running = 0;

#define fullscreen 1

void displayinit(void){
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        exit( EXIT_FAILURE );
    }

    // Open a window and create its OpenGL context
    if( !glfwOpenWindow( swid, shei, 0,0,0,0, 0,0, fullscreen ? GLFW_FULLSCREEN : GLFW_WINDOW_NO_RESIZE ) )
    {
        fprintf( stderr, "Failed to open GLFW window\n" );

        glfwTerminate();
        exit( EXIT_FAILURE );
    }

    displayreset();

    running = 1;

    glfwSetWindowTitle( "fractal" );
}

void displaydistroy(void){
    if(texture)
        free(texture);

    texture = NULL;
    
    if(tmphistogram)
        free(tmphistogram);

    tmphistogram = NULL;

    running = 0;

    glfwTerminate();

    exit(EXIT_SUCCESS);
}

#define MAX(a,b) (a > b ? a : b) 
#define MAX3(a,b,c) MAX(a, MAX(b, c))

void displayreset(void){
    if(texture)
        memset(texture, 0, swid * shei * sizeof(histocell));
    else
        texture = (GLRGB8 *) calloc(swid * shei, sizeof(histocell));
    
    if(tmphistogram)
        memset(tmphistogram, 0, swid * shei * sizeof(histocell));
    else
        tmphistogram = (histocell *) calloc(swid * shei, sizeof(histocell));
}

void updateDisplay(void){
    if(!running)
        return;

    if(glfwGetKey( GLFW_KEY_ESC ) == GLFW_PRESS ||
           !glfwGetWindowParam( GLFW_OPENED ))
           displaydistroy();

    // this holds the maximum number of times a cell in the histogram has been hit by the algorithm
    f32 amax = 1;

    for(u64 y = 0; y < shei; y++){
        for(u64 x = 0; x < swid; x++){
            const u64 cells = x + (y * swid);
            const u64 cellh = (x * ss) + (y * ss * ss * swid);
            tmphistogram[cells] = histoget(cellh);
            amax = MAX(tmphistogram[cells].a, amax);
        }
    }

    const f32 logamax = log(amax);

    for(u64 y = 0; y < shei; y++){
        for(u64 x = 0; x < swid; x++){
            const u64 cell = x + (y * swid);
            const f32 a = log(tmphistogram[cell].a) / logamax;

            const u8 r = (tmphistogram[cell].r / tmphistogram[cell].a) * 0xFF * a;
            const u8 g = (tmphistogram[cell].g / tmphistogram[cell].a) * 0xFF * a;
            const u8 b = (tmphistogram[cell].b / tmphistogram[cell].a) * 0xFF * a;
            
            texture[cell].c[0] = r;
            texture[cell].c[1] = g;
            texture[cell].c[2] = b;
        }
    }
    
    glDrawPixels(swid, shei, GL_RGB, GL_UNSIGNED_BYTE, texture);

    // Swap buffers
    glfwSwapBuffers();
}