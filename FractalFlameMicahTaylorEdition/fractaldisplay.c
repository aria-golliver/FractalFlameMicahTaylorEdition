#include "fractaldisplay.h"
#include "GL\glfw.h"
#include "histogram.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>



typedef struct GLRGB8_t{
    char r, g, b;
} GLRGB8;

static GLRGB8 texture[hwid * hhei];

void displayinit(void){
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        exit( EXIT_FAILURE );
    }

    // Open a window and create its OpenGL context
    if( !glfwOpenWindow( hwid, hhei, 0,0,0,0, 0,0, GLFW_WINDOW ) )
    {
        fprintf( stderr, "Failed to open GLFW window\n" );

        glfwTerminate();
        exit( EXIT_FAILURE );
    }

    glfwSetWindowTitle( "Spinning Triangle" );
}
void updateDisplay(void){
    memset(texture, 0, sizeof(texture));

    int width, height, x, y;
    double t;
    t = glfwGetTime();
    glfwGetMousePos( &x, &y );

    // Get window size (may be different than the requested size)
    glfwGetWindowSize( &width, &height );

    // Special case: avoid division by zero below
    height = height > 0 ? height : 1;

    glViewport( 0, 0, width, height );

    // Clear color buffer to black
    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT );

    // Select and setup the projection matrix
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 65.0f, (GLfloat)width/(GLfloat)height, 1.0f, 100.0f );

    // Select and setup the modelview matrix
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt( 0.0f, 1.0f, 0.0f,    // Eye-position
               0.0f, 20.0f, 0.0f,   // View-point
               0.0f, 0.0f, 1.0f );  // Up-vector

    // Draw a rotating colorful triangle
    glTranslatef( 0.0f, 14.0f, 0.0f );
    glRotatef( 0.3f*(GLfloat)x + (GLfloat)t*100.0f, 0.0f, 0.0f, 1.0f );
    glBegin( GL_TRIANGLES );
      glColor3f( 1.0f, 0.0f, 0.0f );
      glVertex3f( -5.0f, 0.0f, -4.0f );
      glColor3f( 0.0f, 1.0f, 0.0f );
      glVertex3f( 5.0f, 0.0f, -4.0f );
      glColor3f( 0.0f, 0.0f, 1.0f );
      glVertex3f( 0.0f, 0.0f, 6.0f );
    glEnd();

    // Swap buffers
    glfwSwapBuffers();
}