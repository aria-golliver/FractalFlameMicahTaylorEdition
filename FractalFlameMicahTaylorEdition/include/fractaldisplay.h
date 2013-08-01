#ifndef __FRACTALDISPLAY_H__
#define __FRACTALDISPLAY_H__ __FRACTALDISPLAY_H__

#define fullscreen 1

// initializes the opengl display window and temporary buffers
void displayinit(void);

// distroys the opengl display and temporary buffers
void displaydistroy(void);

// updates fractal display
void updateDisplay(void);

// sets all buffers to 0
void displayreset(void);

#endif