/////////////////////////////////////////////////////////////////////////////

// The Software is provided "AS IS" and possibly with faults. 

// Intel disclaims any and all warranties and guarantees, express, implied or

// otherwise, arising, with respect to the software delivered hereunder,

// including but not limited to the warranty of merchantability, the warranty

// of fitness for a particular purpose, and any warranty of non-infringement

// of the intellectual property rights of any third party.

// Intel neither assumes nor authorizes any person to assume for it any other

// liability. Customer will use the software at its own risk. Intel will not

// be liable to customer for any direct or indirect damages incurred in using

// the software. In no event will Intel be liable for loss of profits, loss of

// use, loss of data, business interruption, nor for punitive, incidental,

// consequential, or special damages of any kind, even if advised of

// the possibility of such damages.

//

// Copyright (c) 2003 Intel Corporation

//

// Third-party brands and names are the property of their respective owners

//

///////////////////////////////////////////////////////////////////////////

// Random Number Generation for SSE / SSE2

// Source File

// Version 0.1

// Author Kipp Owens, Rajiv Parikh

////////////////////////////////////////////////////////////////////////



#ifndef RAND_SSE_H

#define RAND_SSE_H RAND_SSE_H

void srand_sse( unsigned int seed , unsigned int th_id);

void rand_sse( unsigned int* , unsigned int th_id);

#endif