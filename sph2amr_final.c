/*
 * SPH2AMR.c
 *
 * Reads in Gadget2 data and projects the data onto a uniform grid suitable for
 * Orion2, which can be read in during initialization
 * Written by Athena Stacy and Aaron Lee, 2014
 *
 * Takes two inputs: the number of cells in each dimension (integer) and the width
 * of the box to carve out, in units of parsecs.
 *
 * To use, you need think about / set a few parameters before compiling:
 *
 *
 * preprocessor BREAD (line ~38ish): 0 or 1, whether you do not or do want to use
 *          the bread model for distributing data to the processors.
 *          0 : If no, then each processor allocates memory for the entire grid,
 *          and loops over a fraction of the particles. Uses a lot of memory for
 *          large grids.
 *          1 : Each processor gets a fraction of the grids and loops over all
 *          the particles.
 *
 * Filename to read in (line ~216ish): Basename of file, the number of files, etc.
 * Filename to print out (line ~351ish): Filename to print out...
 *
 */


/* PreProcessor Directives =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */

/* Debugging */
#define DEBUGGING 0

/* Diagnostic Info? (Automatically included if DEBUGGING) */
#define CHATTY 1


/* End PreProcessor Directives =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */


#include "sph2amr_final.h"



int main(int argc, char **argv)
{
    int i = 0;
    
    return i;
}




