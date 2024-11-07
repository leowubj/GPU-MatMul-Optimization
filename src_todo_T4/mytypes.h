#ifndef MYTYPES_H
#define MYTYPES_H

// tilescale (# of points computed by each thread)
#ifndef TILESCALE_M
#define TILESCALE_M 1 // Enter your own values
#endif
#ifndef TILESCALE_N
#define TILESCALE_N 1 // Enter your own values
#endif
#ifndef TILESCALE_K
#define TILESCALE_K 1 // Enter your own values
#endif

#define TILEDIM_M 128 // number of blocks in the x direction
#define TILEDIM_N 64 // number of blocks in the y direction

// matrix A loads
// with warps along the horiziontal axis (K)
// so to get good coalescaed loads, we want TILEDIM_K to be >= 32
//
#define TILEDIM_K 1 // Enter your own values

#define BLOCK_DOM_X 32 // Enter your own values
#define BLOCK_DOM_Y 8 // Enter your own values
#define NUM_WORK_Y (TILEDIM_M/BLOCK_DOM_Y)
#define NUM_WORK_X (TILEDIM_N/BLOCK_DOM_X)
#define TOTAL_WORK (NUM_WORK_Y*NUM_WORK_X)

// step size in each dimension
#define TILESTEP_N 1 // Enter your own values
#define TILESTEP_K 1 // Enter your own values
#define TILESTEP_M 1 // Enter your own values
 
#endif
