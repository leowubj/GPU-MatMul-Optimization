// ;-*- mode: c;-*-
// Matrix multiply device code
#include <assert.h>
#include <math.h>
#include "../src/utils.h"
#include "../src/types.h"
#include "mytypes.h"
using namespace std;

#include <stdio.h>

#define a(i, j, ld) a[ (i)*(ld) + (j) ]
#define b(i, j, ld) b[ (i)*(ld) + (j) ]
#define c(i, j, ld) c[ (i)*(ld) + (j) ]

#ifdef NAIVE
__global__ void matMul(int N, _FTYPE_ *C, _FTYPE_ *A, _FTYPE_ *B) {

    int I =  blockIdx.y*blockDim.y + threadIdx.y;
    int J =  blockIdx.x*blockDim.x + threadIdx.x;

    if((I < N) && (J < N)){
        _FTYPE_ _c = 0;
        for (unsigned int k = 0; k < N; k++) {
            _FTYPE_ a = A[I * N + k];
            _FTYPE_ b = B[k * N + J];
            _c += a * b;
        }
        C[I * N + J] = _c;
    }
}

#else
//You should be changing the kernel here for the non naive implementation.

//Shared memory
extern __shared__ _FTYPE_ smem[];

// matMul with only shared memory
/* __global__ void matMul(int N, _FTYPE_ *C, _FTYPE_ *A, _FTYPE_ *B) {
  int ty = threadIdx.y;
  int tx = threadIdx.x;
  int by = blockIdx.y;
  int bx = blockIdx.x;

  int abs_row_idx =  by*blockDim.y + ty;
  int abs_col_idx =  bx*blockDim.x + tx;

  _FTYPE_ * __restrict__ As = &smem[0];
  _FTYPE_ * __restrict__ Bs = &smem[TILEDIM_M * TILEDIM_M];
  _FTYPE_ Cij = 0;

  
  #pragma unroll
  for (int kk = 0; kk < gridDim.x; kk++){
    if (abs_row_idx*N + kk*TILEDIM_M + tx < N*N){
      As[ty*TILEDIM_M + tx] = A[abs_row_idx*N + kk*TILEDIM_M + tx];
    }
    else{
      As[ty*TILEDIM_M + tx] = 0;
    }
    if ((kk*TILEDIM_M+ty)*N + abs_col_idx < N*N){  
      Bs[ty*TILEDIM_M + tx] = B[(kk*TILEDIM_M+ty)*N + abs_col_idx];
    }
    else{
      Bs[ty*TILEDIM_M + tx] = 0;
    }
    __syncthreads();

    for (int k=0; k<TILEDIM_M; k++)
      Cij += As[ty*TILEDIM_M + k] * Bs[k * TILEDIM_M + tx];

    __syncthreads();
  }
  if((abs_row_idx < N) && (abs_col_idx < N)){
    C[abs_row_idx*N + abs_col_idx] = Cij;
  }
} */

__global__ void matMul(int N, _FTYPE_ *C, _FTYPE_ *A, _FTYPE_ *B) {
int ty = threadIdx.y;
  int tx = threadIdx.x;
  int by = blockIdx.y;
  int bx = blockIdx.x;

  int abs_row_idx =  by*blockDim.y + ty;
  int abs_col_idx =  bx*blockDim.x + tx;

  _FTYPE_ * __restrict__ As = &smem[0];
  _FTYPE_ * __restrict__ Bs = &smem[TILEDIM_M * TILEDIM_N];
  register _FTYPE_ Cij[TOTAL_WORK] = {0};

  if((abs_row_idx < N) && (abs_col_idx < N)){
  #pragma unroll
  for (int a = by*TILEDIM_M*N, b = bx*TILEDIM_N; a < N + by*TILEDIM_M*N; a += TILEDIM_N, b += TILEDIM_N*N) {
    int i;
    if (N % TILEDIM_N == 0){
      #pragma unroll
      for (i = 0; i < NUM_WORK_Y ; ++i) {
        #pragma unroll
        for(int j = 0; j < NUM_WORK_X; j++){
          As[(ty + i*BLOCK_DOM_Y)*TILEDIM_N + tx + j * BLOCK_DOM_X] = A[a + (ty + i*BLOCK_DOM_Y)*N + tx + j * BLOCK_DOM_X];
          Bs[(ty + i*BLOCK_DOM_Y)*TILEDIM_N + tx + j * BLOCK_DOM_X] = B[b + (ty + i*BLOCK_DOM_Y)*N + tx + j * BLOCK_DOM_X]; 
        }
      }
    }
    else{
      for (int i = 0; i < NUM_WORK_Y; ++i) {
        for(int j = 0; j < NUM_WORK_X; j++){
          int ity = ty + i * BLOCK_DOM_Y;

          int itx = tx + j * BLOCK_DOM_X;
          if ((a + ity * N >= N * N) || (a + itx >= (by * TILEDIM_M + 1) * N)) {
              As[ity * TILEDIM_N + itx] = 0;
          } else {
              As[ity * TILEDIM_N + itx] = A[a + ity * N + itx];
          }

          if ((b + ity * N >= N * N) || (bx * TILEDIM_N + itx >= N)) {
              Bs[ity * TILEDIM_N + itx] = 0;
          } else {
              Bs[ity * TILEDIM_N + itx] = B[b + ity * N + itx];
          }
        }
    }


    }

    __syncthreads(); 
    for (int k = 0; k < TILEDIM_N; k++) {
      #pragma unroll
      for (i = 0; i < NUM_WORK_Y; ++i) {
        #pragma unroll
        for(int j = 0; j < NUM_WORK_X; j++){
          Cij[i + j * NUM_WORK_Y] += As[(ty + i*BLOCK_DOM_Y)*TILEDIM_N + k] * Bs[k*TILEDIM_N + tx + j * BLOCK_DOM_X];
        }
      }
    }
    __syncthreads(); 
  }
    int c = N*TILEDIM_M*by + TILEDIM_N*bx;
    int row_idx = TILEDIM_M*by+ty;
    int col_idx = TILEDIM_N*bx+tx;
    
    if (N % TILEDIM_N != 0)
    {
      #pragma unroll
      for(int j = 0; j < NUM_WORK_X; j++){
        if (col_idx + j * BLOCK_DOM_X<N){
          #pragma unroll
          for (int i = 0; i < NUM_WORK_Y; ++i) {
            if (row_idx + i*BLOCK_DOM_Y < N){
              C[c + N * (ty + i*BLOCK_DOM_Y) + tx + j * BLOCK_DOM_X] = Cij[j*NUM_WORK_Y+i];
            }
          }
        }
      }
    }
    else{
      #pragma unroll
      for (int i = 0; i < NUM_WORK_Y; ++i) {
        for(int j = 0; j < NUM_WORK_X; j++){
          C[c + N * (ty + i*BLOCK_DOM_Y) + tx + j * BLOCK_DOM_X] = Cij[j * NUM_WORK_Y + i];
        }
      }
    }
  }
    
}

  



#endif
