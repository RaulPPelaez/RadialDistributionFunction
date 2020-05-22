/*Raul P. Pelaez 2017. GPU Radial Distribution Function computer kernels

Contains an NBody pair transverser. Passes through all possible index pairs in pos O(N^2).
 Computes the distance between each pair and sums 1 to its corresponding bin in the distance histogram



 */
#ifndef RDFGPU_KERNEL_CUH
#define RDFGPU_KERNEL_CUH

#include"vector_algebra.cuh"
#include"utils.cuh"

namespace gdr{

  //Goes through all pair of particles in the pos array, O(N^2).
  //With each pair, computes distance and sums 1 to the bin corresponding to that distance.
  /*Reference: Fast N-Body Simulation with CUDA. Chapter 31 of GPU Gems 3*/
  template<class PairCounter>
  __global__ void nBody_rdfKernel(const real4* __restrict__ pos,
				  int numTiles,
				  uint N,
				  PairCounter pairCounter
				  ){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    bool active = true;
    if(id>=N) active = false;
    extern __shared__ char shMem[];
    real4 *shPos = (real4*) shMem;
    real4 pi;
    if(active) {
      pi = pos[id];
    }
    for(int tile = 0; tile<numTiles; tile++){
      const int i_load = tile*blockDim.x+threadIdx.x;
      if(i_load<N){
	shPos[threadIdx.x] = pos[i_load];
      }
      __syncthreads();
#pragma unroll 8
      for(uint counter = 0; counter<blockDim.x; counter++){
	if(!active) break;
	int cur_j = tile*blockDim.x+counter;
	if(cur_j<N && cur_j>id && cur_j != id){
	  pairCounter(pi, shPos[counter]);
	}
      }
      __syncthreads();
    }
  }
}

#endif
