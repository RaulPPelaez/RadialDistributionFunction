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
  template<typename BoxType, typename vectorLoadType>
  __global__ void nBody_rdfKernel(const vectorLoadType* __restrict__ pos,
				  int numTiles, /*Thread paralellism level, 
						  controls how many elements are stored in 
						  shared memory and
						  computed in parallel between synchronizations*/
				  uint N,       //Number of particles
				  BoxType box,  //Box object for PBC, can be Box3D of Box2D
				  real rcut,    //Maximum distance
				  real binSize,
				  ullint* __restrict__ pairDistanceCount
				  ){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    /*All threads must pass through __syncthreads, 
      but when N is not a multiple of 32 some threads are assigned a particle i>N.
      This threads cant return, so they are masked to not do any work*/
    bool active = true;
    if(id>=N) active = false;
    
    /*Each thread handles the interaction between particle id and all the others*/
    /*Storing blockDim.x positions in shared memory and processing all of them in parallel*/
    extern __shared__ char shMem[];

    real3 *shPos = (real3*) shMem;    
  
    real3 pi;
    if(active) {
      pi = make_real3(pos[id]); /*My position*/
    }    
    /*Distribute the N particles in numTiles tiles.
      Storing in each tile blockDim.x positions in shared memory*/
    /*This way all threads are accesing the same memory addresses at the same time*/
    for(int tile = 0; tile<numTiles; tile++){
      /*Load this tiles particles positions to shared memory*/
      const int i_load = tile*blockDim.x+threadIdx.x;
      if(i_load<N){ /*Even if im not active,
		      my thread may load a position each tile to shared memory.*/	
	shPos[threadIdx.x] = make_real3(pos[i_load]);
      }
      /*Wait for all threads to arrive*/
      __syncthreads();
      /*Go through all the particles in the current tile*/
#pragma unroll 8
      for(uint counter = 0; counter<blockDim.x; counter++){
	if(!active) break; /*An out of bounds thread must be masked*/
	int cur_j = tile*blockDim.x+counter; 
	if(cur_j<N && cur_j>id && cur_j != id){/*If the current particle exists, compute and accumulate*/
	  /*Compute and accumulate the current pair*/
	  real3 rij = pi - shPos[counter];
	  box.apply_pbc(rij);
	  real r = sqrtf(dot(rij, rij));
	  if(r<rcut){
	    int bin = floorf(r/binSize);
	    atomicAdd(&pairDistanceCount[bin], 2);
	  }
	
	}
      }/*End of particles in tile loop*/
      __syncthreads();

    }/*End of tile loop*/
  }
}

#endif
