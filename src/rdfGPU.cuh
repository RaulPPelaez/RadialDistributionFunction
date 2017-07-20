/*Raul P. Pelaez 2017. GPU Radial Distribution Function computer

  A class that computes rdf on the GPU, can process a snapshot and return the normalized rdf.


  Usage:

  Create with RadialDistributionFunctionGPU rdfComputerGPU;

  call 
  rdfComputerGPU.processSnapshot for all snapshots just once

  call 
  rdfComputerGPU.getRadialDistributionFunction anytime you want the last version of the rdf


 */
#ifndef RDFGPU_CUH
#define RDFGPU_CUH
#include"vector_algebra.cuh"
#include"config.h"
#include"utils.cuh"
#include"rdfGPU_kernel.cuh"
#include"rdf_common.h"

#include<vector>
#include<thrust/device_vector.h>

namespace gdr{

  class RadialDistributionFunctionGPU{
    thrust::device_vector<real4> posGPU;
    thrust::device_vector<ullint> pairDistanceCountGPU;
    int processedSnapshots = 0;
  public:
  
    RadialDistributionFunctionGPU(){}

    //Compute the pair distance histogram (proportional to rdf), meant to be called once per snapshot
    //Each time it is called the histogram is summed to the previous one
    template<class vecType>
    void processSnapshot(const vecType *posCPU, const Configuration &config);

    //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
    void getRadialDistributionFunction(real *rdfCPU, real *stdCPU, const Configuration &config);
    void reset(){
      this->processedSnapshots = 0;
    }
    
  };


  //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
  void RadialDistributionFunctionGPU::getRadialDistributionFunction(real *rdfCPU, real *stdCPU, const Configuration &config){
    auto pairDistanceCountGPUPtr=thrust::raw_pointer_cast(pairDistanceCountGPU.data());
    //Downloads pair histogram
    thrust::host_vector<ullint> countCPU = pairDistanceCountGPU;
   
    ullint *countCPUPtr=thrust::raw_pointer_cast(countCPU.data());
    //pair distance count to radial function distribution
    normalizeRadialDistributionFunction(rdfCPU, stdCPU, countCPUPtr, config, this->processedSnapshots);
  }

  

  
}
#endif


namespace gdr{

  //Compute the pair distance histogram, meant to be called once per snapshot
  //Each time it is called the histogram is summed to the previous one
  template<class vecType>
  void RadialDistributionFunctionGPU::processSnapshot(const vecType *posCPU, const Configuration &config){
    if(!posCPU){cerr<<"ERROR: position pointer is NULL!! in gdr GPU"<<endl;return; }

    //Recover parameters
    int N = config.numberParticles;
    Box3D box(config.boxSize);
    real rcut = config.maxDistance;
    int numberBins= config.numberBins;
    real binSize = rcut/numberBins;

    if(posGPU.size() != N) posGPU.resize(N);
    if(pairDistanceCountGPU.size() != N) pairDistanceCountGPU.resize(N, 0);

    //Get raw pointers to device memory
    auto posGPUPtr=thrust::raw_pointer_cast(posGPU.data());
    auto pairDistanceCountGPUPtr=thrust::raw_pointer_cast(pairDistanceCountGPU.data());
      
    cudaMemcpy(posGPUPtr, posCPU, N*sizeof(vecType), cudaMemcpyHostToDevice);

    //Configure and lauch kernel
    int BLOCKSIZE=128;
    int Nthreads = BLOCKSIZE<N?BLOCKSIZE:N;
    int Nblocks  = (N+Nthreads-1)/Nthreads;
    int numTiles = (N + Nthreads-1)/Nthreads;        
  
    size_t sharedMemorySize =  Nthreads*(sizeof(real3));
      
    nBody_rdfKernel<<<Nblocks, Nthreads, sharedMemorySize>>>(posGPUPtr,
							     numTiles,
							     N,
							     box,
							     rcut,
							     binSize,
							     pairDistanceCountGPUPtr);
    processedSnapshots++;
  }
}