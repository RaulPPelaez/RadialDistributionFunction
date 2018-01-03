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
#include"NeighbourListGPU.cuh"
#include<vector>
#include<thrust/device_vector.h>

namespace gdr{

  class RadialDistributionFunctionGPU{
    thrust::device_vector<real4> posGPU;
    thrust::device_vector<ullint> pairDistanceCountGPU;
    int processedSnapshots = 0;
    shared_ptr<CellList> nl;
  public:
  
    RadialDistributionFunctionGPU(){}

    //Compute the pair distance histogram (proportional to rdf), meant to be called once per snapshot
    //Each time it is called the histogram is summed to the previous one
    template<class vecType>
    void processSnapshot(const vecType *posCPU, const Configuration &config);
    template<class vecType>
    void computeWithNBody(const vecType *posGPU, const Configuration &config);
    template<class vecType>
    void computeWithNeighbourList(const vecType *posGPU, const Configuration &config);
    
    //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
    void getRadialDistributionFunction(real *rdfCPU, real *stdCPU, const Configuration &config);
    void reset(){
      this->processedSnapshots = 0;
    }
    
  };


  //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
  void RadialDistributionFunctionGPU::getRadialDistributionFunction(real *rdfCPU, real *stdCPU, const Configuration &config){
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

    if(posGPU.size() != config.numberParticles) posGPU.resize(config.numberParticles);
    if(pairDistanceCountGPU.size() != config.numberBins) pairDistanceCountGPU.resize(config.numberBins, 0);

    //Get raw pointers to device memory
    auto posGPUPtr=thrust::raw_pointer_cast(posGPU.data());      
    cudaMemcpy(posGPUPtr, posCPU, config.numberParticles*sizeof(vecType), cudaMemcpyHostToDevice);

    real rcut = config.maxDistance;
    real3 L = config.boxSize;
    int nx = int(L.x/rcut + 0.5);
    int ny = int(L.y/rcut + 0.5);
    int nz = int(L.z/rcut + 0.5);
    if(nx<3 && ny<3 && nz<3){
      this->computeWithNBody(posGPUPtr, config);
    }
    else{
      this->computeWithNeighbourList(posGPUPtr, config);
    }

    
    processedSnapshots++;
  }

  template<class vecType>
  void RadialDistributionFunctionGPU::computeWithNBody(const vecType *posGPU, const Configuration &config){
    //Recover parameters
    int N = config.numberParticles;
    Box3D box(config.boxSize);
    real rcut = config.maxDistance;
    int numberBins= config.numberBins;
    real binSize = rcut/numberBins;

    //Get raw pointers to device memory
    auto pairDistanceCountGPUPtr=thrust::raw_pointer_cast(pairDistanceCountGPU.data());
      
    //Configure and lauch kernel
    int BLOCKSIZE=128;
    int Nthreads = BLOCKSIZE<N?BLOCKSIZE:N;
    int Nblocks  = (N+Nthreads-1)/Nthreads;
    int numTiles = (N + Nthreads-1)/Nthreads;        
  
    size_t sharedMemorySize =  Nthreads*(sizeof(vecType));
      
    nBody_rdfKernel<<<Nblocks, Nthreads, sharedMemorySize>>>(posGPU,
							     numTiles,
							     N,
							     box,
							     rcut,
							     binSize,
							     pairDistanceCountGPUPtr);


  }

  struct PairCounterTransverser{
    Box3D box;
    real rcut;
    real binSize;
    ullint *pairDistanceCounterGPUPtr;
    PairCounterTransverser(ullint *pairDistanceCounterGPUPtr,
			   Box3D box, real rcut, real binSize):
      box(box), rcut(rcut), binSize(binSize),
      pairDistanceCounterGPUPtr(pairDistanceCounterGPUPtr){}
    __device__ void operator ()(real4 pi, real4 pj){
      real3 rij = box.apply_pbc(make_real3(pi) - make_real3(pj));
      real r = sqrtf(dot(rij, rij));
      if(r < rcut){
	int bin = floorf(r/binSize);
	atomicAdd(&pairDistanceCounterGPUPtr[bin], 2);
      }
    }

  };
  template<class vecType>
  void RadialDistributionFunctionGPU::computeWithNeighbourList(const vecType *posGPU, const Configuration &config){
    //Recover parameters
    Box3D box(config.boxSize);
    real rcut = config.maxDistance;
    int numberBins= config.numberBins;
    real binSize = rcut/numberBins;

    //Get raw pointers to device memory    
    auto pairDistanceCountGPUPtr = thrust::raw_pointer_cast(pairDistanceCountGPU.data());

    
    if(!nl){
      nl = make_shared<CellList>();
    }
    nl->updateNeighbourList(posGPU, config);

    PairCounterTransverser pairCounter(pairDistanceCountGPUPtr,
				       box,
				       rcut,
				       binSize);
				       
				       
    nl->transverseList(pairCounter);

  }

}