/*Raul P. Pelaez 2017. GPU Radial Distribution Function computer

  A class that computes rdf on the GPU, can process a snapshot and return the normalized rdf.


  Usage:

  Create with RadialDistributionFunctionGPU rdfComputerGPU;

  call
  rdfComputerGPU.processSnapshot for all snapshots just once

  call
  rdfComputerGPU.getRadialDistributionFunction anytime you want the last version of the rdf


 */

#include"vector_algebra.cuh"
#include"config.h"
#include"utils.cuh"
#include"rdfGPU_kernel.cuh"
#include"rdf_common.h"
#include"NeighbourListGPU.cuh"
#include<vector>
#include<thrust/device_vector.h>
#include<limits>
#include"atomic.cuh"
namespace gdr{
  template<bool fixBinBIAS>
  struct PairCounterTransverser{
    using pairCountType = typename pairDistanceCounterType<fixBinBIAS>::type;
    Box3D box;
    real rcut;
    real binSize;
    pairCountType *pairDistanceCounterGPUPtr;
    bool is3D;
    PairCounterTransverser(pairCountType *pairDistanceCounterGPUPtr,
			   Box3D box, real rcut, real binSize, bool is3D):
      box(box), rcut(rcut), binSize(binSize),
      pairDistanceCounterGPUPtr(pairDistanceCounterGPUPtr), is3D(is3D){}
    template<class vecType>
    __device__ void operator ()(vecType pi, vecType pj){
      const real3 rij = box.apply_pbc(make_real3(pi) - make_real3(pj));
      const real r = sqrtf(dot(rij, rij));
      if(r < rcut){
	const int bin = floorf(r/binSize);
	//This mode takes into account that the pair distance
	// has a different weight depending on where it lies within a bin
	if(fixBinBIAS){
	  const real rbin = (bin+0.5)*binSize;
	  real norm;
	  if(is3D){
	    norm = rbin*rbin/(r*r);
	  }
	  else{
	    norm = rbin/r;
	  }
	  atomicAdd((real*)&pairDistanceCounterGPUPtr[bin], real(2.0)*norm);
	}
	else{
	  atomicAdd((ullint*)&pairDistanceCounterGPUPtr[bin], 2);
	}
      }
    }

  };

  template<bool fixBinBIAS>
  class RadialDistributionFunctionGPU{
    using pairDistanceCountType = typename pairDistanceCounterType<fixBinBIAS>::type;
    thrust::device_vector<real4> posGPU;
    thrust::device_vector<pairDistanceCountType> pairDistanceCountGPU;
    thrust::host_vector<pairDistanceCountType> pairDistanceCountCPU;
    std::vector<real2> rdf_mean_and_var; //Current mean and variance of the rdf
    int processedSnapshots = 0;
    std::shared_ptr<CellList> nl;
    std::vector<double> count2rdf; //Conversion factor between pairDistanceCount and radial distribution function
  public:

    RadialDistributionFunctionGPU(){}

    //Compute the pair distance histogram (proportional to rdf), meant to be called once per snapshot
    //Each time it is called the histogram is summed to the previous one
    template<class vecType>
    void processSnapshot(const vecType *posCPU, const Configuration &config){
      if(!posCPU){std::cerr<<"ERROR: position pointer is NULL!! in gdr GPU"<<std::endl;return; }
      if(posGPU.size() != config.numberParticles) posGPU.resize(config.numberParticles);
      if(pairDistanceCountCPU.size() != config.numberBins) pairDistanceCountCPU.resize(config.numberBins, 0);
      pairDistanceCountGPU = pairDistanceCountCPU;
      auto posGPUPtr=thrust::raw_pointer_cast(posGPU.data());
      cudaMemcpy(posGPUPtr, posCPU, config.numberParticles*sizeof(vecType), cudaMemcpyHostToDevice);
      real rcut = config.maxDistance;
      real3 L = config.boxSize;
      int nx = int(L.x/rcut + 0.5);
      int ny = int(L.y/rcut + 0.5);
      int nz = int(L.z/rcut + 0.5);
      if(nx<3 || ny<3 || nz<3 || config.numberParticles < 1e4){
	this->computeWithNBody(posGPUPtr, config);
      }
      else{
	this->computeWithNeighbourList(posGPUPtr, config);
      }
      if(count2rdf.size()!=config.numberBins){
	count2rdf.resize(config.numberBins, 0);
	computeCount2rdf(config, count2rdf.data());
      }
      if(rdf_mean_and_var.size()!= config.numberBins) rdf_mean_and_var.resize(config.numberBins, real2());
      pairDistanceCountCPU = pairDistanceCountGPU;
      int time = processedSnapshots;
      for(int i = 0; i<config.numberBins; i++){
	double rdf = pairDistanceCountCPU[i]*count2rdf[i];
	double mean = rdf_mean_and_var[i].x;
	rdf_mean_and_var[i].x += (rdf - mean)/double(time + 1);
	rdf_mean_and_var[i].y += time*pow(mean - rdf,2)/double(time+1);
	pairDistanceCountCPU[i] = 0;
      }
      processedSnapshots++;
    }

    template<class vecType>
    void computeWithNBody(const vecType *posGPU, const Configuration &config){
      int N = config.numberParticles;
      Box3D box(config.boxSize);
      real rcut = config.maxDistance;
      int numberBins= config.numberBins;
      real binSize = rcut/numberBins;
      auto pairDistanceCountGPUPtr=thrust::raw_pointer_cast(pairDistanceCountGPU.data());
      PairCounterTransverser<fixBinBIAS> pairCounter(pairDistanceCountGPUPtr,
						     box,
						     rcut,
						     binSize,
						     config.dimension==Configuration::dimensionality::D3);
      int BLOCKSIZE=128;
      int Nthreads = BLOCKSIZE<N?BLOCKSIZE:N;
      int Nblocks  = (N+Nthreads-1)/Nthreads;
      int numTiles = (N + Nthreads-1)/Nthreads;
      size_t sharedMemorySize =  Nthreads*(sizeof(vecType));
      nBody_rdfKernel<<<Nblocks, Nthreads, sharedMemorySize>>>(posGPU,
							       numTiles,
							       N,
							       pairCounter);
    }

    template<class vecType>
    void computeWithNeighbourList(const vecType *posGPU, const Configuration &config){
      Box3D box(config.boxSize);
      real rcut = config.maxDistance;
      int numberBins= config.numberBins;
      real binSize = rcut/numberBins;
      auto pairDistanceCountGPUPtr = thrust::raw_pointer_cast(pairDistanceCountGPU.data());
      if(!nl){
	nl = std::make_shared<CellList>();
      }
      nl->updateNeighbourList(posGPU, config);
      PairCounterTransverser<fixBinBIAS> pairCounter(pairDistanceCountGPUPtr,
						     box,
						     rcut,
						     binSize,
						     config.dimension==Configuration::dimensionality::D3);
      nl->transverseList(pairCounter);
    }

    void getRadialDistributionFunction(real *rdfCPU, real *stdCPU, const Configuration &config){
      int T = processedSnapshots;
      for(int i=0; i<config.numberBins; i++){
	rdfCPU[i] = rdf_mean_and_var[i].x;
	if(T==1)
	  stdCPU[i] = std::numeric_limits<real>::quiet_NaN();
	else
	  stdCPU[i] = sqrt(rdf_mean_and_var[i].y)/sqrt(T*max(T-1,1));
      }
    }

    void reset(){
      this->processedSnapshots = 0;
      std::fill(rdf_mean_and_var.begin(), rdf_mean_and_var.end(), real2());
    }

  };
}

