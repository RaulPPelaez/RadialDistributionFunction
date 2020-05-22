/*Raul P. Pelaez 2017-2020. GPU Radial Distribution Function computer

  A class that computes rdf on the GPU, can process snapshots and return the normalized rdf.


  Usage:

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
#include<algorithm>
namespace gdr{
  template<bool fixBinBIAS>
  class RadialDistributionFunctionGPU{
    using pairDistanceCountType = typename pairDistanceCounterType<fixBinBIAS>::type;
    thrust::device_vector<real4> posGPU;
    thrust::device_vector<pairDistanceCountType> pairDistanceCountGPU;
    thrust::host_vector<pairDistanceCountType> pairDistanceCountCPU;
    std::vector<real2> rdf_mean_and_var; //Current mean and variance of the rdf
    int processedSnapshots = 0;
    std::shared_ptr<CellList> nl;
    Configuration config;
    std::vector<int> numberParticlesPerType;
  public:

    RadialDistributionFunctionGPU(Configuration config, std::vector<int>& numberParticlesPerType):
      config(config), numberParticlesPerType(numberParticlesPerType){
      int Ntypes = numberParticlesPerType.size();
      int numberRDFs = Ntypes*(Ntypes+1)/2;
      rdf_mean_and_var.resize(config.numberBins*numberRDFs, real2());
    }

    void processSnapshot(const real4 *posCPU){
      resetContainers();
      storePositions(posCPU);
      computeDistanceCount();
      updateRDFFromDistanceCount();
      processedSnapshots++;
    }

    std::vector<real2> getRadialDistributionFunction(){
      int T = processedSnapshots;
      std::vector<real2> rdf(rdf_mean_and_var.size());
      std::transform(rdf_mean_and_var.begin(), rdf_mean_and_var.end(),
		     rdf.begin(),
		     [T](real2 avgRDF){
		       real2 rdf;
		       rdf.x = avgRDF.x;
		       if(T == 0){
			 rdf.y = std::numeric_limits<real>::quiet_NaN();
		       }
		       else{
			 rdf.y = sqrt(avgRDF.y)/sqrt(T*std::max(T-1,1));
		       }
		       return rdf;
		     });
      return std::move(rdf);
    }

    void reset(){
      this->processedSnapshots = 0;
      std::fill(rdf_mean_and_var.begin(), rdf_mean_and_var.end(), real2());
    }

  private:

    void resetContainers(){
      int Ntypes = numberParticlesPerType.size();
      int numberRDFs = Ntypes*(Ntypes+1)/2;
      posGPU.resize(config.numberParticles);
      pairDistanceCountCPU.resize(config.numberBins*numberRDFs);
      std::fill(pairDistanceCountCPU.begin(), pairDistanceCountCPU.end(), 0);
      pairDistanceCountGPU = pairDistanceCountCPU;
    }

    void storePositions(const real4 *posCPU){
      thrust::copy(posCPU, posCPU + config.numberParticles, posGPU.begin());
    }

    bool shouldUseNeighbourList(){
      real rcut = config.maxDistance;
      real3 L = config.boxSize;
      int nx = int(L.x/rcut + 0.5);
      int ny = int(L.y/rcut + 0.5);
      int nz = int(L.z/rcut + 0.5);
      if(nx<3 or ny<3 or nz<3 or config.numberParticles < 1e4){
	return false;
      }
      return true;
    }

    void computeDistanceCount(){
      bool useNeighbourList = shouldUseNeighbourList();
      if(not useNeighbourList){
	this->computeWithNBody();
      }
      else{
	this->computeWithNeighbourList();
      }
    }

    void updateAverageRDF(pairDistanceCountType* pairCount, real2 *averageRDF, int typei, int typej){
      const int time = processedSnapshots;
      const double prefactor = computeCount2rdf(config);
      const int Ni = numberParticlesPerType[typei];
      const int Nj = numberParticlesPerType[typej];
      for(int i = 0; i < config.numberBins; i++){
	double R = binToDistance(i, config);
	double normalization = 2.0*prefactor/(Ni*Nj*R);
	if(config.dimension == Configuration::dimensionality::D3){
	  normalization /= R;
	}
	double rdf = pairCount[i]*normalization;
	double mean = averageRDF[i].x;
	averageRDF[i].x += (rdf - mean)/double(time + 1);
	averageRDF[i].y += time*pow(mean - rdf,2)/double(time+1);
      }
    }

    void updateRDFFromDistanceCount(){
      int Ntypes = numberParticlesPerType.size();
      pairDistanceCountCPU = pairDistanceCountGPU;
      for(int typei = 0; typei < Ntypes; typei++){
	for(int typej = typei; typej < Ntypes; typej++){
	  int typeIndex = triangularIndex(typei, typej, Ntypes);
	  int firstElement = typeIndex*config.numberBins;
	  updateAverageRDF(&pairDistanceCountCPU[firstElement], &rdf_mean_and_var[firstElement], typei, typej);
	}
      }
    }

    PairCounterTransverser<fixBinBIAS> createPairCounter(){
      int Ntypes = numberParticlesPerType.size();
      Box3D box(config.boxSize);
      real rcut = config.maxDistance;
      int numberBins= config.numberBins;
      real binSize = rcut/numberBins;
      auto pairDistanceCountGPUPtr = thrust::raw_pointer_cast(pairDistanceCountGPU.data());
      PairCounterTransverser<fixBinBIAS> pairCounter(pairDistanceCountGPUPtr,
						     Ntypes, config.numberBins,
						     box, rcut, binSize,
						     config.dimension==Configuration::dimensionality::D3);
      return pairCounter;
    }

    void computeWithNBody(){
      int N = config.numberParticles;
      auto pairCounter = createPairCounter();
      int BLOCKSIZE=128;
      int Nthreads = BLOCKSIZE<N?BLOCKSIZE:N;
      int Nblocks  = (N+Nthreads-1)/Nthreads;
      int numTiles = (N + Nthreads-1)/Nthreads;
      size_t sharedMemorySize =  Nthreads*(sizeof(real4));
      auto posGPU_ptr = thrust::raw_pointer_cast(posGPU.data());
      nBody_rdfKernel<<<Nblocks, Nthreads, sharedMemorySize>>>(posGPU_ptr, numTiles, N, pairCounter);
    }

    void computeWithNeighbourList(){
      if(!nl){
	nl = std::make_shared<CellList>();
      }
      auto posGPU_ptr = thrust::raw_pointer_cast(posGPU.data());
      nl->updateNeighbourList(posGPU_ptr, config);
      auto pairCounter = createPairCounter();
      nl->transverseList(pairCounter);
    }

  };
}

