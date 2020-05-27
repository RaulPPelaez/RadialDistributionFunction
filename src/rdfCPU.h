/*Raul P. Pelaez 2017-2020. CPU Radial Distribution Function computer

  A class that computes rdf on the CPU, can process snapshots and return the normalized rdf.

  Usage:

  call
  rdfComputerCPU.processSnapshot for all snapshots just once

  call
  rdfComputerCPU.getRadialDistributionFunction anytime you want the last version of the rdf


 */
#include"vector_algebra.cuh"
#include"config.h"
#include"utils.cuh"
#include"rdf_common.h"
#include<vector>
#include<limits>
#include"NeighbourListCPU.h"
namespace gdr{

  template<bool fixBinBIAS>
  class RadialDistributionFunctionCPU{
    using pairDistanceCountType = typename pairDistanceCounterType<fixBinBIAS>::type;
    std::vector<pairDistanceCountType> pairDistanceCount;
    NeighbourListCPU neighbourList;
    int processedSnapshots = 0;
    std::vector<real2> rdf_mean_and_var;
    Configuration config;
    std::vector<int> numberParticlesPerType;
  public:

    RadialDistributionFunctionCPU(Configuration config, std::vector<int> numberParticlesPerType):
      config(config), numberParticlesPerType(numberParticlesPerType){
      int Ntypes = numberParticlesPerType.size();
      int numberRDFs = Ntypes*(Ntypes+1)/2;
      rdf_mean_and_var.resize(config.numberBins*numberRDFs, real2());
    }

    void processSnapshot(const real4 *pos){
      resetContainers();
      computeDistanceCount(pos);
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
    }

  private:

    void resetContainers(){
      int Ntypes = numberParticlesPerType.size();
      int numberRDFs = Ntypes*(Ntypes+1)/2;
      pairDistanceCount.resize(config.numberBins*numberRDFs);
      std::fill(pairDistanceCount.begin(), pairDistanceCount.end(), 0);
    }

    void computeDistanceCount(const real4* pos){
      if(not neighbourList.shouldUse(config)){
	this->computeWithNBody(pos);
      }
      else{
	this->computeWithNeighbourList(pos);
      }
    }

    void updateRDFFromDistanceCount(){
      int Ntypes = numberParticlesPerType.size();
      for(int typei = 0; typei < Ntypes; typei++){
	for(int typej = typei; typej < Ntypes; typej++){
	  int typeIndex = triangularIndex(typei, typej, Ntypes);
	  int firstElement = typeIndex*config.numberBins;
	  updateAverageRDF(&pairDistanceCount[firstElement], &rdf_mean_and_var[firstElement], typei, typej);
	}
      }
    }

    PairCounterTransverser<fixBinBIAS> createPairCounter(){
      int Ntypes = numberParticlesPerType.size();
      Box3D box(config.boxSize);
      real rcut = config.maxDistance;
      real binSize = rcut/double(config.numberBins);
      PairCounterTransverser<fixBinBIAS> pairCounter(pairDistanceCount.data(),
						     Ntypes, config.numberBins,
						     box,
						     rcut,
						     binSize,
						     config.dimension==Configuration::dimensionality::D3);

      return pairCounter;
    }

    void computeWithNBody(const real4* pos){
      auto pairCounter = createPairCounter();
      int N = config.numberParticles;
      for(int i=0; i<N; i++){
	for(int j=i+1; j<N; j++){
	  pairCounter(pos[j], pos[i]);
	}
      }
    }

    void computeWithNeighbourList(const real4* pos){
      neighbourList.makeList(pos, config);
      auto pairCounter = createPairCounter();
      neighbourList.transverseList(pairCounter, config);
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
  };
}
