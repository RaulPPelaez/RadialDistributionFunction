/*Raul P. Pelaez 2017. CPU Radial Distribution Function computer

  A class that computes rdf on the CPU, can process a snapshot and return the normalized rdf.


  Usage:

  Create with RadialDistributionFunctionCPU rdfComputerCPU;

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
#include<thrust/device_vector.h>
#include<limits>
#include"NeighbourListCPU.h"
namespace gdr{
  
  //A functor for the neighbour list that sums to the histogram the distance between two particles
  template<class vecType, bool fixBinBIAS>  
  struct DistanceCounter{
    using pairDistanceCountType = typename pairDistanceCounterType<fixBinBIAS>::type;
    const vecType* pos;
    pairDistanceCountType* pairDistanceCount;
    Box3D box;
    real rcut;
    real binSize;
    bool is3D;
    DistanceCounter(const vecType* pos,
		    pairDistanceCountType *pairDistanceCount,
		    const Configuration &config):
      pos(pos),
      pairDistanceCount(pairDistanceCount),
      box(config.boxSize),
      rcut(config.maxDistance),
      binSize(config.maxDistance/config.numberBins)
    {
	this->is3D = config.dimension == Configuration::dimensionality::D3;
	
    }
    inline void operator()(int index_i, int index_j){
      if(index_i > index_j){
	real3 rij = box.apply_pbc(make_real3(pos[index_i]) - make_real3(pos[index_j]));
	
	real r = sqrtf(dot(rij, rij));
	if(r<rcut){
	  int bin=floorf(r/binSize);

	  if(fixBinBIAS){
	    const real rbin = (bin+0.5)*binSize;
	    real norm;
	    if(is3D){
	      norm = rbin*rbin/(r*r);
	    }
	    else{
	      norm = rbin/r;
	    }
	    pairDistanceCount[bin]+=2.0*norm;
	  }
	  else
	    pairDistanceCount[bin]+=2;
	}
      }
    }

  };

  template<bool fixBinBIAS>
  class RadialDistributionFunctionCPU{    
    using pairDistanceCountType = typename pairDistanceCounterType<fixBinBIAS>::type;
    std::vector<pairDistanceCountType> pairDistanceCount;
    NeighbourListCPU neighbourList;
    int processedSnapshots = 0;
    std::vector<real2> rdf_mean_and_var; //Current mean and variance of the rdf
    std::vector<double> count2rdf; //Conversion factor between pairDistanceCount and radial distribution function
  public:
  
    RadialDistributionFunctionCPU(){  }

    //Compute the pair distance histogram (proportional to rdf), meant to be called once per snapshot
    //Each time it is called the histogram is summed to the previous one
    template<class vecType>
    void processSnapshot(const vecType *pos, const Configuration &config){
    if(!pos){std::cerr<<"ERROR: position pointer is NULL!! in gdr CPU"<<std::endl;return; }

    //Lazy initialization of the arrays
    if(pairDistanceCount.size() != config.numberBins) pairDistanceCount.resize(config.numberBins, 0);       

    DistanceCounter<vecType, fixBinBIAS> distanceCounter(pos, pairDistanceCount.data(), config);

    //Ask neighbourList for advice, if false, use Nbody
    if(neighbourList.shouldUse(config)){
      neighbourList.makeList(pos, config);
      neighbourList.transverseList(pos, distanceCounter, config);
    }
    else{
      int N = config.numberParticles;
      for(int i=0; i<N; i++){
	for(int j=i+1; j<N; j++){
	  distanceCounter(j,i);
	}
      }
    }

    //Compute conversion factor if necessary
    if(count2rdf.size()!=config.numberBins){
      count2rdf.resize(config.numberBins, 0);
      computeCount2rdf(config, count2rdf.data());
    }
    
    if(rdf_mean_and_var.size()!= config.numberBins) rdf_mean_and_var.resize(config.numberBins, real2());

    //Compute mean and variance
    int time = processedSnapshots;
    for(int i = 0; i<config.numberBins; i++){
      //rdf in this snapshot
      double rdf = pairDistanceCount[i]*count2rdf[i];
      
      double mean = rdf_mean_and_var[i].x;
      
      rdf_mean_and_var[i].x += (rdf - mean)/double(time + 1); //Update mean
      rdf_mean_and_var[i].y += time*pow(mean - rdf,2)/double(time+1); //Update variance
      
      //Reset count CPU
      pairDistanceCount[i] = 0;     
    }


    
    processedSnapshots++;
  }

    //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
    void getRadialDistributionFunction(real *rdf, real *std, const Configuration &config){
      // //pair distance count to radial function distribution
      // normalizeRadialDistributionFunction(rdf, std, pairDistanceCount.data() , config, this->processedSnapshots);

      int T = processedSnapshots;

      for(int i=0; i<config.numberBins; i++){
	rdf[i] = rdf_mean_and_var[i].x;
	if(T==1)
	  std[i] = std::numeric_limits<real>::quiet_NaN();
	else
	  std[i] = sqrt(rdf_mean_and_var[i].y)/sqrt(T*max(T-1,1));
      
      }

    }
    void reset(){
      this->processedSnapshots = 0;
    }
    
  };


  

  
}
