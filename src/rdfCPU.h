/*Raul P. Pelaez 2017. CPU Radial Distribution Function computer

  A class that computes rdf on the CPU, can process a snapshot and return the normalized rdf.


  Usage:

  Create with RadialDistributionFunctionCPU rdfComputerCPU;

  call 
  rdfComputerCPU.processSnapshot for all snapshots just once

  call 
  rdfComputerCPU.getRadialDistributionFunction anytime you want the last version of the rdf


 */
#ifndef RDFCPU_CUH
#define RDFCPU_CUH
#include"vector_algebra.cuh"
#include"config.h"
#include"utils.cuh"
#include"rdf_common.h"

#include<vector>
#include<thrust/device_vector.h>

namespace gdr{

  class RadialDistributionFunctionCPU{    
    std::vector<ullint> pairDistanceCount;
    int processedSnapshots = 0;
  public:
  
    RadialDistributionFunctionCPU(){

      cerr<<"CPU version not implemented yet!!"<<endl;
      exit(1);
    }

    //Compute the pair distance histogram (proportional to rdf), meant to be called once per snapshot
    //Each time it is called the histogram is summed to the previous one
    template<class vecType>
    void processSnapshot(const vecType *pos, const Configuration &config);

    //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
    void getRadialDistributionFunction(real *rdf, real *stdCPU, const Configuration &config);
    void reset(){
      this->processedSnapshots = 0;
    }
    
  };


  //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
  void RadialDistributionFunctionCPU::getRadialDistributionFunction(real *rdf, real *std, const Configuration &config){           
    //pair distance count to radial function distribution
    normalizeRadialDistributionFunction(rdf, std, pairDistanceCount.data() , config, this->processedSnapshots);
  }

  

  
}
#endif


namespace gdr{

  //Compute the pair distance histogram, meant to be called once per snapshot
  //Each time it is called the histogram is summed to the previous one
  template<class vecType>
  void RadialDistributionFunctionCPU::processSnapshot(const vecType *pos, const Configuration &config){
    if(!pos){cerr<<"ERROR: position pointer is NULL!! in gdr CPU"<<endl;return; }

    //Recover parameters
    int N = config.numberParticles;
    Box3D box(config.boxSize);
    real rcut = config.maxDistance;
    int numberBins= config.numberBins;
    real binSize = rcut/numberBins;
    
    if(pairDistanceCount.size() != N) pairDistanceCount.resize(N, 0);

    
    
    
    processedSnapshots++;
  }
}
