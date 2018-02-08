/*Raul P. Pelaez 2017. Radial Distribution Function common  (to any GPU/CPU implementation) functions
  
 */
#ifndef RDF_COMMON_H
#define RDF_COMMON_H

#include"vector_algebra.cuh"
#include"config.h"


namespace gdr{
  template<bool fixBinBIAS> struct pairDistanceCounterType;
  template<> struct pairDistanceCounterType<true >{using type= real;};
  template<> struct pairDistanceCounterType<false>{using type= ullint;};
  
  //Computes the conversion (normalization) factor between a pair count and the rdf for each bin
  void computeCount2rdf(const Configuration &config, double * count2rdf){
    real3 L = config.boxSize;
    int N = config.numberParticles;
    real binSize = config.maxDistance/config.numberBins;
    double V = L.x*L.y;
    if(config.dimension==Configuration::dimensionality::D3) V *= L.z;
      
    double prefactor = M_PI;
    if(config.dimension==Configuration::dimensionality::D3) prefactor *= 2.0;
      
    constexpr double countedTwice = 2.0;
      
    double normalization = countedTwice*prefactor*binSize*N*N/V;

    
    double invNormalization = 1.0/normalization;
    for(int i=0; i<config.numberBins; i++){
      double R = (i+0.5)*binSize;
      double invR =1.0/R;
      double count2rdf_i = invR*invNormalization;
      if(config.dimension == Configuration::dimensionality::D3) count2rdf_i *= invR;
      count2rdf[i] = count2rdf_i;
    }

  }
}


#endif
