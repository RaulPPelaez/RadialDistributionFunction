/*Raul P. Pelaez 2017. Radial Distribution Function common  (to any GPU/CPU implementation) functions
  
 */
#ifndef RDF_COMMON_H
#define RDF_COMMON_H

#include"vector_algebra.cuh"
#include"config.h"


namespace gdr{
  //Takes an histogram of pair distances and normlized to get the RDF
  void normalizeRadialDistributionFunction(real *rdf,
					   real *std, //Standard Deviation
					   ullint *pairDistanceCount,
					   const Configuration &config,
					   int numberProcessedSnapshots //Number of snapshots summed in pairDistanceCount
					   ){
    real3 L = config.boxSize;
    int N = config.numberParticles;
    real binSize = config.maxDistance/config.numberBins;
    double V = L.x*L.y;
    if(config.dimension==Configuration::dimensionality::D3) V *= L.z;
      
    double prefactor = M_PI;
    if(config.dimension==Configuration::dimensionality::D3) prefactor *= 2.0;
      
    constexpr double countedTwice = 2.0;
      
    double normalization = numberProcessedSnapshots*countedTwice*prefactor*binSize*N*N/V;


    
    double invNormalization = 1.0/normalization;
    for(int i=1; i<=config.numberBins; i++){
      double R = (i-0.5)*binSize;
      double invR =1.0/R;
      double count2rdf = invR*invNormalization;
      if(config.dimension == Configuration::dimensionality::D3) count2rdf *= invR;
      ullint count = pairDistanceCount[i-1];
      rdf[i-1] = count*count2rdf;
      
      std[i-1] = sqrt(count*(1.0-count/pow(N*numberProcessedSnapshots, 2)))*count2rdf;

    }
  }
}


#endif
