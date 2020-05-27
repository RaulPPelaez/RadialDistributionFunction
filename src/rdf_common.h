/*Raul P. Pelaez 2017-2020. Radial Distribution Function common  (to any GPU/CPU implementation) functions

 */
#ifndef RDF_COMMON_H
#define RDF_COMMON_H

#include"vector_algebra.cuh"
#include"config.h"
#include"atomic.cuh"
namespace gdr{
  #ifdef USE_CUDA
#define HOSTDEVICE __host__ __device__
#else
#define HOSTDEVICE
#endif

  //Given two indices of a square triangular matrix and a size return the corresponding linear index.
  HOSTDEVICE int triangularIndex(int i, int j, int N){
    if(i>j){
      int tmp = i;
      i = j;
      j = tmp;
    }
    return j + i*(N - 1) - i*(i - 1)/2;
  }

  template<bool fixBinBIAS> struct pairDistanceCounterType;
  template<> struct pairDistanceCounterType<true >{using type= real;};
  template<> struct pairDistanceCounterType<false>{using type= ullint;};
  template<bool fixBinBIAS>
  struct PairCounterTransverser{
    using pairCountType = typename pairDistanceCounterType<fixBinBIAS>::type;
    Box3D box;
    real rcut;
    real binSize;
    pairCountType *pairDistanceCounterGPUPtr;
    bool is3D;
    int ntypes, nbins;

    PairCounterTransverser(pairCountType *pairDistanceCounterGPUPtr,
			   int ntypes, int nbins,
			   Box3D box, real rcut, real binSize, bool is3D):
      ntypes(ntypes), box(box), rcut(rcut), binSize(binSize), nbins(nbins),
      pairDistanceCounterGPUPtr(pairDistanceCounterGPUPtr), is3D(is3D){}

    HOSTDEVICE void operator ()(real4 pi, real4 pj){
      const real3 rij = box.apply_pbc(make_real3(pi) - make_real3(pj));
      const real r = sqrt(dot(rij, rij));
      if(r < rcut){
	int ti = pi.w;
	int tj = pj.w;
	if(ti>=tj){
	  int typeIndex = triangularIndex(ti, tj, ntypes);
	  addDistance(r, typeIndex);
	}
      }
    }

  private:
    template<class T>
    HOSTDEVICE void storeValue(T* address, T value){
#ifdef __CUDA_ARCH__
      atomicAdd(address, value);
#else
      *address += value;
#endif
    }

    HOSTDEVICE void addDistance(real r, int typeIndex){
      int bin = int(r/binSize);
      bin = (bin>=nbins)?(nbins-1):bin;
      const int indexInCounter = bin + nbins*typeIndex;
      if(fixBinBIAS){
	const real rbin = (bin + 0.5)*binSize;
	real norm = rbin/r;
	if(is3D){
	  norm *= norm;
	}
	storeValue((real*)pairDistanceCounterGPUPtr + indexInCounter, norm);
      }
      else{
	storeValue((ullint*)pairDistanceCounterGPUPtr + indexInCounter, 1ull);
      }
    }

  };

  real binToDistance(int i, Configuration config){
    real binSize = config.maxDistance/config.numberBins;
    return (i+0.5)*binSize;
  }

  real computeCount2rdf(const Configuration &config){
    real3 L = config.boxSize;
    real binSize = config.maxDistance/config.numberBins;
    double V = L.x*L.y;
    if(config.dimension==Configuration::dimensionality::D3) V *= L.z;
    double prefactor = M_PI;
    if(config.dimension==Configuration::dimensionality::D3) prefactor *= 2.0;
    constexpr double countedTwice = 2.0;
    double normalization = countedTwice*prefactor*binSize/V;
    return 1.0/normalization;
  }
}


#endif
