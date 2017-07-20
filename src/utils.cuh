/*Raul P. Pelaez 2017. Some usefull utilities

  Contains a Box struct that can apply PBC to a position

 */
#ifndef UTILS_CUH
#define UTILS_CUH
#include"vector_algebra.cuh"
#define fori(x,y) for(int i=x; i<y; i++)
#define forj(x,y) for(int j=x; j<y; j++)

namespace gdr{
  
  //A box that applies PBC on positions, can be created for 2D or 3D using real2 or real3 as template argument
  template<class vecType>
  struct Box{
    vecType L, invL;
    //Only compiled in the 3D version
    Box(real3 L):
      L(L),
      invL(make_real3(1.0/L.x, 1.0/L.y, 1.0/L.z)){
      if(L.z==real(0.0))
	invL.z = real(0.0);
    }
    //Only compiled in the 2D version
    Box(real2 L):
      L(L),
      invL(make_real2(1.0/L.x, 1.0/L.y)){
      
    }

    inline __host__  __device__ void apply_pbc(vecType &r) const{    
      r -= floorf(r*invL+real(0.5))*L; //MIC Algorithm
    }
  };
  typedef Box<real3> Box3D;
  typedef Box<real2> Box2D;



}

#endif
