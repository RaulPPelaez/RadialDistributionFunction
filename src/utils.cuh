/*Raul P. Pelaez 2017. Some usefull utilities

  Contains a Box struct that can apply PBC to a position

 */
#ifndef UTILS_CUH
#define UTILS_CUH
#include"vector_algebra.cuh"
#define fori(x,y) for(int i=x; i<y; i++)
#define forj(x,y) for(int j=x; j<y; j++)

#ifdef GPUMODE
#define HOSTDEVICE __host__ __device__
#else
#define HOSTDEVICE
#endif
namespace gdr{

  //A box that applies PBC on positions, can be created for 2D or 3D using real2 or real3 as template argument
  template<class vecType>
  struct Box{
    vecType boxSize, invBoxSize;
    Box():Box(make_real3(0)){}
    //Only compiled in the 3D version
    Box(real3 boxSize):
      boxSize(boxSize),
      invBoxSize(make_real3(1.0/boxSize.x, 1.0/boxSize.y, 1.0/boxSize.z)){
      if(boxSize.z==real(0.0))
	invBoxSize.z = real(0.0);
    }
    //Only compiled in the 2D version
    Box(real2 boxSize):
      boxSize(boxSize),
      invBoxSize(make_real2(1.0/boxSize.x, 1.0/boxSize.y)){

    }

    // inline HOSTDEVICE void apply_pbc(vecType &r) const{
    //   r -= floorf(r*invBoxSize+real(0.5))*boxSize; //MIC Algorithm
    // }
    inline HOSTDEVICE vecType apply_pbc(const vecType &r) const{
      return r - floorf(r*invBoxSize+real(0.5))*boxSize; //MIC Algorithm
    }

  };
  typedef Box<real3> Box3D;
  typedef Box<real2> Box2D;


  struct Grid{
    Box3D box;
    /*A magic vector that transforms cell coordinates to 1D index when dotted*/
    /*Simply: 1, ncellsx, ncellsx*ncellsy*/
    int3 gridPos2CellIndex;

    int3 cellDim; //ncells in each size
    real3 cellSize;
    real3 invCellSize; /*The inverse of the cell size in each direction*/
    Grid(): Grid(Box3D(), make_int3(0,0,0)){}
    Grid(Box3D box, int3 cellDim):
	box(box),
	cellDim(cellDim){

	cellSize = box.boxSize/make_real3(cellDim);
	invCellSize = 1.0/cellSize;
	if(box.boxSize.z == real(0.0)) invCellSize.z = 0;

	gridPos2CellIndex = make_int3( 1,
				       cellDim.x,
				       cellDim.x*cellDim.y);

    }
    template<class VecType>
    inline HOSTDEVICE int3 getCell(const VecType &r) const{
	// return  int( (p+0.5L)/cellSize )
      int3 cell = make_int3((      box.apply_pbc(make_real3(r)) + real(0.5)*box.boxSize)*invCellSize);
	//Anti-Traquinazo guard, you need to explicitly handle the case where a particle
	// is exactly at the box limit, AKA -L/2. This is due to the precision loss when
	// casting int from floats, which gives non-correct results very near the cell borders.
	// This is completly neglegible in all cases, except with the cell 0, that goes to the cell
	// cellDim, which is catastrophic.
	//Doing the previous operation in double precision (by changing 0.5f to 0.5) also works, but it is a bit of a hack and the performance appears to be the same as this.
	//TODO: Maybe this can be skipped if the code is in double precision mode
	if(cell.x==cellDim.x) cell.x = 0;
	if(cell.y==cellDim.y) cell.y = 0;
	if(cell.z==cellDim.z) cell.z = 0;
	return cell;
    }

    inline HOSTDEVICE int getCellIndex(const int3 &cell) const{
	return dot(cell, gridPos2CellIndex);
    }

    inline HOSTDEVICE int3 pbc_cell(const int3 &cell) const{
	int3 cellPBC;
	cellPBC.x = pbc_cell_coord<0>(cell.x);
	cellPBC.y = pbc_cell_coord<1>(cell.y);
	cellPBC.z = pbc_cell_coord<2>(cell.z);
	return cellPBC;
    }

    template<int coordinate>
    inline HOSTDEVICE int pbc_cell_coord(int cell) const{
	int ncells = 0;
	if(coordinate == 0){
	  ncells = cellDim.x;
	}
	if(coordinate == 1){
	  ncells = cellDim.y;
	}

	if(coordinate == 2){
	  ncells = cellDim.z;
	}

	if(cell <= -1) cell += ncells;
	else if(cell >= ncells) cell -= ncells;
	return cell;
    }

  };



}

#endif
