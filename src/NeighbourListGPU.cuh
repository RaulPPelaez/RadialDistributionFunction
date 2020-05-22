/*Raul P. Pelaez 2017-2020. Cell List implementation for the GPU

 */
#ifndef NEIGHBOURLISTGPU_CUH
#define NEIGHBOURLISTGPU_CUH

#include"vector_algebra.cuh"
#include"config.h"
#include"utils.cuh"

#include"ParticleSorter.cuh"
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include"third_party/cub/cub.cuh"

#include<limits>

namespace gdr{

  template<class T, class OutputIterator>
  __global__ void fillWithGPU(OutputIterator array, T value, int N){
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    if(id>=N) return;

    array[id] = value;
  }

  template<class T, class OutputIterator, class Iterator>
  __global__ void fillWithGPU(OutputIterator array, Iterator indexIterator, T value, int N){
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    if(id>=N) return;
    int i = indexIterator[id];
    array[i] = value;
  }

  namespace CellList_ns{
    constexpr int EMPTY_CELL = std::numeric_limits<int>::max();

    template<class InputIterator>
    __global__ void fillCellList(InputIterator sortPos,
				 int *cellStart, int *cellEnd,
				 int N, Grid grid){
      uint id = blockIdx.x*blockDim.x + threadIdx.x;
      if(id<N){
	uint icell, icell2;
	icell = grid.getCellIndex(grid.getCell(make_real3(sortPos[id])));
	if(id>0){
	  icell2 = grid.getCellIndex(grid.getCell(make_real3(sortPos[id-1])));
	}
	else
	  icell2 = 0;
	if(id ==0 || icell != icell2){
	  cellStart[icell] = id;
	  if(id>0)
	    cellEnd[icell2] = id;
	}
	if(id == N-1) cellEnd[icell] = N;
      }

    }


    template<class Transverser, class InputIterator>
    __global__ void transverseCellList(Transverser tr,
  				       InputIterator sortPos,
  				       const int *sortedIndex,
  				       const int * __restrict__ cellStart,
				       const int * __restrict__ cellEnd,
  				       real cutOff2,
  				       int N, Grid grid){
      int id = blockIdx.x*blockDim.x + threadIdx.x;
      if(id>=N) return;

      const real4 myParticle = sortPos[id];
      int3 cellj;
      const int3 celli = grid.getCell(myParticle);

      int zi = -1;
      int zf = 1;
      if(grid.cellDim.z == 1){
  	zi = zf = 0;
      }

      for(int x=-1; x<=1; x++){
  	cellj.x = grid.pbc_cell_coord<0>(celli.x + x);
  	for(int z=zi; z<=zf; z++){
  	  cellj.z = grid.pbc_cell_coord<2>(celli.z + z);
  	  for(int y=-1; y<=1; y++){
  	    cellj.y = grid.pbc_cell_coord<1>(celli.y + y);
  	    const int icell  = grid.getCellIndex(cellj);
  	    const int firstParticle = cellStart[icell];
  	    if(firstParticle != EMPTY_CELL){ /*Continue only if there are particles in this cell*/
  	      /*Index of the last particle in the cell's list*/
  	      const int lastParticle = cellEnd[icell];
  	      const int nincell = lastParticle-firstParticle;

  	      for(int j=0; j<nincell; j++){
  		int cur_j = j + firstParticle;// sortedIndex[j+firstParticle];
  		if(cur_j < N && cur_j > id){
		  tr(myParticle, sortPos[cur_j]);
  		}//endif
  	      }//endfor
  	    }//endif
  	  }//endfor y
  	}//endfor z
      }//endfor x
    }

  }



  class CellList{
  protected:
    thrust::device_vector<int> cellStart, cellEnd;

    thrust::device_vector<real4> sortPos;

    ParticleSorter ps;  //hash sort handler

    Configuration config;
  public:
    CellList(){ }
    ~CellList(){ }

    //Use a transverser to transverse the list using directly the cell list (without constructing a neighbour list)
    template<class Transverser>
    void transverseList(Transverser &tr, cudaStream_t st = 0){
      int numberParticles = config.numberParticles;


      int3 cellDim = make_int3(config.boxSize/config.maxDistance + real(0.5));
      if(cellDim.z == 0) cellDim.z = 1;
      Grid grid(Box3D(config.boxSize), cellDim);

      int Nthreads=128;
      int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

      auto sortPos_ptr = thrust::raw_pointer_cast(sortPos.data());
      auto cellStart_ptr = thrust::raw_pointer_cast(cellStart.data());
      auto cellEnd_ptr = thrust::raw_pointer_cast(cellEnd.data());

      CellList_ns::transverseCellList<<<Nblocks, Nthreads, 0, st>>>(tr,
			sortPos_ptr,
			ps.getSortedIndexArray(numberParticles),
			cellStart_ptr, cellEnd_ptr,
			config.maxDistance*config.maxDistance,
			numberParticles,
			grid);

    }
    template<class vecType>
    void updateNeighbourList(const vecType *posGPU, const Configuration &config,  cudaStream_t st = 0){
      this->config = config;
      int numberParticles = config.numberParticles;
      Box3D box(config.boxSize);
      int3 cellDim = make_int3(config.boxSize/config.maxDistance + real(0.5));
      if(cellDim.z == 0) cellDim.z = 1;
      int ncells = cellDim.x*cellDim.y*cellDim.z;
      if(cellStart.size()!= ncells) cellStart.resize(ncells);
      if(cellEnd.size()!= ncells) cellEnd.resize(ncells);
      auto cellStart_ptr = thrust::raw_pointer_cast(cellStart.data());
      auto cellEnd_ptr = thrust::raw_pointer_cast(cellEnd.data());
      cub::CountingInputIterator<int> it(0);
      int Nthreads=512;
      int Nblocks=ncells/Nthreads + ((ncells%Nthreads)?1:0);
      fillWithGPU<<<Nblocks, Nthreads, 0, st>>>(cellStart_ptr, it,
						CellList_ns::EMPTY_CELL, ncells);
      ps.updateOrderByCellHash<Sorter::MortonHash>(posGPU,
						   numberParticles,
						   box, cellDim, st);
      sortPos.resize(numberParticles);
      ps.applyCurrentOrder(posGPU, sortPos.begin(), numberParticles, st);
      auto sortPos_ptr = thrust::raw_pointer_cast(sortPos.data());
      Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);
      Grid grid(box, cellDim);
      CellList_ns::fillCellList<<<Nblocks, Nthreads, 0, st>>>(sortPos_ptr, //posGroupIterator,
       							      cellStart_ptr,
       							      cellEnd_ptr,
       							      numberParticles,
       							      grid);
    }


  };


}
#endif
