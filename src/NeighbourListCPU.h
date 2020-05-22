/* Raul P. Pelaez 2017-2020. Neighbour list for the CPU.
   This file implements the head and list cell list algorithm.
 */
#ifndef NEIGHBOURLISTCPU_H
#define NEIGHBOURLISTCPU_H

#include"vector_algebra.cuh"
#include <algorithm>
#include<vector>
#include"config.h"
#include"utils.cuh"
namespace gdr{
  class NeighbourListCPU{
    std::vector<int> head, list;
    std::vector<real4> storedPos;
    int3 numberCells;
  public:
    NeighbourListCPU(){}

    bool shouldUse(const Configuration &config){
      if(config.numberParticles<500) return false;
      int3 ncells = make_int3(config.boxSize/config.maxDistance +0.5);
      if(ncells.x<=3 || ncells.y<=3){
	return false;
      }
      else if(config.dimension!=Configuration::dimensionality::D2 && ncells.z<=3){
	return false;
      }
      return true;
    }

    template<class Iterator>
    void makeList(Iterator pos, const Configuration &config){
      int numberParticles = config.numberParticles;
      storedPos.resize(numberParticles);
      std::copy(pos, pos + numberParticles, storedPos.begin());
      real rcut = config.maxDistance;
      Box3D box(config.boxSize);
      int3 ncells = make_int3(config.boxSize/rcut +0.5);
      if(ncells.z==0) ncells.z= 1;
      Grid grid(box, ncells);
      this->numberCells = ncells;
      int totalCells = ncells.x*ncells.y*ncells.z+1;
      if(head.size() != totalCells ) head.resize(totalCells);
      if(list.size() != numberParticles+1) list.resize(numberParticles+1);
      std::fill(head.begin(), head.end(), 0);
      int icell;
      real3 temppos;
      for(int i=1; i<=numberParticles; i++){
	temppos =   box.apply_pbc(make_real3(storedPos[i-1]));
	icell = grid.getCellIndex(grid.getCell(temppos));
	list[i] = head[icell];
	head[icell] = i;
      }
    }

    template<class PairFunctor>
    void transverseList(PairFunctor &transverser, const Configuration &config){
      int numberParticles = config.numberParticles;
      Box3D box(config.boxSize);
      Grid grid(box, numberCells);
      for(int i=0; i<numberParticles;i++){
	real3 posindex;
	posindex =  box.apply_pbc(make_real3(storedPos[i]));
	int3 cell;
	cell = grid.getCell(posindex);
	int j;
	bool is2D = config.dimension==Configuration::dimensionality::D2;
	const int nneighbours = is2D?9:27;
	for(int ic=0; ic<nneighbours; ic++){
	  int3 cellj = cell;
	  cellj.x += ic%3-1;
	  cellj.y += (ic/3)%3-1;
	  cellj.z =  is2D?0:(cell.z + ic/9-1);
	  cellj = grid.pbc_cell(cellj);
	  int jcel = grid.getCellIndex(cellj);
	  j = head[jcel];
	  if(j==0) continue;
	  do{
	    if(i < (j-1)){
	      transverser(storedPos[i], storedPos[j-1]);
	    }
	    j = list[j];
	  }while(j!=0);
	}
      }
    }
  };
}
#endif
