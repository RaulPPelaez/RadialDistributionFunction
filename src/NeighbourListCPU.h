
#ifndef NEIGHBOURLISTCPU_H
#define NEIGHBOURLISTCPU_H

namespace gdr{
  class NeighbourListCPU{
    std::vector<int> head, list;
    int3 numberCells;
  public:
    NeighbourListCPU(){}
    //You should not create a cell list with less than 3 cells per dimension
    //Or with too little particles
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
    template<class vecType>
    void makeList(vecType *pos, const Configuration &config){
      
      int N = config.numberParticles;

      real rcut = config.maxDistance;

      Box3D box(config.boxSize);

      int3 ncells = make_int3(config.boxSize/rcut +0.5);
      if(ncells.z==0) ncells.z= 1;
      Grid grid(box, ncells);
      this->numberCells = ncells;
      int totalCells = ncells.x*ncells.y*ncells.z+1;
      if(head.size() != totalCells ) head.resize(totalCells);
      if(list.size() != N+1) list.resize(N+1);
      
      std::fill(head.begin(), head.end(), 0);

      int icell;   //Cell index in head
      real3 temppos; //position of a particle, we dont apply PBC to pos variable
      /*For every particle in the system (carefull with the indices!)*/
      //We go from 1 to N, but address the particles as 0 to N-1.
      for(int i=1; i<=N; i++){
	/*Save the (i-1)th position to tempos*/
	temppos =   box.apply_pbc(make_real3(pos[i-1]));
	/*And reduce it to the primary box*/

	/*Compute the cell coordinates of particle i-1, see getcell below!*/
	/*Compute the head index of cell[] (Look in the notes!)*/
	icell = grid.getCellIndex(grid.getCell(temppos));
	/*Add particle to head and list (Look in the notes!)*/
	list[i] = head[icell];
	head[icell] = i;
      }    

    }
    template<class PairFunctor, class vectorType>
    void transverseList(const vectorType* pos, PairFunctor &transverser, const Configuration &config){
      int N = config.numberParticles;
      Box3D box(config.boxSize);
      Grid grid(box, numberCells);
      for(int i=0; i<N;i++){
	/*Save the index particle position*/
	real3 posindex;//Temporal position storage
	/*Get it to the primary box*/
	posindex =  box.apply_pbc(make_real3(pos[i]));


	int3 cell;   //Cell coordinates of the index particle
	cell = grid.getCell(posindex);

	int j; //Index of neighbour particle
	int jcel, jcelx, jcely, jcelz; //cell coordinates and cell index for particle j
	int maxZcell = cell.z+1;
	int minZcell = cell.z-1;
	if(config.dimension==Configuration::dimensionality::D2){
	  maxZcell = minZcell = 1;
	}
	/*For every neighbouring cell (26 cells in 3D)*/
	for(int jz=minZcell; jz<=maxZcell;jz++){
	  /*The neighbour cell must take into account pbc! (see pbc_cells!)*/
	  jcelz = grid.pbc_cell_coord<2>(jz);
	  for(int jx=cell.x-1; jx<=cell.x+1;jx++){
	    jcelx = grid.pbc_cell_coord<0>(jx);
	    for(int jy=cell.y-1; jy<=cell.y+1;jy++){
	      jcely = grid.pbc_cell_coord<1>(jy);
	      //See getcell!
	      jcel = grid.getCellIndex(make_int3(jcelx, jcely, jcelz));
	      /*Get the highest index particle in cell jcel*/
	      j = head[jcel];
	      /*If there is no particles go to the next cell*/
	      if(j==0) continue;
	      /*Use list to travel through all the particles, j, in cell jcel*/
	      do{
		/*Add the energy of the pair interaction*/
		//Be careful not to compute one particle with itself!, j-1 because of head and list indexes!
		if(i!=(j-1))
		  transverser(i, j-1);
		j=list[j];
		/*When j=0 (list[j] = 0) then there is no more particles in cell jcel (see the notes!)*/
	      }while(j!=0);
	    } //End jz
	  } //End jy
	} //End jx
      }//End particle loop
    } //End function
    
  };
}
#endif
