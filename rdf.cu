/*Raul P. Pelaez 2017. Radial distribution function

NAME 
rdf -  Computes the Radial Distribution Function (RDF) of a group of positions in a file, averages it for all snapshots in the file.

COMPILE WITH

$ nvcc  -arch=sm_52 -std=c++11 -O3 rdf.cu

SYNOPSYS

rdf [OPTIONS]... [FILE]...

DESCRIPTION
   Compute the Radial Distribution Function.
   
   With no FILE, or when file is - or not specified, reads from standard input.

   Required options:

   -N
       Number of particles, all snapshots must have the same number of particles

   -L, -Lx [lx] -Ly [ly]  -Lz[lz]
       Box size, positions will be folded to a box of size Lx Ly Lz. -L will make Lx= Ly = Lz = L

   -rcut
       Maximum distance in the rdf, distances greater than rcut will be ignored.
   
   -nbins
       Number of bins in the position pair histogram (binSize = rcut/nbins)

   -Nsnapshots 
       Number of snapshots in the file, a snapshot must be separated from the next with a single line

   -dim [=3]
       Dimensionality of the input positions. Affects how the histogram is normalized to compute the rdf.

   -device [=GPU]
       Switch between GPU/CPU implementations of the algorithm. Currently only GPU is implemented	

FILE FORMAT
   The file must have at least "dim" columns (the rest will be ignored) and each snapshot (including the first)
   must be preceded by a line (no matter the content as long as it is a single line). See example.


EXAMPLES:

---pos.dat----
#
1 2 3 4 5
6 7 8 9
10 11 12 13 14
#
15 16 17 18
19 20 21
22 23 24
------

$ cat pos.dat | rdf -N 3 -Nsnapshots 2 -L 25 -nbins 10 -rcut 2 > rdf.dat

rdf will take the file as 2 snapshots with the positions of 3 particles in 3D.

*/

#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<string>
#include<iostream>
#include <unistd.h>
#include<vector>
#include<cmath>
#include<omp.h>
#include<memory>
#include<thrust/device_vector.h>

using namespace std;


#define fori(x,y) for(int i=x; i<y; i++)
#define forj(x,y) for(int j=x; j<y; j++)

#define SINGLE_PRECISION
#include"vector_algebra.h"

namespace gdr{

  
  struct Configuration{
    enum device{GPU, CPU};
    //Default config
    int dimension=3;
    device deviceMode=GPU;
#ifdef SINGLE_PRECISION
    bool doublePrecision = false; 
#else
    bool doublePrecision = true;
#endif
    int numberParticles = -1;
    real3 boxSize = {0,0,0};
    real maxDistance = 0;
    int numberBins = 100;
    int numberSnapshots =-1;
    string inputFileName;
  };


  //Reads numbers from a file or standard input line by line
  class InputParse{
    shared_ptr<istream> input;
    string currentLine;
  public:
    //Read from stdin by default
    InputParse(){ }

    //Take input from cin
    bool open(){
      //The lambda ensures cin is not deleted
      input.reset(&cin, [](...){});
    
      if(!input->good()){
	cerr<<"ERROR: Unable to read from stdin!"<<endl;
	return false;
      }
      return true;
    }
    //Open and read from a file
    bool open(string fileName){
      input.reset(new ifstream(fileName.c_str()));
      if(!input->good()){
	cerr<<"ERROR: Unable to open file!"<<endl;
	return false;
      }
      return true;
    }

    //Reads a line from input
    string goToNextLine(){
      if(!input) cerr<<"ERROR: No open file!"<<endl;
      getline(*input, currentLine);
      return currentLine;
    }
    //Reads a line and
    void parseNextLine(real *numbersInLine, int numberColumnsToRead){
      this->goToNextLine();
      //This could be faster
      stringstream ss;
      ss.str(currentLine);
      for(int i=0; i<numberColumnsToRead; i++){
	ss>>numbersInLine[i];
      }

    }

  };



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



  //Reads an entire frame to the second argument. Meant to be used with real4,real3 or real2
  template<class vecType>
  void readFrame(InputParse &input, vecType *pos, int numberParticles, int nCoordinatesToRead){

    //Ignore frame separator
    input.goToNextLine();

    for(int i= 0; i<numberParticles; i++){

      input.parseNextLine((real*)&pos[i], nCoordinatesToRead);

    }

  }


  //Goes through all pair of particles in the pos array, O(N^2).
  //With each pair, computes distance and sums 1 to the bin corresponding to that distance.
  /*Reference: Fast N-Body Simulation with CUDA. Chapter 31 of GPU Gems 3*/  
  template<typename BoxType, typename vectorLoadType>
  __global__ void nBody_gdrKernel(const vectorLoadType* __restrict__ pos,
				  int numTiles, /*Thread paralellism level, 
						  controls how many elements are stored in 
						  shared memory and
						  computed in parallel between synchronizations*/
				  uint N,       //Number of particles
				  BoxType box,  //Box object for PBC, can be Box3D of Box2D
				  real rcut,    //Maximum distance
				  real binSize,
				  ullint* __restrict__ pairDistanceCount
				  ){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    /*All threads must pass through __syncthreads, 
      but when N is not a multiple of 32 some threads are assigned a particle i>N.
      This threads cant return, so they are masked to not do any work*/
    bool active = true;
    if(id>=N) active = false;
    
    /*Each thread handles the interaction between particle id and all the others*/
    /*Storing blockDim.x positions in shared memory and processing all of them in parallel*/
    extern __shared__ char shMem[];

    real3 *shPos = (real3*) shMem;    
  
    real3 pi;
    if(active) {
      pi = make_real3(pos[id]); /*My position*/
    }    
    /*Distribute the N particles in numTiles tiles.
      Storing in each tile blockDim.x positions in shared memory*/
    /*This way all threads are accesing the same memory addresses at the same time*/
    for(int tile = 0; tile<numTiles; tile++){
      /*Load this tiles particles positions to shared memory*/
      const int i_load = tile*blockDim.x+threadIdx.x;
      if(i_load<N){ /*Even if im not active,
		      my thread may load a position each tile to shared memory.*/	
	shPos[threadIdx.x] = make_real3(pos[i_load]);
      }
      /*Wait for all threads to arrive*/
      __syncthreads();
      /*Go through all the particles in the current tile*/
#pragma unroll 8
      for(uint counter = 0; counter<blockDim.x; counter++){
	if(!active) break; /*An out of bounds thread must be masked*/
	int cur_j = tile*blockDim.x+counter; 
	if(cur_j<N && cur_j != id){/*If the current particle exists, compute and accumulate*/
	  /*Compute and accumulate the current pair*/
	  real3 rij = pi - shPos[counter];
	  box.apply_pbc(rij);
	  real r = sqrtf(dot(rij, rij));
	  if(r<rcut){
	    int bin = floor(r/binSize);
	    atomicAdd(&pairDistanceCount[bin], 1);
	  }
	
	}
      }/*End of particles in tile loop*/
      __syncthreads();

    }/*End of tile loop*/
  }

  class RadialFunctionDistributionGPU{
    thrust::device_vector<real4> posGPU;
    thrust::device_vector<ullint> pairDistanceCountGPU;
    int processedSnapshots = 0;
  public:
  
    RadialFunctionDistributionGPU(){}

    //Compute the pair distance histogram, meant to be called once per snapshot
    //Each time it is called the histogram is summed to the previous one
    template<class vecType>
    void processSnapshot(const vecType *posCPU, const Configuration &config){
      if(!posCPU){cerr<<"ERROR: position pointer is NULL!! in gdr GPU"<<endl;return; }

      //Recover parameters
      int N = config.numberParticles;
      Box3D box(config.boxSize);
      real rcut = config.maxDistance;
      int numberBins= config.numberBins;
      real binSize = rcut/numberBins;

      if(posGPU.size() != N) posGPU.resize(N);
      if(pairDistanceCountGPU.size() != N) pairDistanceCountGPU.resize(N, 0);

      //Get raw pointers to device memory
      auto posGPUPtr=thrust::raw_pointer_cast(posGPU.data());
      auto pairDistanceCountGPUPtr=thrust::raw_pointer_cast(pairDistanceCountGPU.data());


      
      cudaMemcpy(posGPUPtr, posCPU, N*sizeof(vecType), cudaMemcpyHostToDevice);

      //Configure and lauch kernel
      int BLOCKSIZE=128;
      int Nthreads = BLOCKSIZE<N?BLOCKSIZE:N;
      int Nblocks  = (N+Nthreads-1)/Nthreads;
      int numTiles = (N + Nthreads-1)/Nthreads;        
  
      size_t sharedMemorySize =  Nthreads*(sizeof(real3));
      
      nBody_gdrKernel<<<Nblocks, Nthreads, sharedMemorySize>>>(posGPUPtr,
							       numTiles,
							       N,
							       box,
							       rcut,
							       binSize,
							       pairDistanceCountGPUPtr);
      processedSnapshots++;
    }

    //Downloads and normalizes the pair distance histogram to compute the rdf, then overwrites gdrCPU 
    void getRadialDistributionFunction(real *rdfCPU, const Configuration &config){
      auto pairDistanceCountGPUPtr=thrust::raw_pointer_cast(pairDistanceCountGPU.data());
      //Downloads pair histogram
      thrust::host_vector<ullint> countCPU = pairDistanceCountGPU;

      
      real3 L = config.boxSize;
      int N = config.numberParticles;
      real binSize = config.maxDistance/config.numberBins;
      double V = L.x*L.y*L.z;
      if(config.dimension==2) V = L.x*L.y;
      
      double prefactor = 2.0*M_PI;
      if(config.dimension==2) prefactor *= 0.5;
      
      constexpr double countedTwice = 2;
      
      double normalization = this->processedSnapshots*countedTwice*prefactor*binSize*N*N/V;
      double R;


      for(int i=1; i<=config.numberBins; i++){
	R = (i-0.5)*binSize;
	double count2rdf = 1.0/(normalization*R);
	if(config.dimension == 3) count2rdf /= R;
	rdfCPU[i-1] = countCPU[i-1]*count2rdf;

      }
    }

    void reset(){
      this->processedSnapshots = 0;
    }
    
  };
}

using namespace gdr;

void print_help();
void processCommandLineArguments(char *argv[], int argc, Configuration &config);
int main(int argc, char *argv[]){

  //Configuration holds all parameters about the input file and gdr that might be needed.
  Configuration config;

  processCommandLineArguments(argv, argc, config);

  //InputParse handles the transformation of a line from the input file to numbers
  InputParse inputParser;

  //If called without arguments, open uses cin
  if(config.inputFileName.empty())
    inputParser.open();
  else
    inputParser.open(config.inputFileName);

  
  
  int N = config.numberParticles;
  //pos array to read a frame from the file
  std::vector<real4> pos(N, make_real4(0));
  
  std::vector<real> rdf(config.numberBins, 0);

  cerr<<"Computing.."<<endl;
  //Select between GPU/CPU implementations
  if(config.deviceMode == Configuration::device::GPU){
    
    RadialFunctionDistributionGPU rdfComputerGPU;

    for(int i=0; i<config.numberSnapshots; i++){

      readFrame<real4>(inputParser, pos.data(), N, config.dimension);

      rdfComputerGPU.processSnapshot(pos.data(), config);
    }
    //Download and normalize rdf
    rdfComputerGPU.getRadialDistributionFunction(rdf.data(), config);

    double binSize = config.maxDistance/config.numberBins;
    for(int i=0; i<config.numberBins; i++){
      double R = (i+0.5)*binSize;
      cout<<R<<" "<<rdf[i]<<endl;
    }

  }
  else if(config.deviceMode == Configuration::device::CPU){

    cerr<<"ERROR: CPU computer not implemented yet!"<<endl;
  }
  cerr<<"DONE"<<endl;


  

  
  return 0;
}




void processCommandLineArguments(char *argv[], int argc, Configuration &config){
  double Lread = 0.0;
  double Lx=0, Ly=0, Lz=0;
  double L[3];
  
  fori(0,argc){
    /*With -L you can have one or three numbers*/
    if(strcmp(argv[i], "-L")==0){
      Lx = strtod(argv[i+1], NULL);
      if(argc>i+3){
	Ly = strtod(argv[i+2], NULL);
	Lz = strtod(argv[i+3], NULL);
      }
      if(!Ly || !Lz ) Lread = Lx;
    }
    if(strcmp(argv[i], "-Lx")==0)              Lx = strtod(argv[i+1], NULL);
    if(strcmp(argv[i], "-Ly")==0)              Ly = strtod(argv[i+1], NULL);
    if(strcmp(argv[i], "-Lz")==0)              Lz = strtod(argv[i+1], NULL);
    if(strcmp(argv[i], "-rcut")==0)            config.maxDistance = strtod(argv[i+1], NULL);
    if(strcmp(argv[i], "-N")==0)               config.numberParticles = atoi(argv[i+1]);
    if(strcmp(argv[i], "-nbins")==0)           config.numberBins = atoi(argv[i+1]);
    if(strcmp(argv[i], "-Nsnapshots")==0)          config.numberSnapshots = atoi(argv[i+1]);
    if(strcmp(argv[i], "-dim")==0)             config.dimension = atoi(argv[i+1]);
    if(strcmp(argv[i], "-device")==0){
      if(strcmp(argv[i+1], "GPU")==0)          config.deviceMode = Configuration::device::GPU;
      else if(strcmp(argv[i+1], "CPU")==0)     config.deviceMode = Configuration::device::CPU;
      else{ cerr<<"ERROR: Selected an invalid device"<<endl; print_help();}
	
    }
    
    if(strcmp(argv[i], "-h")==0){ print_help(); exit(0); }

  }
  if(!Lread && !(Lx&&Ly&&Lz)){cerr<<"ERROR!! NO VALID BOX SIZE WAS GIVEN!!"<<endl; print_help(); exit(1);}
  if(!config.numberParticles){cerr<<"ERROR!! NO VALID NUMBER PARTICLES WAS GIVEN!!"<<endl; print_help(); exit(1);}
  if(!Lx||!Ly||!Lz)  L[0] = L[1] = L[2] = Lread;
  else{
    L[0] = Lx;
    L[1] = Ly;
    L[2] = Lz;
  }

  config.boxSize = make_real3(L[0], L[1], L[2]);
  
  //Look for a valid input file, in stdin or in a given filename
  if(isatty(STDIN_FILENO)){ //If there is no pipe
    bool good_file = false;
    fori(1,argc){ //exclude the exe file
      shared_ptr<istream> in = make_shared<ifstream>(argv[i]); //There must be a filename somewhere in the cli
      if(in->good()){
	good_file = true;
	config.inputFileName = string(argv[i]);
	break;	
      }

    }
    if(!good_file){cout<< "ERROR!, NO INPUT DETECTED!!"<<endl; print_help(); exit(1);}
  }


}





void print_help(){

printf("  Raul P. Pelaez 2017.                                                                                           \n");
printf("														 \n");
printf("NAME 														 \n");
printf("rdf -  Computes the Radial Distribution Function (RDF) of a group of positions in a file,                        \n");
printf("       averages it for all snapshots in the file.								 \n");
printf("														 \n");
printf("COMPILE WITH													 \n");
printf("														 \n");
printf("$ nvcc  -arch=sm_52 -std=c++11 -O3 rdf.cu									 \n");
printf("														 \n");
printf("SYNOPSYS													 \n");
printf("														 \n");
printf("rdf [OPTIONS]... [FILE]...											 \n");
printf("														 \n");
printf("DESCRIPTION													 \n");
printf("   Compute the Radial Distribution Function.									 \n");
printf("   														 \n");
printf("   With no FILE, or when file is - or not specified, reads from standard input.					 \n");
printf("														 \n");
printf("   Required options:												 \n");
printf("														 \n");
printf("   -N														 \n");
printf("       Number of particles, all snapshots must have the same number of particles				 \n");
printf("														 \n");
printf("   -L, -Lx [lx] -Ly [ly]  -Lz[lz]										 \n");
printf("       Box size, positions will be folded to a box of size Lx Ly Lz. -L will make Lx= Ly = Lz = L		 \n");
printf("														 \n");
printf("   -rcut													 \n");
printf("       Maximum distance in the rdf, distances greater than rcut will be ignored.				 \n");
printf("   														 \n");
printf("   -nbins													 \n");
printf("       Number of bins in the position pair histogram (binSize = rcut/nbins)					 \n");
printf("														 \n");
printf("   -Nsnapshots 													 \n");
printf("       Number of snapshots in the file, a snapshot must be separated from the next with a single line		 \n");
printf("														 \n");
printf("   -dim [=3]													 \n");
printf("       Dimensionality of the input positions. Affects how the histogram is normalized to compute the rdf.	 \n");
printf("														 \n");
printf("   -device [=GPU]												 \n");
printf("       Switch between GPU/CPU implementations of the algorithm. Currently only GPU is implemented		 \n");
printf("														 \n");
printf("FILE FORMAT													 \n");
printf("   The file must have at least \"dim\" columns (the rest will be ignored) and each snapshot (including the first)	 \n");
printf("   must be preceded by a line (no matter the content as long as it is a single line). See example.		 \n");
printf("														 \n");
printf("														 \n");
printf("EXAMPLES:													 \n");
printf("														 \n");
printf("---pos.dat----													 \n");
printf("#														 \n");
printf("1 2 3 4 5													 \n");
printf("6 7 8 9														 \n");
printf("10 11 12 13 14													 \n");
printf("#														 \n");
printf("15 16 17 18													 \n");
printf("19 20 21													 \n");
printf("22 23 24													 \n");
printf("------														 \n");
printf("														 \n");
printf("$ cat pos.dat | rdf -N 3 -Nsnapshots 2 -L 25 -nbins 10 -rcut 2 > rdf.dat					 \n");
printf("														 \n");
printf("rdf will take the file as 2 snapshots with the positions of 3 particles in 3D.                                   \n");



}




