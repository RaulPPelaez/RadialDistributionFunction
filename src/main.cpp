/*Raul P. Pelaez 2017-2020. Radial distribution function

NAME
rdf -  Computes the Radial Distribution Function (RDF) of a group of positions in a file, averages it for all snapshots in the file.

COMPILE WITH

$ mkdir build; cd build; cmake ..; make

INSTALL WITH

$ make install

SYNOPSYS

rdf [OPTIONS]... [FILE]...

DESCRIPTION
   Compute the Radial Distribution Function.

   With no FILE, or when file is - or not specified, reads from standard input.

   Required options:

   -N
       Number of particles, all snapshots must have the same number of particles

   -L [lx ly lz], -Lx [lx] -Ly [ly]  -Lz[lz]
       Box size, positions will be folded to a box of size Lx Ly Lz. -L will make Lx= Ly = Lz = L

   -rcut
       Maximum distance in the rdf, distances greater than rcut will be ignored.

   -nbins
       Number of bins in the position pair histogram (binSize = rcut/nbins)

   -Nsnapshots
       Number of snapshots in the file, a snapshot must be separated from the next with a single line

   -dim [=3D]
       Dimensionality of the input positions. Affects how the histogram is normalized to compute the rdf.
       Can be 3D, 2D or q2D (treat as 3D, but normalize as 2D)

   -device [=GPU]
       Switch between GPU/CPU implementations of the algorithm. Currently only GPU is implemented

   -outputDecimals [=5]
       Number of decimals in the output file, set through cout<<setprecision()

   -fixBIAS
       This will weight the distance of a pair in a bin according to the position inside the bin (instead of weighting all distances as 1).

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

#include<string>
#include<iostream>
#include<memory>


#define SINGLE_PRECISION
#include"vector_algebra.cuh"
#include"utils.cuh"
#include"config.h"
#include"inputFast.h"
#ifdef GPUMODE
#include"rdfGPU.cuh"
#endif
#include"rdfCPU.h"
#include<iomanip>
using namespace gdr;
using std::cerr;
using std::endl;
using std::cout;

template<bool fixBIAS>
void computeWithGPU(InputParse &inputParser, const Configuration &config, int numberCoordinatesPerParticle);

template<class vecType, bool fixBIAS>
void computeWithCPU(InputParse &inputParser, const Configuration &config, int numberCoordinatesToRead);

int main(int argc, char *argv[]){
  Configuration config;
  processCommandLineArguments(argv, argc, config);
  if(config.deviceMode == Configuration::device::none){
#ifdef GPUMODE
    if(config.numberParticles > 500) config.deviceMode = Configuration::device::GPU;
    else config.deviceMode = Configuration::device::CPU;
#else
    config.deviceMode = Configuration::device::CPU;
#endif
  }
  InputParse inputParser;
  if(config.inputFileName.empty())
    inputParser.open();
  else
    inputParser.open(config.inputFileName);
  int numberCoordinatesPerParticle = 3;
  if(config.dimension == Configuration::dimensionality::D2){
    numberCoordinatesPerParticle = 2;
    config.boxSize.z = 0;
  }
  cout<<std::setprecision(config.outputDecimals);
  if(config.deviceMode == Configuration::device::GPU){
    if(config.fixBIAS)
      computeWithGPU<true>(inputParser, config, numberCoordinatesPerParticle);
    else
      computeWithGPU<false>(inputParser, config, numberCoordinatesPerParticle);
  }
  else if(config.deviceMode == Configuration::device::CPU){
    if(config.fixBIAS)
      computeWithCPU<real3, true>(inputParser, config, numberCoordinatesPerParticle);
    else
      computeWithCPU<real3, false>(inputParser, config, numberCoordinatesPerParticle);
  }
  return 0;
}



template<bool fixBIAS>
void computeWithGPU(InputParse &inputParser, const Configuration &config, int numberCoordinatesPerParticle){
  #ifdef GPUMODE
  int N = config.numberParticles;
  std::vector<real> rdf(config.numberBins, 0);
  std::vector<real> std(config.numberBins, 0);
  RadialDistributionFunctionGPU<fixBIAS> rdfComputerGPU;
  std::vector<real4> pos(N, make_real4(0));
  for(int i=0; i<config.numberSnapshots; i++){
    readFrame(inputParser, pos.data(), N, numberCoordinatesPerParticle);
    rdfComputerGPU.processSnapshot(pos.data(), config);
  }
  rdfComputerGPU.getRadialDistributionFunction(rdf.data(), std.data(), config);
  double binSize = config.maxDistance/config.numberBins;
  for(int i=0; i<config.numberBins; i++){
    double R = (i+0.5)*binSize;
    std::cout<<std::setprecision(2*sizeof(real))<<R<<" "<<rdf[i]<<" "<<std[i]<<std::endl;
  }
  #else
  cerr<<"ERROR: Compiled in CPU mode only"<<endl;
  exit(1);
  #endif
}

template<class vecType, bool fixBIAS>
void computeWithCPU(InputParse &inputParser, const Configuration &config, int numberCoordinatesPerParticle){
  int N = config.numberParticles;
  std::vector<real> rdf(config.numberBins, 0);
  std::vector<real> std(config.numberBins, 0);
  RadialDistributionFunctionCPU<fixBIAS> rdfComputerCPU;
  std::vector<vecType> pos(N);
  memset(pos.data(), 0, N*sizeof(vecType));
  for(int i=0; i<config.numberSnapshots; i++){
    readFrame(inputParser, pos.data(), N, numberCoordinatesPerParticle);
    rdfComputerCPU.processSnapshot(pos.data(), config);
  }
  rdfComputerCPU.getRadialDistributionFunction(rdf.data(), std.data(), config);
  double binSize = config.maxDistance/config.numberBins;
  for(int i=0; i<config.numberBins; i++){
    double R = (i+0.5)*binSize;
    cout<<std::setprecision(2*sizeof(real))<<R<<" "<<rdf[i]<<" "<<std[i]<<endl;
  }
}
