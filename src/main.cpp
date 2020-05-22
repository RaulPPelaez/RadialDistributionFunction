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

   -useTypes
       This option will interpret the fourth (third in 2D) column as particle type and will compute and output a RDF for each type pair (Ntype*(Ntype+1)/2 in total). Each RDF will start with "# typei typej"

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
std::vector<real2> computeWithGPU(Configuration config);

template<bool fixBIAS>
std::vector<real2> computeWithCPU(Configuration config);

std::vector<real2> computeRDF(Configuration config){
  if(config.deviceMode == Configuration::device::GPU){
    if(config.fixBIAS)
      return computeWithGPU<true>(config);
    else
      return computeWithGPU<false>(config);
  }
  else if(config.deviceMode == Configuration::device::CPU){
    if(config.fixBIAS)
      return computeWithCPU<true>(config);
    else
      return computeWithCPU<false>(config);
  }
  else{
    cerr<<"Invalid configuration"<<endl;
    exit(1);
    return std::vector<real2>();
  }
}

int Ntypes;
void writeRDF(const std::vector<real2> &rdf, Configuration config){
  cout<<std::setprecision(config.outputDecimals);
  for(int typei=0; typei<Ntypes; typei++){
    for(int typej = typei; typej<Ntypes; typej++){
      std::cout<<"# "<<typei<<" "<<typej<<std::endl;
      int typeIndex = triangularIndex(typei, typej, Ntypes)*config.numberBins;
      for(int i=0; i<config.numberBins; i++){
	double R = binToDistance(i, config);
	std::cout<<R<<" "<<rdf[typeIndex + i].x<<" "<<rdf[typeIndex + i].y<<std::endl;
      }
    }
  }
}

Configuration::device chooseDevice(Configuration config){
  Configuration::device dev = config.deviceMode;
  if(dev == Configuration::device::none){
#ifdef GPUMODE
    if(config.numberParticles > 500) dev = Configuration::device::GPU;
    else dev = Configuration::device::CPU;
#else
    dev = Configuration::device::CPU;
#endif
  }
  return dev;
}

int main(int argc, char *argv[]){
  Configuration config;
  processCommandLineArguments(argv, argc, config);
  config.deviceMode = chooseDevice(config);
  if(config.dimension == Configuration::dimensionality::D2){
    config.boxSize.z = 0;
  }
  auto rdf = computeRDF(config);
  writeRDF(rdf, config);
  return 0;
}

std::vector<int> countParticlesPerType(const real4* pos, Configuration config){
  int foundTypes = 1;
  std::vector<int> numberParticlesPerType(1, 0);
  for(int i = 0; i<config.numberParticles; i++){
    int typei = int(pos[i].w + 0.5);
    if(typei >= foundTypes){
      foundTypes++;
      numberParticlesPerType.resize(foundTypes, 0);
    }
    numberParticlesPerType[typei]++;
  }
  return numberParticlesPerType;
}

template<class RDFComputer>
std::vector<real2> computeWith(Configuration config){
  InputParse inputParser(config.inputFileName);
  int numberCoordinatesPerParticle = 3;
  if(config.dimension == Configuration::dimensionality::D2){
    numberCoordinatesPerParticle = 2;
  }
  int N = config.numberParticles;
  std::vector<real4> pos(N, real4());
  readFrame(inputParser, pos.data(), N, numberCoordinatesPerParticle + config.useTypes);
  auto numberParticlesPerType = countParticlesPerType(pos.data(), config);
  Ntypes = numberParticlesPerType.size();
  RDFComputer rdfComputer(config, numberParticlesPerType);
  rdfComputer.processSnapshot(pos.data());
  for(int i=1; i<config.numberSnapshots; i++){
    readFrame(inputParser, pos.data(), N, numberCoordinatesPerParticle + config.useTypes);
    rdfComputer.processSnapshot(pos.data());
  }
  auto rdf = rdfComputer.getRadialDistributionFunction();
  return std::move(rdf);
}

template<bool fixBIAS>
std::vector<real2> computeWithGPU(Configuration config){
  #ifdef GPUMODE
  return std::move(computeWith<RadialDistributionFunctionGPU<fixBIAS>>(config));
  #else
  cerr<<"ERROR: Compiled in CPU mode only"<<endl;
  exit(1);
  return std::vector<real2>();
  #endif
}

template<bool fixBIAS>
std::vector<real2> computeWithCPU(Configuration config){
  return std::move(computeWith<RadialDistributionFunctionCPU<fixBIAS>>(config));
}

