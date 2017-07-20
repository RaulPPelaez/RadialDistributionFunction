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

   -L [lx ly lz], -Lx [lx] -Ly [ly]  -Lz[lz]
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

#include<string>
#include<iostream>
#include<memory>

using namespace std;

#define SINGLE_PRECISION
#include"vector_algebra.cuh"
#include"utils.cuh"
#include"config.h"
#include"input.h"
#include"rdfGPU.cuh"
#include"rdfCPU.h"
#include<iomanip>
using namespace gdr;

template<class vecType>
void computeWithCPU(InputParse &inputParser, const Configuration &config, int numberCoordinatesToRead);
int main(int argc, char *argv[]){

  //Configuration holds all parameters about the input file and gdr that might be needed.
  Configuration config;

  //Fill config with cli arguments
  processCommandLineArguments(argv, argc, config);

  //If device is automatic, use this rule of hand to select which one to use
  if(config.deviceMode == Configuration::device::none){
    if(config.numberParticles > 2e3) config.deviceMode = Configuration::device::GPU;
    else config.deviceMode = Configuration::device::CPU;
  }
  //InputParse handles the transformation of a line from the input file to numbers
  InputParse inputParser;

  //If called without arguments, open uses cin
  if(config.inputFileName.empty())
    inputParser.open();
  else
    inputParser.open(config.inputFileName);

  int numberCoordinatesPerParticle = 3;
  if(config.dimension == Configuration::dimensionality::D2){
    numberCoordinatesPerParticle = 2;
    config.boxSize.z = 0;
  }  


  cerr<<"Computing.."<<endl;
  cout<<setprecision(config.outputDecimals);
  //Select between GPU/CPU implementations
  if(config.deviceMode == Configuration::device::GPU){
    int N = config.numberParticles;
    std::vector<real> rdf(config.numberBins, 0);
    //Standard deviation
    std::vector<real> std(config.numberBins, 0);
    
    RadialDistributionFunctionGPU rdfComputerGPU;
    //pos array to read a frame from the file. real4 really improves GPU efficiency 
    std::vector<real4> pos(N, make_real4(0));

    for(int i=0; i<config.numberSnapshots; i++){
      //In 2D the 3rd coordinate is never read and thus is always 0.
      readFrame(inputParser, pos.data(), N, numberCoordinatesPerParticle);

      rdfComputerGPU.processSnapshot(pos.data(), config);
    }
    //Download and normalize rdf
    rdfComputerGPU.getRadialDistributionFunction(rdf.data(), std.data(), config);

    //Print
    double binSize = config.maxDistance/config.numberBins;
    for(int i=0; i<config.numberBins; i++){
      double R = (i+0.5)*binSize;
      cout<<R<<" "<<rdf[i]<<" "<<std[i]<<endl;
    }

  }
  else if(config.deviceMode == Configuration::device::CPU){
    computeWithCPU<real3>(inputParser, config, numberCoordinatesPerParticle);
  }
  cerr<<"DONE"<<endl;
  
  return 0;
}



template<class vecType>
void computeWithCPU(InputParse &inputParser, const Configuration &config, int numberCoordinatesPerParticle){
  int N = config.numberParticles;

  std::vector<real> rdf(config.numberBins, 0);
  //Standard deviation
  std::vector<real> std(config.numberBins, 0);

  RadialDistributionFunctionCPU rdfComputerCPU;
    
  //pos array to read a frame from the file. real4 really improves GPU efficiency    
  std::vector<vecType> pos(N);
  memset(pos.data(), 0, N*sizeof(vecType));
  for(int i=0; i<config.numberSnapshots; i++){
    //In 2D the 3rd coordinate is never read and thus is always 0.
    readFrame(inputParser, pos.data(), N, numberCoordinatesPerParticle);

    rdfComputerCPU.processSnapshot(pos.data(), config);
  }
  //Download and normalize rdf
  rdfComputerCPU.getRadialDistributionFunction(rdf.data(), std.data(), config);

  //Print
  double binSize = config.maxDistance/config.numberBins;
  for(int i=0; i<config.numberBins; i++){
    double R = (i+0.5)*binSize;
    cout<<R<<" "<<rdf[i]<<" "<<std[i]<<endl;
  }
}
