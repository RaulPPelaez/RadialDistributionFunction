/*Raul P. Pelaez 2017-2020. A struct with all needed parameters
 */
#ifndef CONFIG_H
#define CONFIG_H

#define SINGLE_PRECISION

#include"vector_algebra.cuh"
#include<unistd.h>
#include"utils.cuh"
#include<memory>
#include<fstream>
#include<iostream>
#include"defines.h"
#include <cstdio>
#include <cstring>
namespace gdr{

  struct Configuration{
    enum device{GPU, CPU, none};
    enum dimensionality{D3, D2, qD2};
    //Default config
    dimensionality dimension=D3;
    device deviceMode=none;
#ifdef SINGLE_PRECISION
    bool doublePrecision = false;
#else
    bool doublePrecision = true;
#endif

    int outputDecimals=5;

    int numberParticles = -1;
    real3 boxSize = {0,0,0};
    real maxDistance = 0;
    int numberBins = 100;
    int numberSnapshots =-1;
    std::string inputFileName;

    //Takes into account that the distances do not have the same weight in a particular bin.
    bool fixBIAS = false;
    bool useTypes = false;
  };

  void print_help();
  //Fill config with cli arguments
  void processCommandLineArguments(char *argv[], int argc, Configuration &config){
    double Lread = 0.0;
    double Lx=0, Ly=0, Lz=0;
    double L[3];

    fori(0,argc){
      /*With -L you can have one or three numbers*/
      if(strcmp(argv[i], "-L")==0){
	Lx = strtod(argv[i+1], nullptr);
	if(argc>i+3){
	  Ly = strtod(argv[i+2], nullptr);
	  Lz = strtod(argv[i+3], nullptr);
	}
	if(!Ly || !Lz ) Lread = Lx;
      }
      if(strcmp(argv[i], "-Lx")==0)              Lx = strtod(argv[i+1], nullptr);
      if(strcmp(argv[i], "-Ly")==0)              Ly = strtod(argv[i+1], nullptr);
      if(strcmp(argv[i], "-Lz")==0)              Lz = strtod(argv[i+1], nullptr);
      if(strcmp(argv[i], "-rcut")==0)            config.maxDistance = strtod(argv[i+1], nullptr);
      if(strcmp(argv[i], "-N")==0)               config.numberParticles = atoi(argv[i+1]);
      if(strcmp(argv[i], "-nbins")==0)           config.numberBins = atoi(argv[i+1]);
      if(strcmp(argv[i], "-Nsnapshots")==0)      config.numberSnapshots = atoi(argv[i+1]);
      if(strcmp(argv[i], "-outputDecimals")==0)  config.outputDecimals = atoi(argv[i+1]);
      if(strcmp(argv[i], "-fixBIAS")==0)         config.fixBIAS = true;
      if(strcmp(argv[i], "-useTypes")==0)         config.useTypes = true;

      if(strcmp(argv[i], "-dim")==0){
	if(strcmp(argv[i+1], "3D")==0)           config.dimension = Configuration::dimensionality::D3;
	else if(strcmp(argv[i+1], "2D")==0)      config.dimension = Configuration::dimensionality::D2;
	else if(strcmp(argv[i+1], "q2D")==0)     config.dimension = Configuration::dimensionality::qD2;
	else{ std::cerr<<"ERROR: INVALID DIMENSIONALITY IN -dim!!!"<<std::endl; print_help(); exit(1);}

      }
      if(strcmp(argv[i], "-device")==0){
	if(strcmp(argv[i+1], "GPU")==0)          config.deviceMode = Configuration::device::GPU;
	else if(strcmp(argv[i+1], "CPU")==0)     config.deviceMode = Configuration::device::CPU;
	else if(strcmp(argv[i+1], "auto")==0){}
	else{ std::cerr<<"ERROR: Selected an invalid device"<<std::endl; print_help();}

      }

      if(strcmp(argv[i], "-h")==0){ print_help(); exit(0); }

    }
    if(!Lread && !(Lx&&Ly&&Lz)){std::cerr<<"ERROR!! NO VALID BOX SIZE WAS GIVEN!!"<<std::endl; print_help(); exit(1);}
    if(!config.numberParticles){std::cerr<<"ERROR!! NO VALID NUMBER PARTICLES WAS GIVEN!!"<<std::endl; print_help(); exit(1);}
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
	std::shared_ptr<std::istream> in = std::make_shared<std::ifstream>(argv[i]); //There must be a filename somewhere in the cli
	if(in->good()){
	  good_file = true;
	  config.inputFileName = std::string(argv[i]);
	  break;
	}

      }
      if(!good_file){std::cerr<< "ERROR!, NO INPUT DETECTED!!"<<std::endl; print_help(); exit(1);}
    }



  }

#include"gitversion.h"
  void print_help(){

    printf(" Raul P. Pelaez 2017-2020.\n");
    printf(" \n");
    printf("RadialDistributionFunction v%s.%s\n",
	   RadialDistributionFunction_VERSION_MAJOR,
	   RadialDistributionFunction_VERSION_MINOR);
    printf("Compiled from git commit: %s\n", GITVERSION);
    printf("NAME \n");
    printf("rdf -  Computes the Radial Distribution Function (RDF) of a group of positions in a file,\n");
    printf(" averages it for all snapshots in the file. \n");
    printf(" \n");
    printf("SYNOPSYS \n");
    printf(" \n");
    printf("rdf [OPTIONS]... [FILE]... \n");
    printf(" \n");
    printf("DESCRIPTION \n");
    printf(" Compute the Radial Distribution Function. \n");
    printf(" \n");
    printf(" With no FILE, or when file is - or not specified, reads from standard input. \n");
    printf(" \n");
    printf(" Required options: \n");
    printf(" \n");
    printf(" -N \n");
    printf(" Number of particles, all snapshots must have the same number of particles \n");
    printf(" \n");
    printf(" -L  [lx ly lz], -Lx [lx] -Ly [ly]  -Lz[lz] \n");
    printf(" Box size, positions will be folded to a box of size Lx Ly Lz. -L will make Lx= Ly = Lz = L \n");
    printf(" \n");
    printf(" -rcut \n");
    printf(" Maximum distance in the rdf, distances greater than rcut will be ignored. \n");
    printf(" \n");
    printf(" -nbins \n");
    printf(" Number of bins in the position pair histogram (binSize = rcut/nbins) \n");
    printf(" \n");
    printf(" -Nsnapshots  \n");
    printf(" Number of snapshots in the file, a snapshot must be separated from the next with a single line \n");
    printf(" \n");
    printf(" -dim [=3D] \n");
    printf(" Dimensionality of the input positions. Affects how the histogram is normalized to compute the rdf.\n");
    printf(" Can be 3D, 2D or q2D (treat as 3D, but normalize as 2D)\n");
    printf(" \n");
    printf(" -device [=auto] \n");
    printf(" Switch between GPU/CPU implementations of the algorithm. By default rdf chooses the best according to N\n");
    printf(" \n");
    printf(" -outputDecimals [=5] \n");
    printf(" Number of decimals in the output file, set through cout<<setprecision() \n");
printf("-fixBIAS\n");
printf("    This will weight the distance of a pair in a bin according to the position inside the bin (instead of weighting all distances as 1).\n");
 printf("   -useTypes\n");
 printf("       This option will interpret the fourth (third in 2D) column as particle type and will compute and output a RDF for each type pair (Ntype*(Ntype+1)/2 in total). Each RDF will start with \"# typei typej\"\n");
    printf(" \n");
    printf("FILE FORMAT \n");
    printf(" The file must have at least \"dim\" columns (the rest will be ignored) and each snapshot (including the first) \n");
    printf(" must be preceded by a line (no matter the content as long as it is a single line). See example. \n");
    printf(" \n");
    printf(" \n");
    printf("EXAMPLES: \n");
    printf(" \n");
    printf("---pos.dat---- \n");
    printf("# \n");
    printf("1 2 3 4 5 \n");
    printf("6 7 8 9 \n");
    printf("10 11 12 13 14 \n");
    printf("# \n");
    printf("15 16 17 18 \n");
    printf("19 20 21 \n");
    printf("22 23 24 \n");
    printf("------ \n");
    printf(" \n");
    printf("$ cat pos.dat | rdf -N 3 -Nsnapshots 2 -L 25 -nbins 10 -rcut 2 > rdf.dat \n");
    printf(" \n");
    printf("rdf will take the file as 2 snapshots with the positions of 3 particles in 3D.\n");

  }

}
#endif
