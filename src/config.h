/*Raul P. Pelaez 2017. A struct with all neded parameters
 */
#ifndef CONFIG_H
#define CONFIG_H

#define SINGLE_PRECISION

#include"vector_algebra.cuh"
#include<unistd.h>
#include<fstream>
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
    string inputFileName;
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
      if(strcmp(argv[i], "-Nsnapshots")==0)      config.numberSnapshots = atoi(argv[i+1]);
      if(strcmp(argv[i], "-outputDecimals")==0)  config.outputDecimals = atoi(argv[i+1]);
      if(strcmp(argv[i], "-dim")==0){
	if(strcmp(argv[i+1], "3D")==0)           config.dimension = Configuration::dimensionality::D3;
	else if(strcmp(argv[i+1], "2D")==0)      config.dimension = Configuration::dimensionality::D2;
	else if(strcmp(argv[i+1], "q2D")==0)     config.dimension = Configuration::dimensionality::qD2;
	else{ cerr<<"ERROR: INVALID DIMENSIONALITY IN -dim!!!"<<endl; print_help(); exit(1);}

      }
      if(strcmp(argv[i], "-device")==0){
	if(strcmp(argv[i+1], "GPU")==0)          config.deviceMode = Configuration::device::GPU;
	else if(strcmp(argv[i+1], "CPU")==0)     config.deviceMode = Configuration::device::CPU;
	else if(strcmp(argv[i+1], "auto")==0){}
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

    printf(" Raul P. Pelaez 2017.\n");
    printf(" \n");
    printf("NAME \n");
    printf("rdf -  Computes the Radial Distribution Function (RDF) of a group of positions in a file,\n");
    printf(" averages it for all snapshots in the file. \n");
    printf(" \n");
    printf("COMPILE WITH \n");
    printf(" \n");
    printf("$ nvcc  -arch=sm_52 -std=c++11 -O3 rdf.cu \n");
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
