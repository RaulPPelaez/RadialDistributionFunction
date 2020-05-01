/*Raul P. Pelaez 2017-2020. Contains utility classes/functions to read snapshots from a file
  This file differs from input.h in the way a line is processed, this uses strtod instead of stringstream, it speeds up the reading process by 300%, which results in an overall speed up of 600% when using GPU (when the computations is highly IO bound, small rcut, high number of particles).
 */
#ifndef INPUT_H
#define INPUT_H

#include"vector_algebra.cuh"
#include<fstream>
#include<sstream>
#include<string>
#include"superIO.h"
#include <memory>
#include<iostream>
namespace gdr{
  //Reads numbers from a file or standard input line by line
  class InputParse{
    std::shared_ptr<superIO::superInputFile> in;
    char * currentLine = nullptr;
    size_t linesize = 0;
  public:

    InputParse(){ }

    bool open(){
      in = std::make_shared<superIO::superInputFile>();
      return true;
    }

    bool open(std::string fileName){
      in = std::make_shared<superIO::superInputFile>(fileName);
      return true;
    }

    bool good(){
      return in->good();
    }

    char* goToNextLine(){
      linesize = in->getNextLine(currentLine);
      return currentLine;
    }

    template<class Iterator>
    void parseNextLine(Iterator numbersInLine, int numberColumnsToRead){
      this->goToNextLine();
      superIO::string2numbers(currentLine, linesize, numberColumnsToRead, numbersInLine);
    }

  };


}
#endif

namespace gdr{

  template<class vecType>
  void readFrame(InputParse &input, vecType *pos, int numberParticles, int nCoordinatesToRead){
    if(!input.good()){
      std::cerr<<"ERROR: No open file!"<<std::endl;
      exit(1);
    }
    input.goToNextLine();
    for(int i= 0; i<numberParticles; i++){
      input.parseNextLine((real*)&pos[i], nCoordinatesToRead);
    }
  }
}
