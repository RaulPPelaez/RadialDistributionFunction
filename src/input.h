/*Raul P. Pelaez 2017. Contains utility classes/functions to read snapshots from a file

 */
#ifndef INPUT_H
#define INPUT_H

#include"vector_algebra.cuh"
#include<fstream>
#include<sstream>
#include<string>
#include<memory>
#include<iostream>

namespace gdr{
    //Reads numbers from a file or standard input line by line
  class InputParse{
    std::shared_ptr<std::istream> input;
    std::string currentLine;
  public:
    //Read from stdin by default
    InputParse(){ }

    //Take input from cin
    bool open(){
      //The lambda ensures cin is not deleted
      input.reset(&std::cin, [](...){});
    
      if(!input->good()){
	std::cerr<<"ERROR: Unable to read from stdin!"<<std::endl;
	return false;
      }
      return true;
    }
    //Open and read from a file
    bool open(std::string fileName){
      input.reset(new std::ifstream(fileName.c_str()));
      if(!input->good()){
	std::cerr<<"ERROR: Unable to open file!"<<std::endl;
	return false;
      }
      return true;
    }

    //Reads a line from input
    std::string goToNextLine(){
      if(!input) std::cerr<<"ERROR: No open file!"<<std::endl;
      getline(*input, currentLine);
      return currentLine;
    }
    //Reads a line and
    template<class Iterator>
    void parseNextLine(Iterator numbersInLine, int numberColumnsToRead){
      this->goToNextLine();
      //This could be faster
      std::stringstream ss;
      ss.str(currentLine);
      for(int i=0; i<numberColumnsToRead; i++){
	ss>>numbersInLine[i];
      }

    }

  };


}
#endif

namespace gdr{
  //Reads an entire frame to the second argument. Meant to be used with real4,real3 or real2
  template<class vecType>
  void readFrame(InputParse &input, vecType *pos, int numberParticles, int nCoordinatesToRead){

    //Ignore frame separator
    input.goToNextLine();

    for(int i= 0; i<numberParticles; i++){

      input.parseNextLine((real*)&pos[i], nCoordinatesToRead);

    }

  }
}
