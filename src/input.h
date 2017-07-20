/*Raul P. Pelaez 2017. Contains utility classes/functions to read snapshots from a file

 */
#ifndef INPUT_H
#define INPUT_H

#include"vector_algebra.cuh"
#include<fstream>
#include<sstream>
#include<string>
namespace gdr{
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
    template<class Iterator>
    void parseNextLine(Iterator numbersInLine, int numberColumnsToRead){
      this->goToNextLine();
      //This could be faster
      stringstream ss;
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
