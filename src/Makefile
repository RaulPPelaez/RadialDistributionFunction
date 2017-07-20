

all: main.o
	mkdir -p ../bin
	nvcc  -arch=sm_52 -std=c++11 -O3 $^ -o ../bin/rdf

depend: .depend

.depend: main.cu
	rm -f .depend
	g++ -MM -xc++ $^ -MF .depend

include .depend

%.o:%.cu
	nvcc  -arch=sm_52 -std=c++11 -O3 $< -c -o $@

