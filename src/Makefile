

ARCH=30

all: main.o
	mkdir -p ../bin
	nvcc -g -I . -arch=sm_$(ARCH) -std=c++11 -O3 $^ -o ../bin/rdf

depend: .depend

.depend: main.cu
	rm -f .depend
	g++ -MM -xc++ $^ -MF .depend

include .depend

%.o:%.cu
	nvcc -g  -I . -arch=sm_$(ARCH) -std=c++11 -O3 $< -c -o $@

