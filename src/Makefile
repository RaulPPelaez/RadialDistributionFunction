

ARCH=52

#Uncomment if you have boost, it will greatly improve performance
#USE_BOOST=-DUSE_BOOST
all: #main.o
	mkdir -p ../bin
	nvcc -g -I . -arch=sm_$(ARCH) -std=c++11 -O3 main.cu -o ../bin/rdf $(USE_BOOST)

depend: .depend

.depend: main.cu
	rm -f .depend
	g++ -MM -xc++ $^ -MF .depend

include .depend

# %.o:%.cu
# 	nvcc -g  -I . -arch=sm_$(ARCH) -std=c++11 -O3 $(USE_BOOST) $< -c -o $@

