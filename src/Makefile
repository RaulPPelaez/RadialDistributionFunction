

ARCH=35

#Uncomment if you have boost, it will greatly improve performance
#USE_BOOST=-DUSE_BOOST
all: gitversion.h #main.o
	mkdir -p ../bin
	nvcc -g -I . -arch=sm_$(ARCH) -std=c++11 -O3 main.cu -o ../bin/rdf $(USE_BOOST)

gitversion.h: ../.git/HEAD ../.git/index
	echo "#ifndef GITVERSION" > $@
	echo "#define GITVERSION \"$(shell git rev-parse HEAD)\"" >> $@
	echo "#endif" >> $@

clean:
	rm -f gitversion.h
# %.o:%.cu
# 	nvcc -g  -I . -arch=sm_$(ARCH) -std=c++11 -O3 $(USE_BOOST) $< -c -o $@

