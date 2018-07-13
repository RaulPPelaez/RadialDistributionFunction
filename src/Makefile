#Uncomment if you have boost, it will greatly improve performance
#USE_BOOST=-DUSE_BOOST

#The target CUDA compute capability(s), if not set it will be autodetected
#ARCH ?= 30 35 50 52 60 61
ARCH=

#Uncomment to compile in double precision mode, single by default
#DOUBLE_PRECISION=-DDOUBLE_PRECISION

#C++ compiler, I tested up to clang++-5.0
#CXX=clang++-5.0
CXX=g++

#Cuda version (assumed to be in /usr/local/cuda*) You can change this in CUDA_ROOT
#If not set it will be autodetected
CUDA_VER=

ifeq ($(CUDA_VER),)
CUDA_VER:=$(shell ls -d /usr/local/cuda*/ | grep -Eo '\-[0-9]\.[0-9]' | cut -d- -f2 | sort -grk1 | head -1)
endif
CUDA_ROOT=/usr/local/cuda-$(CUDA_VER)

#Flags to $(CXX)
CPU= -O3 -funroll-loops  -fno-math-errno -fno-signed-zeros  -march=native -fPIC -Wall -Wextra -Wno-unused-parameter -Wno-unused-function -Wno-sign-compare


#If arch was not set, autodetect all GPUs in the system
ifeq ($(ARCH),)
GENCODE_FLAGS:=$(shell printf '\#include<cstdio>\n int main(){int nD;cudaGetDeviceCount(&nD);for(int i=0;i<nD;i++){cudaDeviceProp dp;cudaGetDeviceProperties(&dp, i);std::printf("%%d\\n", dp.major*10+dp.minor);} return 0;}' | $(CUDA_ROOT)/bin/nvcc -Wno-deprecated-gpu-targets -x cu - -o /tmp/listarch --run | sort -g -k1 | uniq | awk 'END{system("rm -f /tmp/listarch")}{print "-gencode arch=compute_"$$1",code=sm_"$$1}')
else
$(foreach sm,$(ARCH),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
endif


INCLUDEFLAGS= -I$(CUDA_ROOT)/include
BASIC_LINE= $(CUDA_ROOT)/bin/nvcc   -O3 -std=c++11 -x cu    $(INCLUDEFLAGS) $(OPTIONS)  -ccbin="$(CXX)" -Xcompiler="$(CPU)" $(GENCODE_FLAGS) -L$(CUDA_ROOT)/lib64



all: gitversion.h #main.o
	mkdir -p ../bin
	$(BASIC_LINE) -I . main.cu -o ../bin/rdf $(USE_BOOST)

gitversion.h: ../.git/HEAD ../.git/index
	echo "#ifndef GITVERSION" > $@
	echo "#define GITVERSION \"$(shell git rev-parse HEAD)\"" >> $@
	echo "#endif" >> $@

clean:
	rm -f gitversion.h

# %.o:%.cu
# 	nvcc -g  -I . -arch=sm_$(ARCH) -std=c++11 -O3 $(USE_BOOST) $< -c -o $@

