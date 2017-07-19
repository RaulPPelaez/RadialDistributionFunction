

all:
	nvcc  -arch=sm_52 -std=c++11 -O3 rdf.cu -o rdf
clean:
	rm -f rdf
