

EXE=build/bin/rdf
all:
	mkdir -p build && cd build && cmake ../..  && make; cd -
	seq 1 1e8 | awk '{print rand()*10, rand()*10, rand()*10}' > inipos
	echo "CPU: ">>timing.dat
	bash test.bash $(EXE)  "-device CPU"
	mkdir -p CPU
	mv *gdr CPU

	echo "GPU: ">>timing.dat
	bash test.bash $(EXE)  "-device GPU"	
	mkdir -p GPU
	mv *gdr GPU
	rm -f inipos octave-workspace

	bash checkResults.bash

clean:
	rm -rf CPU GPU timing.dat
