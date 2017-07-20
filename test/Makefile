


all:
	octave -qf random_positions.m > inipos 2> /dev/null

	bash test.bash ../bin/rdf  "-device CPU"
	mkdir -p CPU
	mv *gdr CPU

	bash test.bash ../bin/rdf  "-device GPU"
	mkdir -p GPU
	mv *gdr GPU
	rm -f inipos octave-workspace

	bash checkResults.bash

clean:
	rm -rf CPU GPU 
