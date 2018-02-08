

all:
	$(MAKE) -C src

install: all
	mkdir -p ~/bin
	mv bin/rdf ~/bin/

test: all
	$(MAKE) -C test

clean:
	$(MAKE) -C test clean
	$(MAKE) -C src clean
	rm -f bin/rdf src/main.o

