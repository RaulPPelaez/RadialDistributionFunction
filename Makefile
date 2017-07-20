

all:
	$(MAKE) -C src

install: all
	mkdir -p ~/bin
	mv bin/rdf ~/bin/

test: all
	$(MAKE) -C test
clean:
	rm -f bin/rdf src/main.o
