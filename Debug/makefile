## This is a dummy makefile, for package creation
# use R CMD build "directory"

rFEM.o: src/main.cpp
	g++ -Wall -pedantic -std=c++0x -c -I/usr/share/R/include -DNDEBUG -fpic  -O3 -pipe  -g  -c src/main.cpp -o rFEM.o

all: rFEM.o
	g++ -g -o rFEM rFEM.o