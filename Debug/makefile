## This is a dummy makefile, for package creation
# use R CMD build "directory"

FEMr.o: src/main.cpp
	g++ -Wall -pedantic -std=c++0x -c -I/usr/share/R/include -DNDEBUG -fpic  -O3 -pipe  -g  -c src/main.cpp -o FEMr.o

all: FEMr.o
	g++ -g -o FEMr FEMr.o