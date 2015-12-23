## This is a dummy makefile, for package creation
# use R CMD build "directory"

main.o: src/main.cpp
	g++ -Wall -pedantic -std=c++0x -c -DNDEBUG -fpic -pipe  -g  -c src/main.cpp -o main.o
integration.o: ../src/integration.cpp
	g++ -Wall -pedantic -std=c++0x -c -DNDEBUG -fpic -pipe  -g  -c ../src/integration.cpp -o integration.o
mesh_objects.o: ../src/mesh_objects.cpp
	g++ -Wall -pedantic -std=c++0x -c -DNDEBUG -fpic -pipe  -g  -c ../src/mesh_objects.cpp -o mesh_objects.o

all: main.o integration.o mesh_objects.o
	g++ -o fdaPDE main.o integration.o mesh_objects.o