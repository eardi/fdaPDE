rFEM: ../src/rFEM.cpp
	R CMD build ../

all: clean rFEM
	R CMD INSTALL --no-lock rFEM_0.9-1.tar.gz
	
clean:
	rm -f *.tar.gz *.o
	rm -f ../src/*.o ../src/*.so 
	rm -f ../Debug/src/*.o ../Debug/src/*.so 