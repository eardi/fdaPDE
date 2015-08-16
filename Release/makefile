FEMr: ../src/FEMr.cpp
	R CMD build ../

all: clean FEMr
	R CMD INSTALL --no-lock FEMr_0.9-1.tar.gz
	
clean:
	rm -f *.tar.gz *.o
	rm -f ../src/*.o ../src/*.so 
	rm -f ../Debug/src/*.o ../Debug/src/*.so 