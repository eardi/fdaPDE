fdaPDE: ../src/fdaPDE.cpp
	R CMD build ../

all: clean fdaPDE
	R CMD INSTALL --no-lock fdaPDE_0.1-2.tar.gz
	
clean:
	rm -f *.tar.gz *.o
	rm -f ../src/*.o ../src/*.so 
	rm -f ../Debug/src/*.o ../Debug/src/*.so 
