

CXX = g++
LAPACK = /Users/mmccull/lib
OPTS = -O3 -ftree-vectorize

gb: gb.cpp psflib.h psflib.cpp dcdlib.h dcdlib.c stringlib.h stringlib.cpp pdblib.h pdblib.c
	$(CXX) -c gb.cpp stringlib.cpp pdblib.c psflib.cpp dcdlib.c $(OPTS) 
	$(CXX) gb.o stringlib.o pdblib.o psflib.o dcdlib.o $(OPTS) -L$(LAPACK) -lclapack -lcblas -lf2c -lm -o gb.x

