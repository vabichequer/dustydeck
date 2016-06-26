# Makefile to build dustyFortran
#
# Vicenzo Abichequer, Computer Engineering
# Mercer University, 2016
# 

FC = gfortran
CC = g++
FLAGS = -cpp -funroll-loops -DFOUR -DONE -DTWO -O2 -march=nocona -DSTRIDEACCESS
ALLDOUBLEPRECISION = -DNOPRECISIONLOSS
VECTORIZER = -ftree-vectorize -fdump-tree-vect -ftree-vectorizer-verbose=5
GCOV = -fprofile-arcs -ftest-coverage
NOTUSED = -ffast-math 

all: mainf mainc dusty

evropt: cpp fortran dustyAll

dusty: 
	$(FC) -O2 -c dusty.f
	$(FC) -O2 dusty.f -o dusty
	
dustyAll:	
	$(FC) -c dusty.f
	$(FC) dusty.f -o dusty1
	$(FC) -O2 -c dusty.f
	$(FC) -O2 dusty.f -o dusty2

baseline:
	$(FC) -cpp  -O1 -c mainf.f90
	$(FC) -cpp  -O1 mainf.f90 -o mainf1
	$(FC) -cpp  -O3 -c mainf.f90
	$(FC) -cpp  -O3 mainf.f90 -o mainf2
	$(CC) -cpp -O1 -c mainc.cpp
	$(CC) -cpp -O1 mainc.cpp -o mainc1
	$(CC) -cpp -O3 -c mainc.cpp
	$(CC) -cpp -O3 mainc.cpp -o mainc2
	$(FC) -O1 -c dusty.f
	$(FC) -O1 dusty.f -o dusty1
	$(FC) -O3 -c dusty.f
	$(FC) -O3 dusty.f -o dusty2

ebase:
	@echo "************ C++ START ************"
	@echo "*** -O1 ***"
	mainc1
	@echo "*** -O3 ***"
	mainc2
	@echo "************ C++ END ************"
	@echo "************ FORTAN START ************"
	@echo "*** -O1 ***"
	mainf1
	@echo "*** -O3 ***"
	mainf2
	@echo "************ FORTAN END ************"
	@echo "************ DUSTY START ************"
	@echo "*** -O1 ***"
	dusty1
	@echo "*** -O3 ***"
	dusty2
	@echo "************ DUSTY END ************"
	
fortran:
	$(FC) -cpp  -c mainf.f90
	$(FC) -cpp  mainf.f90 -o mainf1
	$(FC) -cpp  -O2 -c mainf.f90
	$(FC) -cpp  -O2 mainf.f90 -o mainf2
	$(FC) -cpp -O2 -funroll-loops -c mainf.f90
	$(FC) -cpp -O2 -funroll-loops mainf.f90 -o mainf3
	$(FC) -cpp -O2 -DSTRIDEACCESS -c mainf.f90
	$(FC) -cpp -O2 -DSTRIDEACCESS mainf.f90 -o mainf4
	$(FC) -cpp -O2 -DONE -c mainf.f90
	$(FC) -cpp -O2 -DONE mainf.f90 -o mainf5
	$(FC) -cpp -O2 -DTWO -c mainf.f90
	$(FC) -cpp -O2 -DTWO mainf.f90 -o mainf6
	$(FC) -cpp -O2 -DFOUR -c mainf.f90
	$(FC) -cpp -O2 -DFOUR mainf.f90 -o mainf7
	$(FC) -cpp -O2 -march=nocona -c mainf.f90
	$(FC) -cpp -O2 -march=nocona mainf.f90 -o mainf8
	$(FC) -cpp -O2 $(VECTORIZER) -c mainf.f90
	$(FC) -cpp -O2 $(VECTORIZER) mainf.f90 -o mainf9
cpp:	
	$(CC) -cpp -c mainc.cpp
	$(CC) -cpp mainc.cpp -o mainc1
	$(CC) -cpp -O2 -c mainc.cpp
	$(CC) -cpp -O2 mainc.cpp -o mainc2
	$(CC) -cpp -O2 -funroll-loops -c mainc.cpp
	$(CC) -cpp -O2 -funroll-loops mainc.cpp -o mainc3
	$(CC) -cpp -O2 -DSTRIDEACCESS -c mainc.cpp
	$(CC) -cpp -O2 -DSTRIDEACCESS mainc.cpp -o mainc4
	$(CC) -cpp -O2 -DONE -c mainc.cpp
	$(CC) -cpp -O2 -DONE mainc.cpp -o mainc5
	$(CC) -cpp -O2 -DTWO -c mainc.cpp
	$(CC) -cpp -O2 -DTWO mainc.cpp -o mainc6
	$(CC) -cpp -O2 -DFOUR -c mainc.cpp
	$(CC) -cpp -O2 -DFOUR mainc.cpp -o mainc7
	$(CC) -cpp -O2 -march=nocona -c mainc.cpp
	$(CC) -cpp -O2 -march=nocona mainc.cpp -o mainc8
	$(CC) -cpp -O2 $(VECTORIZER) -c mainc.cpp
	$(CC) -cpp -O2 $(VECTORIZER) mainc.cpp -o mainc9

mainc: mainc.o
	$(CC) $(FLAGS) $(VECTORIZER) mainc.cpp -o mainc

mainc.o: mainc.cpp
	$(CC) $(FLAGS) $(VECTORIZER) -c mainc.cpp

mainf: mainf.o
	$(FC) $(FLAGS) $(VECTORIZER) mainf.f90 -o mainf

mainf.o: mainf.f90
	$(FC) $(FLAGS) $(VECTORIZER) -c mainf.f90
	
gcovc:
	$(CC) $(FLAGS) $(GCOV) $(VECTORIZER) -c mainc.cpp
	$(CC) $(FLAGS) $(GCOV) $(VECTORIZER) mainc.cpp -o mainc
	mainc
	gcov mainc.cpp

gcovf:
	$(FC) $(FLAGS) $(GCOV) $(VECTORIZER) -c mainf.f90
	$(FC) $(FLAGS) $(GCOV) $(VECTORIZER) mainf.f90 -o mainf
	mainf
	gcov mainf.f90
	
clean:
	rm *.o
	rm *.a

clvect:
	rm *.vect
	
clgcov:
	rm *.gcov
	rm *.gcda
	rm *.gcno
	
clqsub:
	rm vicenzo.*
	
pristine:
	rm *.o
	rm *.a
	touch *.f90 *.cpp

ctags:
	ctags *.f90 *.cpp

efortran:
	@echo "************ FORTAN START ************"
	@echo "*** Stock ***"
	mainf1
	@echo "*** -O2 ***"
	mainf2
	@echo "*** -funroll-loops ***"
	mainf3
	@echo "*** -DSTRIDEACCESS ***"
	mainf4
	@echo "*** -DONE ***"
	mainf5
	@echo "*** -DTWO ***"
	mainf6
	@echo "*** -DFOUR ***"
	mainf7
	@echo "*** -march=nocona ***"
	mainf8
	@echo "*** Vectorizer ***"
	mainf9
	@echo "************ FORTRAN END ************"
	
ecpp:
	@echo "************ C++ START ************"
	@echo "*** Stock ***"
	mainc1
	@echo "*** -O2 ***"
	mainc2
	@echo "*** -funroll-loops ***"
	mainc3
	@echo "*** -DONE ***"
	mainc5
	@echo "*** -DTWO ***"
	mainc6
	@echo "*** -DFOUR ***"
	mainc7
	@echo "*** -march=nocona ***"
	mainc8
	@echo "*** Vectorizer ***"
	mainc9
	@echo "************ C++ END ************"
	
edusty:
	@echo "************ DUSTY START ************"
	@echo "*** Stock ***"
	dusty1
	@echo "*** -O2 ***"
	dusty2
	@echo "************ DUSTY END ************"
	
execall: ecpp efortran edusty
