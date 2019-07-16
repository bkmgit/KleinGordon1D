#define the complier
COMPILER = gfortran
# compilation settings, optimization, precision, parallelization
	FLAGS = -O3 -fopenmp

# libraries
LIBS = -L/usr/lib64 -lfftw3 -lfftw3_threads -lm 
# source list for main program
SOURCES =  KgSemiImp1d.f90 initialdata.f90 savedata.f90 getgrid.f90 \
	   errorcalc.f90 storeold.f90 saveresults.f90 enercalc.f90
		   
test: $(SOURCES)
		${COMPILER} -o kg $(FLAGS) $(SOURCES) $(LIBS)   
		
clean:
	rm kg *.dat *.o 
