#define the complier
COMPILER = gfortran
# compilation settings, optimization, precision, parallelization
        FLAGS = -O3 -fopenmp

# source list for main program
SOURCES =  KgSemiImp1d.f90 initialdata.f90 savedata.f90 getgrid.f90 \
           errorcalc.f90 storeold.f90 saveresults.f90 enercalc.f90
                   
test: $(SOURCES)
                ${COMPILER} -o kg $(FLAGS) $(SOURCES)    

clean:
        rm kg *.o 
