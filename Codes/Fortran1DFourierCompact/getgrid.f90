	SUBROUTINE getgrid(Nx,Lx,pi,name_config,x,kx)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine gets grid points and fourier frequencies for a
	! pseudospectral simulation of the 1D nonlinear Klein-Gordon equation
	!
	! u_{tt}-u_{xx}+u=Es*u^3
	!
	! The boundary conditions are u(x=0)=u(2*Lx*\pi) 
	!
	! INPUT
	!
	! .. Scalars ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  pi				= 3.142....
	!  Lx				= width of box in x direction
	! OUTPUT
	!
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  x				= x locations
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!
	! REFERENCES
	!
	! ACKNOWLEDGEMENTS
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	! Check that the initial iterate is consistent with the 
	! boundary conditions for the domain specified
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! OpenMP library
	USE omp_lib		 	   
        IMPLICIT NONE
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)			:: Nx
	REAL(kind=8), INTENT(IN)			:: Lx,pi
	REAL(KIND=8), DIMENSION(1:NX), INTENT(OUT) 	:: x
	COMPLEX(KIND=8), DIMENSION(1:NX), INTENT(OUT)	:: kx
	CHARACTER*100, INTENT(OUT)			:: name_config
	INTEGER(kind=4)					:: i
		
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0d0,kind(0d0))*REAL(i-1,kind(0d0))/Lx  			
	END DO
	!$OMP END PARALLEL DO
	kx(1+Nx/2)=0.0d0
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
	!$OMP END PARALLEL DO
		
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
	END DO
	!$OMP END PARALLEL DO
	! Save x grid points in text format
	name_config = 'xcoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO i=1,Nx
		WRITE(11,*) x(i)
	END DO
	CLOSE(11)
	
	END SUBROUTINE getgrid
