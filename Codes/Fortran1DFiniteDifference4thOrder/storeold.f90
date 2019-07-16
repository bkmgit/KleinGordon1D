	SUBROUTINE storeold(Nx,unew,u,uold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine copies arrays for a
	! pseudospectral simulation of the 1D nonlinear Klein-Gordon equation
	!
	! u_{tt}-u_{xx}+u=Es*u^3
	!
	! INPUT
	!
	! .. Parameters ..
	!  Nx				= number of x grid points
	!  .. Arrays ..
	!  unew 			= approximate solution
	!  u 				= approximate solution
	!  uold 			= approximate solution
	!
	! OUTPUT
	!
	!  u 				= approximate solution
	!  uold 			= approximate solution
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
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! OpenMP library
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)			:: Nx
	REAL(KIND=8), DIMENSION(1:NX), INTENT(OUT)	:: uold
	REAL(KIND=8), DIMENSION(1:NX), INTENT(INOUT)    :: u
	REAL(KIND=8), DIMENSION(1:NX), INTENT(IN)	:: unew
	INTEGER(kind=4)					:: i

	!$OMP PARALLEL PRIVATE(i)
	!$OMP DO SCHEDULE(static)
	DO i=1,Nx
		uold(i)=u(i)
	END DO
	!$OMP END DO 
	!$OMP DO SCHEDULE(static)
	DO i=1,Nx
		u(i)=unew(i)
	END DO
	!$OMP END DO 
	!$OMP END PARALLEL 
	
	END SUBROUTINE storeold
