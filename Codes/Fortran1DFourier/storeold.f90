	SUBROUTINE storeold(Nx,unew,u,uold,vnew,v,vold)
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
	!  Nx				= number of modes in x - power of 2 for FFT
	!  .. Arrays ..
	!  unew 			= approximate solution
	!  vnew 			= Fourier transform of approximate solution
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
	!
	! OUTPUT
	!
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
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
	!USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)			:: Nx
	COMPLEX(KIND=8), DIMENSION(1:NX), INTENT(OUT)	:: vold,uold
	COMPLEX(KIND=8), DIMENSION(1:NX), INTENT(INOUT) :: u,v
	COMPLEX(KIND=8), DIMENSION(1:NX), INTENT(IN)	:: unew,vnew
	INTEGER(kind=4)					:: i

	!!$OMP PARALLEL PRIVATE(i)
	
	!!$OMP DO SCHEDULE(static)
	DO i=1,Nx
		vold(i)=v(i)
	END DO
	!!$OMP END DO NOWAIT

	!!$OMP DO SCHEDULE(static)
	DO i=1,Nx
		uold(i)=u(i)
	END DO
	!!$OMP END DO NOWAIT
        !!$OMP BARRIER	
	!!$OMP DO SCHEDULE(static)
	DO i=1,Nx
		u(i)=unew(i)
	END DO
	!!$OMP END DO NOWAIT
	
	!!$OMP DO SCHEDULE(static)
	DO i=1,Nx
		v(i)=vnew(i)
	END DO
	!!$OMP END DO NOWAIT
	
	!!$OMP END PARALLEL 
	
	END SUBROUTINE storeold
