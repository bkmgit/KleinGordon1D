	SUBROUTINE initialdata(Nx,dt,c,x,u,uold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine gets initial data for nonlinear Klein-Gordon equation
	! in 1 dimensions
	! u_{tt}-u_{xx}+u=Es*u^3
	!
	! The boundary conditions are u(x=-Lx*\pi)=u(x=Lx*\pi), 
	! The initial condition is u=sqrt(2)*sech((x(i)-c*t)/sqrt(1-c**2)
	!
	! INPUT
	!
	! .. Parameters ..
	!  Nx	= number of modes in x - power of 2 for FFT
        !  dt   = timestep
        !  c    = wavespeed
	! .. Vectors ..
	!  x	= x locations
	!
	! OUTPUT
	!
	! .. Arrays ..
	!  u 	= initial solution
	!  uold = approximate solution at previous timestep
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i	= loop counter in x direction
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
	INTEGER(KIND=4), INTENT(IN)				:: Nx
	REAL(KIND=8), DIMENSION(1:NX), INTENT(IN) 		:: x
	REAL(KIND=8), DIMENSION(1:NX), INTENT(OUT)		:: u,uold
	INTEGER(kind=4)						:: i
        REAL(KIND=8), INTENT(IN)                                :: c,dt
        REAL(KIND=8)                                            :: t=0.0d0
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		u(i)=sqrt(2.0d0)/(cosh((x(i)-c*t)/sqrt(1.0d0-c**2)))
                !u(i)=sqrt(2.0d0)/cosh(x(i)+t)
	END DO
	!$OMP END PARALLEL DO
        t=-1.0d0*dt;
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		uold(i)=sqrt(2.0d0)/(cosh((x(i)-c*t)/sqrt(1.0d0-c**2)))
                !uold(i)=sqrt(2.0d0)/cosh(x(i)+t)
	END DO
	!$OMP END PARALLEL DO
	
	END SUBROUTINE initialdata
