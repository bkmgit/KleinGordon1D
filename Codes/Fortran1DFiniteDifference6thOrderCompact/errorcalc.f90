	SUBROUTINE errorcalc(Nt,Nx,dt,c,x,u,err)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine calculates the error in the approximate solution
        ! of the nonlinear Klein-Gordon equation in 1 dimensions
	! u_{tt}-u_{xx}+u=Es*u^3+
	!
	! The boundary conditions are u(x=-Lx*\pi)=u(x=Lx*\pi), 
	! The initial condition is u=sqrt(2)*sech((x(i)-c*t)/sqrt(1-c**2)
	!
	! INPUT
	!
	! .. Parameters ..
	!  Nx	= number of modes in x - power of 2 for FFT
        !  Nt   = number of timesteps taken
        !  dt   = timestep
        !  c    = wavespeed
	! .. Vectors ..
	!  x	= x locations
	!  u    = solution at current timestep
        !
	! OUTPUT
	!
	! .. L2 norm of error ..
	! err = L2 norm of error
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
	INTEGER(KIND=4), INTENT(IN)				:: Nx,Nt
	REAL(KIND=8), DIMENSION(1:NX), INTENT(IN) 		:: x
	REAL(KIND=8), DIMENSION(1:NX), INTENT(IN)		:: u
	INTEGER(kind=4)						:: i
        REAL(KIND=8), INTENT(IN)                                :: c,dt
        REAL(KIND=8)                                            :: t
        REAL(KIND=8), INTENT(OUT)                               :: err
        t=REAL(Nt,kind=8)*dt
        err=0.0d0
	!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:err) SCHEDULE(static)
	DO i=1,Nx
		err=err+&
                (u(i)-sqrt(2.0d0)/(cosh((x(i)-c*t)/sqrt(1.0d0-c**2))) )**2
        !        err=err+&
        !        (u(i)-sqrt(2.0d0)/cosh(x(i)+t))**2
	END DO
	!$OMP END PARALLEL DO
        PRINT *,'The final L2 error is ',sqrt(err/REAL(Nx,kind=8))
	END SUBROUTINE errorcalc
