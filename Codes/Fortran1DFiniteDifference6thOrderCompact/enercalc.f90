	SUBROUTINE enercalc(Nx,dt,Es,enkin,enstr,enpot,en,dx,u,uold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine program calculates the energy for the nonlinear 
	! Klein-Gordon equation in 1 dimension
	! u_{tt}-u_{xx}+u=Es*|u|^2u
	!
	! The energy density is given by 
	! 0.5u_t^2+0.5u_x^2+0.5u^2+Es*0.25u^4
	!
	! INPUT 
	!
	! .. Scalars ..
	!  Nx				= number of grid points in x
	!  dt				= timestep
	!  Es				= +1 for focusing, -1 for defocusing
	!  dx                           = spatial discretization
        ! .. Vectors ..
	!  u 				= approximate solution
	!  uold 			= approximate solution
	!
	! OUTPUT
	!
	! .. Scalars ..
	!  enkin			= Kinetic energy
	!  enstr			= Strain energy
	!  enpot			= Potential energy
	!  en				= Total energy
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
	! 	OpenMP library
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)			:: Nx
	REAL(KIND=8), INTENT(IN)			:: dt,Es,dx
	REAL(KIND=8), DIMENSION(1:Nx),INTENT(IN)	:: u,uold
	REAL(KIND=8), INTENT(OUT)			:: enkin,enstr
	REAL(KIND=8), INTENT(OUT)			:: enpot,en	
	INTEGER(KIND=4)					:: i,ip

	
	!.. Strain energy ..
        enstr=0.0d0
	!!$OMP PARALLEL DO PRIVATE(i,ip) SCHEDULE(static)
	DO i=1,Nx
                ip=1+mod(i+1+Nx-1,Nx)
		enstr=enstr+0.5d0*abs(  0.5d0*(u(ip)+uold(ip)-u(i)-uold(i))/dx  )**2
	END DO
	!!$OMP END PARALLEL DO

	! .. Kinetic Energy ..
        enkin=0.0d0
	!!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		enkin=enkin+0.5d0*( abs(u(i)-uold(i))/dt )**2
	END DO
	!!$OMP END PARALLEL DO
	
	! .. Potential Energy ..
        enpot=0.0d0
	!!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		enpot=enpot+0.5d0*(abs((u(i)+uold(i))*0.50d0))**2 &
			   -0.125d0*Es*(abs(u(i))**4+abs(uold(i))**4)
	END DO
	!!$OMP END PARALLEL DO
        enpot=enpot/Nx
        enkin=enkin/Nx
        enstr=enstr/Nx	
	en=enpot+enkin+enstr
	
	END SUBROUTINE enercalc
