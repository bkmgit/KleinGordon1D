	SUBROUTINE enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                            enkin,enstr,enpot,en,kx,temp1,temp2,&
                            v,vold,u,uold,tabledz,tablezd,work)

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
	!  Nx				= number of modes in x - power of 2 for FFT
	!  planfx			= Forward 2d fft plan
	!  planbx			= Backward 2d fft plan
	!  dt				= timestep
	!  Es				= +1 for focusing, -1 for defocusing
	! .. Arrays ..
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
	!  temp1			= array to hold temporary values
	!  temp2			= array to hold temporary values
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
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
	! 	FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)
	! 	OpenMP library
	!USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)			     :: Nx
	REAL(KIND=8), INTENT(IN)			     :: dt,Es
	COMPLEX(KIND=8), DIMENSION(1:Nx/2 + 1), INTENT(IN)   :: kx
	REAL(KIND=8), DIMENSION(1:Nx), INTENT(IN)	     :: u,uold
        COMPLEX(KIND=8), DIMENSION(1:Nx/2 + 1), INTENT(IN)   :: v,vold
	REAL(KIND=8), DIMENSION(1:Nx), INTENT(INOUT)	     :: temp1
        COMPLEX(KIND=8), DIMENSION(1:Nx/2 + 1), INTENT(INOUT):: temp2
	REAL(KIND=8), INTENT(OUT)			     :: enkin,enstr
	REAL(KIND=8), INTENT(OUT)			     :: enpot,en
	INTEGER(kind=4), INTENT(IN)                          :: ntable,nwork
        INTEGER(kind=4), INTENT(INOUT)                       :: isys
        INTEGER(kind=4), INTENT(IN)                          :: isignb, isignf
        REAL(kind=8)                                         :: scale
        REAL(kind=8), DIMENSION(1:ntable), INTENT(INOUT)     :: tablezd,tabledz
        REAL(kind=8), DIMENSION(1:nwork), INTENT(INOUT)      :: work
        INTEGER(KIND=4)			                     :: i
	REAL(kind=8), PARAMETER                   :: c=0.5d0,t=0.0d0
        scale=1.0d0/(REAL(Nx,kind(0d0)))
	!.. Strain energy ..
	DO i=1,Nx/2 + 1
		temp2(i)=0.5d0*kx(i)*(vold(i)+v(i))
	END DO
        CALL ZDFFT(isignb,Nx,scale,temp2,temp1,tablezd,work,isys)
        enstr=0.0d0
	DO i=1,Nx
		enstr=enstr+0.5d0*abs(temp1(i))**2
	END DO
	! .. Kinetic Energy ..
        enkin=0.0d0
	DO i=1,Nx
		enkin=enkin+0.5d0*( abs(u(i)-uold(i))/dt )**2
	END DO
	
	! .. Potential Energy ..
        enpot=0.0D0
	DO i=1,Nx
		enpot=enpot+0.5d0*(abs((u(i)+uold(i))*0.50d0))**2&
				-0.125d0*Es*(abs(u(i))**4+abs(uold(i))**4)
	END DO
        !enpot=enpot/Nx
        !enkin=enkin/Nx
        !enstr=enstr/Nx	
	en=enpot+enkin+enstr
	
	END SUBROUTINE enercalc
