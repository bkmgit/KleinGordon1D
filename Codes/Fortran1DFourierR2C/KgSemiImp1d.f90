	!------------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves nonlinear Klein-Gordon equation in 1 dimension
	! u_{tt}-u_{xx}+u=Es*|u|^2u
	! using a second order implicit-explicit time stepping scheme.
	!
	! The boundary conditions are u(x=0)=u(2*Lx*\pi), 
	! The initial condition is u=sqrt(2)*sech((x-c*t)/sqrt(1-c^2));
	!
	! .. Parameters ..
	!  Nx                   = number of modes in x - power of 2 for FFT
	!  Nt			= number of timesteps to take
	!  Tmax			= maximum simulation time
	!  plotgap		= number of timesteps between plots
	!  FFTW_IN_PLACE 	= value for FFTW input 
	!  FFTW_MEASURE 	= value for FFTW input
	!  FFTW_EXHAUSTIVE 	= value for FFTW input
	!  FFTW_PATIENT 	= value for FFTW input    
	!  FFTW_ESTIMATE 	= value for FFTW input
	!  FFTW_FORWARD         = value for FFTW input
	!  FFTW_BACKWARD	= value for FFTW input	
	!  pi = 3.14159265358979323846264338327950288419716939937510d0
	!  Lx			= width of box in x direction
	!  ES			= +1 for focusing and -1 for defocusing
	! .. Scalars ..
	!  i			= loop counter in x direction
	!  n			= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  start		= variable to record start time of program
	!  finish		= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfx		= Forward 2d fft plan 
	!  planbx		= Backward 2d fft plan
	!  dt			= timestep
	!  ierr			= error code
	!  plotnum		= number of plot
	! .. Arrays ..
	!  unew 		= approximate solution
	!  vnew 		= Fourier transform of approximate solution
	!  u 			= approximate solution
	!  v 			= Fourier transform of approximate solution
	!  uold 		= approximate solution
	!  vold 		= Fourier transform of approximate solution
	!  nonlin 		= nonlinear term, u^3
	!  nonlinhat 		= Fourier transform of nonlinear term, u^3
	! .. Vectors ..
	!  kx			= fourier frequencies in x direction
	!  x			= x locations
	!  time			= times at which save data
	!  en			= total energy	
	!  enstr		= strain energy
	!  enpot		= potential energy
	!  enkin		= kinetic energy
	!  name_config		= array to store filename for data to be saved    		
	!  fftfx		= array to setup 1D Fourier transform
	!  fftbx		= array to setup 1D Fourier transform
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
	!	getgrid.f90	-- Get initial grid of points
	! 	initialdata.f90 -- Get initial data
	!	enercalc.f90 -- Subroutine to calculate the energy
	!	savedata.f90 -- Save initial data
	!	storeold.f90 -- Store old data
	! External libraries required
	! 	FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)
	! 	MPI library
		
	PROGRAM Kg
	USE mpi	 	   
	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4), PARAMETER		:: Nx=512
	INTEGER(kind=4), PARAMETER		:: Nt=128
	INTEGER(kind=4), PARAMETER		:: plotgap=32	
	REAL(kind=8), PARAMETER			:: &
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER			::  Lx=9.0d0
	REAL(kind=8), PARAMETER			::  Es=1.0d0	
        REAL(kind=8), PARAMETER         	::  c=0.5
	REAL(kind=8)				::  dt=0.50d0/REAL(Nt,kind(0d0))     
	COMPLEX(kind=8), DIMENSION(1:1+Nx/2)	::  kx 
	REAL(kind=8),  	 DIMENSION(1:Nx)	::  x
	REAL(kind=8), DIMENSION(1:Nx)           ::  u,nonlin
	COMPLEX(kind=8), DIMENSION(1:1+Nx/2)    ::  v,nonlinhat
	REAL(kind=8), DIMENSION(1:Nx)           ::  uold
	COMPLEX(kind=8), DIMENSION(1:1+Nx/2)    ::  vold
	REAL(kind=8), DIMENSION(1:Nx)           ::  unew
	COMPLEX(kind=8), DIMENSION(1:1+Nx/2)    ::  vnew
	REAL(kind=8), DIMENSION(1:Nt+1)         ::  time,enkin,enstr,enpot,en
        REAL(kind=8)                            ::  runtime,start,finish,err
	INTEGER(kind=4)				::  ierr,i,n,nn
	INTEGER(kind=4)				::  plotnum
	INTEGER(kind=4), PARAMETER		::  ntable=2*Nx+64, nwork=4*Nx
	INTEGER(kind=4), PARAMETER		::  isign=0, isignf=-1, isignb=1
	REAL(kind=8),PARAMETER                  ::  scale=1.0d0,scalef=1.0d0,scaleb=1.0d0/DBLE(Nx)
        REAL(kind=8), DIMENSION(1:ntable), SAVE ::  tabledz,tablezd
        REAL(kind=8), DIMENSION(1:nwork)        ::  work
        INTEGER(kind=4), SAVE                   ::  isys
        CHARACTER*100				::  name_config
	
        CALL MPI_INIT(ierr)
	! set up fft
        CALL DZFFT(isign,Nx,scale,u,v,tabledz,work,isys)
        CALL ZDFFT(isign,Nx,scale,v,u,tablezd,work,isys)

	PRINT *,'Setup FFTs'
	! setup fourier frequencies
	CALL getgrid(Nx,Lx,pi,name_config,x,kx)
	PRINT *,'Setup grid and fourier frequencies'
	CALL initialdata(Nx,dt,c,x,u,uold)
	plotnum=1	
	name_config = 'data/u' 
	!CALL savedata(Nx,plotnum,name_config,u)
        CALL DZFFT(isignf,Nx,scalef,u,   v,   tabledz,work,isys)
        CALL DZFFT(isignf,Nx,scalef,uold,vold,tabledz,work,isys)
        PRINT *,'Transformed initialdata'
	CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
			enkin(plotnum),enstr(plotnum),&
			enpot(plotnum),en(plotnum),&
			kx,nonlin,nonlinhat,&
			v,vold,u,uold,tabledz,tablezd,work)

			
	PRINT *,'Got initial data, starting timestepping'
	time(plotnum)=0.0d0
        start = MPI_WTIME()
        nn=1
	DO n=1,(Nt/3)					
                ! Iteration 1
		DO i=1,Nx
			nonlin(i)=u(i)**3
		END DO
                CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,tabledz,work,isys)
		DO i=1,Nx/2 + 1
			vnew(i)=( 0.25d0*(kx(i)*kx(i) - 1.0d0) &
				*(2.0d0*v(i)+vold(i)) &
				+(2.0d0*v(i)-vold(i))/(dt*dt) &
				+Es*nonlinhat(i) ) &
				/(1.0d0/(dt*dt) - 0.25d0*(kx(i)*kx(i) - 1.0d0))
		END DO
                CALL ZDFFT(isignb,Nx,scaleb,vnew,unew,tablezd,work,isys)
                nn=nn+1
		IF (mod(nn,plotgap)==0) then
			plotnum=plotnum+1
			time(plotnum)=n*dt
			PRINT *,'time',n*dt
			CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
				enkin(plotnum),enstr(plotnum),&
				enpot(plotnum),en(plotnum),kx,&
				nonlin,nonlinhat,vnew,v,unew,u,&
                                tabledz,tablezd,work)
		END IF
                
                ! Iteration 2
                DO i=1,Nx
                        nonlin(i)=unew(i)**3
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        vold(i)=( 0.25d0*(kx(i)*kx(i) - 1.0d0) &
                                *(2.0d0*vnew(i)+v(i)) &
                                +(2.0d0*vnew(i)-v(i))/(dt*dt) &
                                +Es*nonlinhat(i) ) &
                                /(1.0d0/(dt*dt) - 0.25d0*(kx(i)*kx(i)-1.0d0))
                END DO
                call ZDFFT(isignb,Nx,scaleb,vold,uold,tablezd,work,isys)
                nn=nn+1
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=n*dt
                        PRINT *,'time',n*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                enkin(plotnum),enstr(plotnum),&
                                enpot(plotnum),en(plotnum),kx,&
                                nonlin,nonlinhat,vold,vnew,uold,unew,&
                                tabledz,tablezd,work)
                END IF
                
                ! Iteration 3
                DO i=1,Nx
                        nonlin(i)=uold(i)**3
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        v(i)=( 0.25d0*(kx(i)*kx(i) - 1.0d0) &
                                *(2.0d0*vold(i)+vnew(i)) &
                                +(2.0d0*vold(i)-vnew(i))/(dt*dt) &
                                +Es*nonlinhat(i) ) &
                                /(1.0d0/(dt*dt) - 0.25d0*(kx(i)*kx(i)-1.0d0))
                END DO
                CALL ZDFFT(isignb,Nx,scaleb,v,u,tablezd,work,isys)
                nn=nn+1
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=n*dt
                        PRINT *,'time',n*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                enkin(plotnum),enstr(plotnum),&
                                enpot(plotnum),en(plotnum),kx,&
                                nonlin,nonlinhat,v,vold,u,uold,&
                                tabledz,tablezd,work)
                END IF
	END DO	
        ! clean up steps
        IF (mod(Nt,3).GT.0) THEN
                DO i=1,Nx
                        nonlin(i)=u(i)**3
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        vnew(i)=( 0.25d0*(kx(i)*kx(i) - 1.0d0) &
                                *(2.0d0*v(i)+vold(i)) &
                                +(2.0d0*v(i)-vold(i))/(dt*dt) &
                                +Es*nonlinhat(i) ) &
                                /(1.0d0/(dt*dt) - 0.25d0*(kx(i)*kx(i)-1.0d0))
                END DO
                CALL ZDFFT(isignb,Nx,scaleb,vnew,unew,tablezd,work,isys)
                nn=nn+1
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=n*dt
                        PRINT *,'time',n*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                enkin(plotnum),enstr(plotnum),&
                                enpot(plotnum),en(plotnum),kx,&
                                nonlin,nonlinhat,vnew,v,unew,u,&
                                tabledz,tablezd,work)
                END IF
                
        END IF
        IF (mod(Nt,3).GT.1) THEN
                DO i=1,Nx
                        nonlin(i)=unew(i)**3
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,tabledz,work,isys)
                DO i=1,Nx/2 +1
                        vold(i)=( 0.25d0*(kx(i)*kx(i) - 1.0d0) &
                                *(2.0d0*vnew(i)+v(i)) &
                                +(2.0d0*vnew(i)-v(i))/(dt*dt) &
                                +Es*nonlinhat(i) ) &
                                /(1.0d0/(dt*dt) -0.25d0*(kx(i)*kx(i)-1.0d0))
                END DO
                CALL ZDFFT(isignb,Nx,scaleb,vold,uold,tablezd,work,isys)
                nn=nn+1
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=n*dt
                        PRINT *,'time',n*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                      enkin(plotnum),enstr(plotnum),&
                                      enpot(plotnum),en(plotnum),kx,&
                                      nonlin,nonlinhat,vold,vnew,uold,unew,&
                                      tabledz,tablezd,work)
                END IF
        END IF	
	PRINT *,'Finished time stepping'
        finish = MPI_WTIME()
        runtime=REAL(finish-start,kind(0d0))
	PRINT*,'Program took ',runtime,'for Time stepping'
	CALL saveresults(Nt,plotgap,time(1:1+n/plotgap),en(1:1+n/plotgap),&
			enstr(1:1+n/plotgap),enkin(1:1+n/plotgap),&
                        enpot(1:1+n/plotgap))	
			
	! Save times at which output was made in text format
	PRINT *,'Saved data'
        IF(mod(Nt,3)==0) CALL errorcalc(Nt,Nx,dt,c,x,unew,err) 
        IF(mod(Nt,3)==1) CALL errorcalc(Nt,Nx,dt,c,x,uold,err)
        IF(mod(Nt,3)==2) CALL errorcalc(Nt,Nx,dt,c,x,u,err)
        name_config = 'data/errortime'
        CALL savetimeanderror(Nx,Nt,err,runtime,name_config)

	PRINT *,'Program execution complete'
        CALL MPI_FINALIZE(ierr)
	END PROGRAM Kg
