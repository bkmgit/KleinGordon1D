	!------------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves nonlinear Klein-Gordon equation in 1 dimension
	! u_{tt}-u_{xx}+u=Es*|u|^2u
	! using a fourth order implicit time stepping scheme.
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
	! 	OpenMP library
		
	PROGRAM Kg
	USE mpi		 	   
	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4), PARAMETER		:: Nx=1024
	INTEGER(kind=4), PARAMETER		:: Nt=1024 
	INTEGER(kind=4), PARAMETER		:: plotgap=32	
        INTEGER(kind=4), PARAMETER              :: maxiter=50
	REAL(kind=8), PARAMETER			:: &
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER			:: Lx=9.0d0
	REAL(kind=8), PARAMETER			:: Es=1.0d0	
        REAL(kind=8), PARAMETER         	:: c=0.5
        REAL(kind=8), PARAMETER                 :: tol=0.1d0**12
	REAL(kind=8)				:: dt=5.0d0/REAL(Nt,kind(0d0))	
	COMPLEX(kind=8), DIMENSION(1:Nx/2+1)	:: kx 
	REAL(kind=8),  	 DIMENSION(1:Nx)	:: x
	REAL(kind=8), DIMENSION(1:Nx)           :: u,uold,unew,utemp,nonlin,nonlinfix
	COMPLEX(kind=8), DIMENSION(1:Nx/2+1)    :: v,vold,vnew,vtemp,nonlinhat,nonlinfixhat
	REAL(kind=8), DIMENSION(1:Nt+1)         :: time,enkin,enstr,enpot,en
        INTEGER(kind=4), DIMENSION(1:Nt+1)      :: itercount
        REAL(kind=8)                            :: chg,err,runtime,start,finish
	INTEGER(kind=4)				:: ierr,i,j,n,nn,iter
	INTEGER(kind=4)				:: count_rate,plotnum
	INTEGER(kind=4), PARAMETER              :: ntable=2*Nx+64,nwork=4*Nx
        INTEGER(kind=4),PARAMETER               :: isign=0,isignf=-1,isignb=1
        REAL(kind=8),PARAMETER                  :: scalef=1.0d0,scaleb=1.0d0/DBLE(Nx)
        REAL(kind=8), DIMENSION(1:ntable), SAVE :: tabledz,tablezd
        REAL(kind=8), DIMENSION(1:nwork)        :: work
        INTEGER(kind=4), SAVE                   :: isys
        CHARACTER*100                           :: name_config
        CALL MPI_INIT(ierr)        
        ! set up fft
        CALL DZFFT(isign,Nx,scalef,u,v,tabledz,work,isys)
        CALL ZDFFT(isign,Nx,scaleb,v,u,tablezd,work,isys)

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
        DO i=1,Nx/2 + 1
           nonlinhat(i)=v(i)
        END DO
        CALL ZDFFT(isignb,Nx,scaleb,nonlinhat,nonlin,tablezd,work,isys) 
        PRINT *,'Transformed initialdata'
        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                        enkin(plotnum),enstr(plotnum),&
                        enpot(plotnum),en(plotnum),&
                        kx,nonlin,nonlinhat,&
                        v,vold,u,uold,tabledz,tablezd,work)
        DO i=1,Nx
           unew(i)=u(i)
           utemp(i)=u(i)
        END DO
			
	PRINT *,'Got initial data, starting timestepping'
	time(plotnum)=0.0d0
        start = MPI_WTIME()
	DO n=1,(Nt/3)	
               ! Iteration 1		
                nn=nn+1		
		DO i=1,Nx
		       nonlinfix(i)=Es*u(i)**3&
                             +(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                              (3.0d0*(u(i)**2)*(-2.0d0*u(i)+uold(i))+&
                               6.0d0*u(i)*(-1.0d0*u(i))*(u(i)-uold(i)))
		END DO
                CALL DZFFT(isignf,Nx,scalef,nonlinfix,nonlinfixhat,&
                           tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        nonlinfixhat(i)=nonlinfixhat(i)+&
                             (kx(i)*kx(i) -1.0d0)*v(i) + &
                             (2.0d0*v(i)-vold(i))*(1.0d0/(dt*dt)-&
                             (kx(i)*kx(i) -1.0d0)*2.0d0/(4.0d0*3.0d0*2.0d0)) 
                END DO

                chg=1.0d0
                iter=1
                DO WHILE ((chg .gt. tol) .and. (iter .lt. maxiter))
                   DO i=1,Nx
                        utemp(i)=unew(i)
                        nonlin(i)=(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                                  (3.0d0*(u(i)**2)*unew(i)+&
                                   6.0d0*u(i)*unew(i)*(u(i)-uold(i)))
                   END DO
                   CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,&
                              tabledz,work,isys)
		   DO i=1,Nx/2 + 1
			vnew(i)=(nonlinhat(i)+nonlinfixhat(i) ) &
			/(1.0d0/(dt*dt) - (2.0d0/(4.0d0*3.0d0*2.0d0))*(kx(i)*kx(i) - 1.0d0))
		   END DO
                   CALL ZDFFT(isignb,Nx,scaleb,vnew,unew,tablezd,work,isys)
                   chg=0.0d0
		   ! normalize result
		   DO i=1,Nx
                        chg=chg+(unew(i)-utemp(i))**2
		   END DO
                   chg=SQRT(chg/REAL(Nx,kind(0d0)))
                   iter=iter+1
                END DO		
		IF (mod(nn,plotgap)==0) then
			plotnum=plotnum+1
			time(plotnum)=nn*dt
			PRINT *,'time',nn*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                      enkin(plotnum),enstr(plotnum),&
                                      enpot(plotnum),en(plotnum),&
                                      kx,nonlin,nonlinhat,&
                                      vnew,v,unew,u,tabledz,tablezd,work)
			!CALL savedata(Nx,plotnum,name_config,unew)
		END IF
                ! Iteration 2
                nn=nn+1
                DO i=1,Nx
                       nonlinfix(i)=Es*unew(i)**3&
                             +(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                              (3.0d0*(unew(i)**2)*(-2.0d0*unew(i)+u(i))+&
                               6.0d0*unew(i)*(-1.0d0*unew(i))*(unew(i)-u(i)))
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlinfix,nonlinfixhat,&
                           tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        nonlinfixhat(i)=nonlinfixhat(i)+&
                             (kx(i)*kx(i) -1.0d0)*vnew(i) + &
                             (2.0d0*vnew(i)-v(i))*(1.0d0/(dt*dt)-&
                             (kx(i)*kx(i)-1.0d0)*2.0d0/(4.0d0*3.0d0*2.0d0))
                END DO

                chg=1.0d0
                iter=1
                DO WHILE ((chg .gt. tol) .and. (iter .lt. maxiter))
                   DO i=1,Nx
                        utemp(i)=uold(i)
                        nonlin(i)=(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                                  (3.0d0*(unew(i)**2)*uold(i)+&
                                   6.0d0*unew(i)*uold(i)*(unew(i)-u(i)))
                   END DO
                   CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,&
                              tabledz,work,isys)
                   DO i=1,Nx/2 + 1
                        vold(i)=(nonlinhat(i)+nonlinfixhat(i) ) &
                        /(1.0d0/(dt*dt) - (2.0d0/(4.0d0*3.0d0*2.0d0))*(kx(i)*kx(i) - 1.0d0))
                   END DO
                   CALL ZDFFT(isignb,Nx,scaleb,vold,uold,tablezd,work,isys)
                   chg=0.0d0
                   ! normalize result
                   DO i=1,Nx
                        chg=chg+(uold(i)-utemp(i))**2
                   END DO
                   chg=SQRT(chg/REAL(Nx,kind(0d0)))
                   iter=iter+1
                END DO          
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=nn*dt
                        PRINT *,'time',nn*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                      enkin(plotnum),enstr(plotnum),&
                                      enpot(plotnum),en(plotnum),&
                                      kx,nonlin,nonlinhat,&
                                      vold,vnew,uold,unew,tabledz,tablezd,work)
                        !CALL savedata(Nx,plotnum,name_config,unew)
                END IF
                ! Iteration 3
                nn=nn+1
                DO i=1,Nx
                       nonlinfix(i)=Es*uold(i)**3&
                             +(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                              (3.0d0*(uold(i)**2)*(-2.0d0*uold(i)+unew(i))+&
                               6.0d0*uold(i)*(-1.0d0*uold(i))*(uold(i)-unew(i)))
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlinfix,nonlinfixhat,&
                           tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        nonlinfixhat(i)=nonlinfixhat(i)+&
                             (kx(i)*kx(i) -1.0d0)*vold(i) + &
                             (2.0d0*vold(i)-vnew(i))*(1.0d0/(dt*dt)-&
                             (kx(i)*kx(i)-1.0d0)*2.0d0/(4.0d0*3.0d0*2.0d0))
                END DO

                chg=1.0d0
                iter=1
                DO WHILE ((chg .gt. tol) .and. (iter .lt. maxiter))
                   DO i=1,Nx
                        utemp(i)=u(i)
                        nonlin(i)=(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                                  (3.0d0*(uold(i)**2)*u(i)+&
                                   6.0d0*uold(i)*u(i)*(uold(i)-unew(i)))
                   END DO
                   CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,&
                              tabledz,work,isys)
                   DO i=1,Nx/2 + 1
                        v(i)=(nonlinhat(i)+nonlinfixhat(i) ) &
                        /(1.0d0/(dt*dt) - (2.0d0/(4.0d0*3.0d0*2.0d0))*(kx(i)*kx(i) - 1.0d0))
                   END DO
                   CALL ZDFFT(isignb,Nx,scaleb,v,u,tablezd,work,isys)
                   chg=0.0d0
                   ! normalize result
                   DO i=1,Nx
                        chg=chg+(u(i)-utemp(i))**2
                   END DO
                   chg=SQRT(chg/REAL(Nx,kind(0d0)))
                   iter=iter+1
                END DO
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=nn*dt
                        PRINT *,'time',nn*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                      enkin(plotnum),enstr(plotnum),&
                                      enpot(plotnum),en(plotnum),&
                                      kx,nonlin,nonlinhat,&
                                      v,vold,u,uold,tabledz,tablezd,work)
                        !CALL savedata(Nx,plotnum,name_config,unew)
                END IF
	END DO	
        IF(mod(Nt,3)>0) THEN
               ! Iteration 1            
                nn=nn+1         
                DO i=1,Nx
                       nonlinfix(i)=Es*u(i)**3&
                             +(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                              (3.0d0*(u(i)**2)*(-2.0d0*u(i)+uold(i))+&
                               6.0d0*u(i)*(-1.0d0*u(i))*(u(i)-uold(i)))
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlinfix,nonlinfixhat,&
                           tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        nonlinfixhat(i)=nonlinfixhat(i)+&
                             (kx(i)*kx(i) -1.0d0)*v(i) + &
                             (2.0d0*v(i)-vold(i))*(1.0d0/(dt*dt)-&
                             (kx(i)*kx(i)-1.0d0)*2.0d0/(4.0d0*3.0d0*2.0d0))
                END DO

                chg=1.0d0
                iter=1
                DO WHILE ((chg .gt. tol) .and. (iter .lt. maxiter))
                   DO i=1,Nx
                        utemp(i)=unew(i)
                        nonlin(i)=(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                                  (3.0d0*(u(i)**2)*unew(i)+&
                                   6.0d0*u(i)*unew(i)*(u(i)-uold(i)))
                   END DO
                   CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,&
                              tabledz,work,isys)
                   DO i=1,Nx/2 + 1
                        vnew(i)=(nonlinhat(i)+nonlinfixhat(i) ) &
                        /(1.0d0/(dt*dt) - (2.0d0/(4.0d0*3.0d0*2.0d0))*(kx(i)*kx(i) - 1.0d0))
                   END DO
                   CALL ZDFFT(isignb,Nx,scaleb,vnew,unew,tablezd,work,isys)
                   chg=0.0d0
                   ! normalize result
                   DO i=1,Nx
                        chg=chg+(unew(i)-utemp(i))**2
                   END DO
                   chg=SQRT(chg/REAL(Nx,kind(0d0)))
                   iter=iter+1
                END DO          
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=nn*dt
                        PRINT *,'time',nn*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                      enkin(plotnum),enstr(plotnum),&
                                      enpot(plotnum),en(plotnum),&
                                      kx,nonlin,nonlinhat,&
                                      vnew,v,unew,u,tabledz,tablezd,work)
                        !CALL savedata(Nx,plotnum,name_config,unew)
                END IF
        ENDIF
        IF(mod(Nt,3)>1) THEN
                ! Iteration 2
                nn=nn+1
                DO i=1,Nx
                       nonlinfix(i)=Es*unew(i)**3&
                             +(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                              (3.0d0*(unew(i)**2)*(-2.0d0*unew(i)+u(i))+&
                               6.0d0*unew(i)*(-1.0d0*unew(i))*(unew(i)-u(i)))
                END DO
                CALL DZFFT(isignf,Nx,scalef,nonlinfix,nonlinfixhat,&
                           tabledz,work,isys)
                DO i=1,Nx/2 + 1
                        nonlinfixhat(i)=nonlinfixhat(i)+&
                             (kx(i)*kx(i) -1.0d0)*vnew(i) + &
                             (2.0d0*vnew(i)-v(i))*(1.0d0/(dt*dt)-&
                             (kx(i)*kx(i)-1.0d0)*2.0d0/(4.0d0*3.0d0*2.0d0))
                END DO

                chg=1.0d0
                iter=1
                DO WHILE ((chg .gt. tol) .and. (iter .lt. maxiter))
                   DO i=1,Nx
                        utemp(i)=uold(i)
                        nonlin(i)=(2.0d0/(4.0d0*3.0d0*2.0d0))*&
                                  (3.0d0*(unew(i)**2)*uold(i)+&
                                   6.0d0*unew(i)*uold(i)*(unew(i)-u(i)))
                   END DO
                   CALL DZFFT(isignf,Nx,scalef,nonlin,nonlinhat,&
                              tabledz,work,isys)
                   DO i=1,Nx/2 + 1
                        vold(i)=(nonlinhat(i)+nonlinfixhat(i) ) &
                        /(1.0d0/(dt*dt) -(2.0d0/(4.0d0*3.0d0*2.0d0))*(kx(i)*kx(i) - 1.0d0))
                   END DO
                   CALL ZDFFT(isignb,Nx,scaleb,vold,uold,tablezd,work,isys)
                   chg=0.0d0
                   ! normalize result
                   DO i=1,Nx
                        chg=chg+(uold(i)-utemp(i))**2
                   END DO
                   chg=SQRT(chg/REAL(Nx,kind(0d0)))
                   iter=iter+1
                END DO
                IF (mod(nn,plotgap)==0) then
                        plotnum=plotnum+1
                        time(plotnum)=nn*dt
                        PRINT *,'time',nn*dt
                        CALL enercalc(Nx,isignf,isignb,isys,ntable,nwork,dt,Es,&
                                      enkin(plotnum),enstr(plotnum),&
                                      enpot(plotnum),en(plotnum),&
                                      kx,nonlin,nonlinhat,&
                                      vold,vnew,uold,unew,tabledz,tablezd,work)
                        !CALL savedata(Nx,plotnum,name_config,unew)
                END IF
        END IF	
        finish=MPI_WTIME()
        PRINT*,'Finished timestepping'
        runtime=REAL(finish-start,kind(0d0))
	PRINT*,'Program took ',runtime,'for Time stepping'       
	CALL saveresults(Nt,plotgap,time(1:1+n/plotgap),en(1:1+n/plotgap),&
			enstr(1:1+n/plotgap),enkin(1:1+n/plotgap),&
                        enpot(1:1+n/plotgap))
			
	! Save times at which output was made in text format
	PRINT *,'Saved data'
        IF(mod(Nt,3)==0) CALL errorcalc(Nt,Nx,dt,c,x,u,err)
        IF(mod(Nt,3)==1) CALL errorcalc(Nt,Nx,dt,c,x,unew,err)
        IF(mod(Nt,3)==2) CALL errorcalc(Nt,Nx,dt,c,x,uold,err)

        name_config = 'data/errortime'
        CALL savetimeanderror(Nx,Nt,err,runtime,name_config)
        
	PRINT *,'Program execution complete'
        CALL MPI_FINALIZE(ierr)
	END PROGRAM Kg
