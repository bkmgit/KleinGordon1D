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
	!  dt			= timestep
	!  ierr			= error code
	!  plotnum		= number of plot
	! .. Arrays ..
	!  unew 		= approximate solution
	!  u 			= approximate solution
	!  uold 		= approximate solution
	!  nonlin 		= nonlinear term, u^3
	! .. Vectors ..
	!  x			= x locations
	!  time			= times at which save data
	!  en			= total energy	
	!  enstr		= strain energy
	!  enpot		= potential energy
	!  enkin		= kinetic energy
	!  name_config		= array to store filename for data to be saved    		
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
	! 	OpenMP library
		
	PROGRAM Kg
	USE omp_lib		 	   
	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4), PARAMETER		:: Nx=1024
	INTEGER(kind=4), PARAMETER		:: Nt=500 
	INTEGER(kind=4), PARAMETER		:: plotgap=32	
        INTEGER(kind=4), PARAMETER		:: maxiter=50
	REAL(kind=8), PARAMETER			:: &
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER			::  Lx=9.0d0
	REAL(kind=8), PARAMETER			::  Es=1.0d0	
        REAL(kind=8), PARAMETER         	::  c=0.5
        REAL(kind=8), PARAMETER 		::  tol=0.1**12
	REAL(kind=8)				::  dt=5.0d0/REAL(Nt,kind(0d0))	
	REAL(kind=8), DIMENSION(1:Nx)	        ::  r,rhs,Mp,p 
	REAL(kind=8), DIMENSION(1:Nx)	        ::  x	
	REAL(kind=8), DIMENSION(1:Nx)           ::  u
	REAL(kind=8), DIMENSION(1:Nx)           ::  uold
	REAL(kind=8), DIMENSION(1:Nx)           ::  unew
	REAL(kind=8), DIMENSION(1:Nx)           ::  savearray
	REAL(kind=8), DIMENSION(1:Nt+1)         ::  time,enkin,enstr,enpot,en
        REAL(kind=8)                            ::  rsum,rr,rrold,alpha,beta,dx,err,runtime,start,finish
	INTEGER(kind=4)				::  ierr,i,j,n,allocatestatus,numthreads
	INTEGER(kind=4)				::  count_rate, plotnum
        INTEGER(kind=4)                         ::  ip,ipp,ippp,ipppp,im,imm,immm,immmm,iter
	CHARACTER*100				::  name_config
	! Start short parallel region to count threads 
	!numthreads=omp_get_max_threads()
	!PRINT *,'There are ',numthreads,' threads.'
		
	! setup grid 
	CALL getgrid(Nx,Lx,pi,name_config,x)
        dx=abs(x(2)-x(1))
	PRINT *,'Setup grid'
	CALL initialdata(Nx,dt,c,x,u,uold)
	plotnum=1	
	name_config = 'data/u' 
	!savearray=REAL(u)
	!CALL savedata(Nx,plotnum,name_config,savearray)

        PRINT *,'Transformed initial data'
	CALL enercalc(Nx,dt,Es,&
			enkin(plotnum),enstr(plotnum),&
			enpot(plotnum),en(plotnum),&
			dx,u,uold)
			
	PRINT *,'Got initial data, starting timestepping'
	time(plotnum)=0.0d0
  	!CALL system_clock(start,count_rate)	
        start = OMP_get_wtime()
	DO n=1,Nt					
                rr=0.0d0
                rsum=0.0d0
		!$OMP PARALLEL DO PRIVATE(i,ipp,ip,im,imm) reduction(+:rr,rsum) SCHEDULE(static)
		DO i=1,Nx
                        ipppp=1+mod(i+4+Nx-1,Nx)
                        ippp =1+mod(i+3+Nx-1,Nx)
                        ipp  =1+mod(i+2+Nx-1,Nx)
                        ip   =1+mod(i+1+Nx-1,Nx)
                        im   =1+mod(i-1+Nx-1,Nx)
                        imm  =1+mod(i-2+Nx-1,Nx) 
                        immm =1+mod(i-3+Nx-1,Nx)
                        immmm=1+mod(i-4+Nx-1,Nx)
			rhs(i)=Es*(u(i)**3) &
                              +(2.0d0*u(i)-uold(i))/(dt*dt) &
                              +0.25d0*( -(1.0d0/560.0d0)*(2.0d0*u(immmm)+uold(immmm)) &
                                        +(8.0d0/315.0d0)*(2.0d0*u(immm)+uold(immm)) &
                                        -(1.0d0/5.0d0)*(2.0d0*u(imm)+uold(imm)) &
                                        +(8.0d0/5.0d0)*(2.0d0*u(im)+uold(im)) &
                                        -(205.0d0/72.0d0)*(2.0d0*u(i)+uold(i)) &
                                        +(8.0d0/5.0d0)*(2.0d0*u(ip)+uold(ip)) &
                                        -(1.0d0/5.0d0)*(2.0d0*u(ipp)+uold(ipp)) &
                                        +(8.0d0/315.0d0)*(2.0d0*u(ippp)+uold(ippp)) & 
                                        -(1.0d0/560.0d0)*(2.0d0*u(ipppp)+uold(ipppp)))/(dx*dx) &
                              -0.25d0*(2.0d0*u(i)+uold(i))
                        r(i)=rhs(i)-u(i)*(1.0d0/(dt*dt) + 0.25d0) &
                             +0.25d0*(-(1.0d0/560.0d0)*u(immmm) &
                                      +(8.0d0/315.0d0)*u(immm) &
                                      -(1.0d0/5.0d0)*u(imm) &
                                      +(8.0d0/5.0d0)*u(im) &
                                      -(205.0d0/72.0d0)*u(i) &
                                      +(8.0d0/5.0d0)*u(ip) &
                                      -(1.0d0/5.0d0)*u(ipp)&
                                      +(8.0d0/315.0d0)*u(ippp) &
                                      -(1.0d0/560.0d0)*u(ipppp) )/(dx*dx) 
                        p(i)=r(i)
                        unew(i)=u(i)
                        rsum=rsum+abs(r(i))
                        rr=rr+r(i)*r(i)
		END DO
		!$OMP END PARALLEL DO
                iter=0
                DO WHILE (( iter .lt. maxiter ) .and. ( rsum .gt. tol))
		  !$OMP PARALLEL DO PRIVATE(i,im,imm,ip,ipp) SCHEDULE(static)
		  DO i=1,Nx
                        ipppp=1+mod(i+4+Nx-1,Nx)
                        ippp =1+mod(i+3+Nx-1,Nx)
                        ipp  =1+mod(i+2+Nx-1,Nx)
                        ip   =1+mod(i+1+Nx-1,Nx)
                        im   =1+mod(i-1+Nx-1,Nx)
                        imm  =1+mod(i-2+Nx-1,Nx)
                        immm =1+mod(i-3+Nx-1,Nx)
                        immmm=1+mod(i-4+Nx-1,Nx)
			Mp(i)=p(i)*(1.0d0/(dt*dt) + 0.25d0) &
                             -0.25d0*(-(1.0d0/560.0d0)*p(immmm) &
                                      +(8.0d0/315.0d0)*p(immm) &
                                      -(1.0d0/5.0d0)*p(imm) &
                                      +(8.0d0/5.0d0)*p(im) &
                                      -(205.0d0/72.0d0)*p(i) &
                                      +(8.0d0/5.0d0)*p(ip) &
                                      -(1.0d0/5.0d0)*p(ipp) &
                                      +(8.0d0/315.0d0)*p(ippp) &
                                      -(1.0d0/560.0d0)*p(ipppp) )/(dx*dx)
		  END DO
		  !$OMP END PARALLEL DO
                  alpha=0.0d0
                  !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:alpha) SCHEDULE(static)
                  DO i=1,Nx
                    alpha=alpha+Mp(i)*p(i)
                  END DO
                  !$OMP END PARALLEL DO
                  alpha=rr/alpha

                  !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
                  DO i=1,Nx
                    unew(i)=unew(i)+alpha*p(i)
                    r(i)=r(i)-alpha*Mp(i)
                  END DO
                  !$OMP END PARALLEL DO

                  rrold=rr
                  rr=0.0d0
                  !$OMP PARALLEL DO PRIVATE(i,im,imm,ip,ipp) REDUCTION(+:rr) SCHEDULE(static)
                  DO i=1,Nx               
                    rr=rr+r(i)*r(i)
                  END DO
                  !$OMP END PARALLEL DO
                  beta=rr/rrold
                  rsum=0.0d0
                  !$OMP PARALLEL DO PRIVATE(i,im,imm,ip,ipp) REDUCTION(+:rsum) SCHEDULE(static)
                  DO i=1,Nx
                    p(i)=r(i)+beta*p(i)
                    rsum=rsum+abs(r(i))
                  END DO
                  !$OMP END PARALLEL DO
                  iter=iter+1    
                END DO
                ! .. Update old values ..
                CALL storeold(Nx,unew,u,uold)
		IF (mod(n,plotgap)==0) then
			plotnum=plotnum+1
			time(plotnum)=n*dt
			PRINT *,'time',n*dt,' iterations ',iter
			CALL enercalc(Nx,dt,Es,&
				enkin(plotnum),enstr(plotnum),&
				enpot(plotnum),en(plotnum),dx,&
				u,uold)
			!CALL savedata(Nx,plotnum,name_config,unew)
		END IF
	END DO		
	PRINT *,'Finished time stepping'
	!CALL system_clock(finish,count_rate)
	finish = OMP_get_wtime()
        runtime=REAL(finish-start,kind(0d0))
	PRINT*,'Program took ',runtime,'for Time stepping'
	CALL saveresults(Nt,plotgap,time(1:1+n/plotgap),en(1:1+n/plotgap),&
			enstr(1:1+n/plotgap),enkin(1:1+n/plotgap),&
                        enpot(1:1+n/plotgap))	
			
	! Save times at which output was made in text format
	PRINT *,'Saved data'
        CALL errorcalc(Nt,Nx,dt,c,x,u,err)
        name_config = 'data/errortime'
        CALL savetimeanderror(Nx,Nt,err,runtime,name_config)
	PRINT *,'Program execution complete'
	END PROGRAM Kg
