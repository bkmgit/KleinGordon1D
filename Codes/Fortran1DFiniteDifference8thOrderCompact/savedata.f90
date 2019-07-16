	SUBROUTINE savedata(Nx,plotnum,name_config,field)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine saves a one dimensional real array in binary 
	! format
	!
	! INPUT
	!
	! .. Scalars ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  plotnum			= number of plot to be made
	! .. Arrays ..
	!  field 			= real data to be saved
	!  name_config			= root of filename to save to 
	!
	! .. Output	..	
	! plotnum			= number of plot to be saved
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  count			= counter
	!  iol				= size of file
	! .. Arrays ..
	! 	number_file		= array to hold the number of the plot
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
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)			:: Nx
	INTEGER(KIND=4), INTENT(IN)			:: plotnum
	REAL(KIND=8), DIMENSION(1:NX), INTENT(IN)	:: field
	CHARACTER*100, INTENT(IN)			:: name_config
	INTEGER(kind=4)					:: i,iol,count,ind
	CHARACTER*100					:: number_file

	! create character array with full filename
	ind = index(name_config,' ') - 1
	WRITE(number_file,'(i0)') 9999999+plotnum
	number_file = name_config(1:ind)//number_file
	ind = index(number_file,' ') - 1
	number_file = number_file(1:ind)//'.datbin'	
	INQUIRE( iolength=iol )	field(1)
	OPEN(unit=11,FILE=number_file,form="unformatted",&
		access="direct",recl=iol) 	
	count=1
	DO i=1,Nx
		WRITE(11,rec=count) field(i)
		count=count+1
	END DO
	CLOSE(11)
	
	END SUBROUTINE savedata
