! Cedrick Ansorge, cedrick@posteo.de 15 Mar 2015  
! FORTRAN CODE FOR VARIANCE DECOMPOSITION AMONG DIMENSIONS OF DATA
! ANalysis Of Variance 
! 
PROGRAM ANOVA_MAIN 
  USE ANOVA_CONSTANTS
  USE ANOVA
  USE nc_readwrite
  USE timer 

  IMPLICIT NONE  

  INTEGER :: ivar,idim,i,j,nx,ny,nz,nt 
  CHARACTER nc_filename_input*128 
  CHARACTER :: vname*128,outbase*256
  CHARACTER(LEN=20) :: xdim,ydim,zdim,tdim,dummy
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: d_arr  
  REAL(KIND=8) :: si_residual
  REAL :: si_max 
  INTEGER :: fixed_pos =1
  INTEGER, PARAMETER :: verbose =1 
  LOGICAL, DIMENSION(DLAST) :: write_decomp
  
  TYPE(T_ANOVA)  anova3, anova4  
  
  IF (iargc() .LT. 6) THEN  
     WRITE(*,*) 'ERROR (ANOVA_MAIN): FILE OR VARIABLE NOT PROVIDED' 
     WRITE(*,*) 'usage ./anova.x <filename> <varname> <xdim> <ydim> <zdim> <tdim>'  
     STOP 'ERROR (ANOVA_MAIN)' 
  ELSEIF (iargc() .GT. 6 ) THEN 
     WRITE(*,*) 'WARNING (ANOVA_MAIN): following arguments ignored:' 
     DO i=7,iargc()  
        CALL GETARG(i,dummy)
        WRITE(*,*) 'WARNING (ANOVA_MAIN): ', i, dummy
     ENDDO
  ENDIF

  CALL GETARG(1,nc_filename_input) 
  CALL GETARG(2,vname)
  CALL GETARG(3,xdim)
  CALL GETARG(4,ydim)
  CALL GETARG(5,zdim)
  CALL GETARG(6,tdim)

  CALL ncrw_init(nc_filename_input,verbose)

  IF ( ncrw_inq_dimtranspose(vname,ydim,zdim) ) THEN 
     WRITE(*,*) 'WARNING (ANOVA_MAIN): ', 'Need to exchange ', TRIM(ydim),' and ',TRIM(zdim), ' to avoid tranposition' 
     dummy=zdim 
     zdim=ydim 
     ydim=dummy
  ENDIF

  write_decomp(:) = .TRUE. 


  ! INITIALIZE, CALCULATE AND OUTPUT 3D ANOVA 
  CALL ANOVA_INIT(vname,xdim,ydim,zdim,tdim,anova3,3) 
  CALL ANOVA_DECOMP3D(vname,xdim,ydim,zdim,tdim,fixed_pos,anova3)  
  outbase=TRIM(vname)//'.si3d'
  CALL ANOVA_OUTPUT_ASC(anova3,outbase)  
  outbase=TRIM(vname)//'.anova3d.'
  CALL ANOVA_OUTPUT_BIN(anova3,outbase,write_decomp) 

  CALL ANOVA_INIT(vname,xdim,ydim,zdim,tdim,anova4,4) 
  CALL ANOVA_DECOMP4D(vname,si_residual,xdim,ydim,zdim,tdim,anova4)  
  outbase=TRIM(vname)//'.si4d'
  CALL ANOVA_OUTPUT_ASC(anova4,outbase)
  outbase=TRIM(vname)//'.anova4d.'
  CALL ANOVA_OUTPUT_BIN(anova4,outbase,write_decomp)

END PROGRAM ANOVA_MAIN

