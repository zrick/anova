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

  INTEGER :: ivar,iarg,idim,i,j,nx,ny,nz,nt,nvar,narg
  CHARACTER nc_filename_input*128 
  CHARACTER :: vname*128,outbase*256
  CHARACTER(LEN=20) :: dim1,    dim2,    dim3,    dim4,    dummy, trange
  CHARACTER(LEN=20) :: xdim_loc,ydim_loc,zdim_loc,tdim_loc 
  CHARACTER(LEN=64), DIMENSION(:), ALLOCATABLE :: var_list 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: d_arr  
  REAL(KIND=8) :: si_residual
  REAL(KIND=8), DIMENSION(:,:),ALLOCATABLE :: si_prof 
  REAL :: si_max 
  INTEGER :: fixed_pos, fixed_len, it_srt, it_end
  INTEGER, PARAMETER :: verbose =0
  LOGICAL, DIMENSION(DLAST) :: write_decomp
  
  TYPE(T_ANOVA)  anova3, anova4  

  narg=iargc()
  IF (narg .LT. 6) THEN  
     WRITE(*,*) 'ERROR (ANOVA_MAIN): FILE OR VARIABLE NOT PROVIDED' 
     WRITE(*,*) 'usage ./anova.x <filename> <xdim> <ydim> <zdim> <tdim> <var1name> [<var2name> ... <varnname>]'  
     STOP 'ERROR (ANOVA_MAIN)' 
  ENDIF

  nvar=narg-5
  ALLOCATE(var_list(nvar))

  CALL GETARG(1,nc_filename_input) 
  CALL GETARG(2,dim1)
  CALL GETARG(3,dim2)
  CALL GETARG(4,dim3)
  CALL GETARG(5,dim4) 

  
  ivar=1
  DO iarg = 6,narg
     CALL GETARG(iarg,var_list(ivar))
     ivar=ivar+1 
  ENDDO

  WRITE(*,*) 'PROCESSING', nvar, 'VARIABLES:', var_list

  CALL ncrw_init(nc_filename_input,verbose)
  WRITE(*,*) 'INITIALIZED MODULE NC_READWRITE'
  write_decomp(:) = .TRUE. 
  
  IF ( ncrw_inq_dimtranspose(var_list(1),dim2,dim3) ) THEN
     WRITE(*,*) 'WARNING (ANOVA_MAIN): ', 'Exchanging ', &
          TRIM(dim2),' and ',TRIM(dim3), ' to avoid tranposition for 3D decomp'
     xdim_loc=dim1; ydim_loc=dim3; zdim_loc=dim2; tdim_loc=dim4
  ELSE
     xdim_loc=dim1; ydim_loc=dim2; zdim_loc=dim3; tdim_loc=dim4
  ENDIF 
  WRITE(*,*) 'INITIALIZING 3D ANOVA' 
  CALL ANOVA_INIT(xdim_loc,ydim_loc,zdim_loc,tdim_loc,anova3,3,1,5)  
  ! 
  fixed_len=ncrw_getdimlen(tdim_loc)
  ALLOCATE(si_prof(DLAST,fixed_len))
  !
  DO ivar=1,nvar
     vname = var_list(ivar) 
     DO fixed_pos=anova3%it_srt,anova3%it_end ! fixed_len
        WRITE(*,*) 'Calling ANOVA_DECOMP3D for VAR:',TRIM(ADJUSTL(vname)),' POS:',fixed_pos
        CALL ANOVA_DECOMP3D(vname,fixed_pos,anova3) 
        si_prof(:,fixed_pos) = anova3%si(:)
     ENDDO
     outbase=TRIM(vname)//'.si3d.'//TRIM(tdim_loc)//'-prof'
     CALL ANOVA_OUTPUT_ASCP(anova3,outbase,si_prof)
  ENDDO

  ! INITIALIZE, CALCULATE AND OUTPUT 4D ANOVA
  IF ( ncrw_inq_dimtranspose(var_list(1),dim2,dim3) ) THEN
     WRITE(*,*) 'WARNING (ANOVA_MAIN): ', 'Exchanging ', &
          TRIM(dim2),' and ',TRIM(dim3), ' to avoid tranposition for 4D decomp'
     xdim_loc=dim1; ydim_loc=dim3; zdim_loc=dim2; tdim_loc=dim4
  ELSE
     xdim_loc=dim1; ydim_loc=dim2; zdim_loc=dim3; tdim_loc=dim4
  ENDIF

  CALL ANOVA_INIT(xdim_loc,ydim_loc,zdim_loc,tdim_loc,anova4,4,5,20)
  DO ivar=1,nvar 
     vname=var_list(ivar)  
     CALL ANOVA_DECOMP4D(vname,si_residual,anova4)
     outbase=TRIM(vname)//'.si4d'
     CALL ANOVA_OUTPUT_ASC(anova4,outbase)
     outbase=TRIM(vname)//'.anova4d'
     CALL ANOVA_OUTPUT_BIN(anova4,outbase,write_decomp)
  ENDDO

END PROGRAM ANOVA_MAIN
