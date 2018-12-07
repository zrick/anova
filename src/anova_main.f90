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
  CHARACTER(LEN=64), DIMENSION(NCRW_MAXVAR) :: var_list 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: d_arr  
  REAL(KIND=8) :: si_residual
  REAL(KIND=8), DIMENSION(:,:),ALLOCATABLE :: si_prof 
  REAL :: si_max 
  INTEGER                  :: fixed_pos, fixed_len 
  INTEGER, DIMENSION(2)    :: range1_3d,range4_3d,range1_4d,range4_4d
  INTEGER, PARAMETER       :: verbose=1
  LOGICAL, DIMENSION(DLAST):: write_decomp 
  LOGICAL                  :: decomp_3d, decomp_4d 

  NAMELIST /anova_global/ nc_filename_input,dim1, dim2, dim3, dim4, & 
       nvar, decomp_3d, decomp_4d 
  NAMELIST /anova_ranges/ range1_3d,range4_3d,range1_4d,range4_4d
  NAMELIST /anova_vars/   var_list 

  TYPE(T_ANOVA)  anova3, anova4   


  OPEN(10,file='anova.nml') 
  READ(10,nml=anova_global)  
  READ(10,nml=anova_vars)  
  READ(10,nml=anova_ranges)


  WRITE(*,*) 'PROCESSING', nvar, 'VARIABLES:', var_list(1:nvar)

  CALL ncrw_init(nc_filename_input,verbose)
  WRITE(*,*) 'INITIALIZED MODULE NC_READWRITE'
  write_decomp(:) = .TRUE. 

  IF ( decomp_3d .EQV. .TRUE. ) THEN
     IF ( ncrw_inq_dimtranspose(var_list(1),dim2,dim3) ) THEN
        WRITE(*,*) 'WARNING (ANOVA_MAIN): ', 'Exchanging ', &
             TRIM(dim2),' and ',TRIM(dim3), ' to avoid tranposition for 3D decomp'
        xdim_loc=dim1; ydim_loc=dim3; zdim_loc=dim2; tdim_loc=dim4
     ELSE
        xdim_loc=dim1; ydim_loc=dim2; zdim_loc=dim3; tdim_loc=dim4
     ENDIF
     WRITE(*,*) 'INITIALIZING 3D ANOVA' 
     WRITE(*,*) 'SUBRANGE IN 3D:', range4_3d(1), range4_3d(2) 
     CALL ANOVA_INIT(xdim_loc,ydim_loc,zdim_loc,tdim_loc,anova3,3,&
          ix_srt=range1_3d(1),ix_end=range1_3d(2), &
          it_srt=range4_3d(1),it_end=range4_3d(2)) 
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
  ENDIF

  IF ( decomp_4d .EQV. .TRUE. ) THEN 
     ! INITIALIZE, CALCULATE AND OUTPUT 4D ANOVA
     IF ( ncrw_inq_dimtranspose(var_list(1),dim2,dim3) ) THEN
        WRITE(*,*) 'WARNING (ANOVA_MAIN): ', 'Exchanging ', &
             TRIM(dim2),' and ',TRIM(dim3), ' to avoid tranposition for 4D decomp'
        xdim_loc=dim1; ydim_loc=dim3; zdim_loc=dim2; tdim_loc=dim4
     ELSE
        xdim_loc=dim1; ydim_loc=dim2; zdim_loc=dim3; tdim_loc=dim4
     ENDIF
     
     CALL ANOVA_INIT(xdim_loc,ydim_loc,zdim_loc,tdim_loc,anova4,4,& 
          ix_srt=range1_4d(1),ix_end=range1_4d(2), &
          it_srt=range4_4d(1),it_end=range4_4d(2))

     DO ivar=1,nvar 
        vname=var_list(ivar)  
        CALL ANOVA_DECOMP4D(vname,si_residual,anova4)
        outbase=TRIM(vname)//'.si4d'
        CALL ANOVA_OUTPUT_ASC(anova4,outbase)
        outbase=TRIM(vname)//'.anova4d'
        CALL ANOVA_OUTPUT_BIN(anova4,outbase,write_decomp)
     ENDDO
  ENDIF

END PROGRAM ANOVA_MAIN
