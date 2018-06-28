MODULE ANOVA  

  USE ANOVA_CONSTANTS
  USE NC_READWRITE 
  USE TIMER 

  TYPE T_ANOVA 
     INTEGER :: nx,ny,nz,nt,ndim  
     CHARACTER(LEN=20)                          :: xdim,ydim,zdim,tdim,vname
     CHARACTER(LEN=20)                          :: fname 
     REAL(KIND=8)                               :: f_empty,si_residual
     REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: f_x,f_y,f_z,f_t 
     REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: f_yx,f_yz,f_yt,f_xz,f_xt,f_zt
     REAL(KIND=8), DIMENSION(:,:,:),ALLOCATABLE :: f_yxz,f_yxt,f_yzt,f_xzt 
     REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: si 
     LOGICAL,      DIMENSION(:),    ALLOCATABLE :: f_save 
  END TYPE T_ANOVA

  CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ANOVA BINARY OUTPUT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE ANOVA_OUTPUT_BIN(a,fbase,wrt)     
    
    IMPLICIT NONE 
    
    ! PARAMETERS 
    TYPE(T_ANOVA), INTENT(IN) :: a 
    LOGICAL, DIMENSION(*) :: wrt 
    CHARACTER(LEN=20) :: xdim,ydim,zdim,tdim  
    INTEGER, PARAMETER :: funit=28  
    INTEGER :: nd
    CHARACTER(LEN=*) :: fbase
    CHARACTER(LEN=200) :: fname 

    xdim=a%xdim(:3); ydim=a%ydim(:3); zdim=a%zdim(:3); tdim=a%tdim(:3) 
    nd=a%ndim
    IF ( wrt(DX) .EQV. .TRUE. ) THEN   
       fname=TRIM(fbase)//TRIM(xdim)
       WRITE(*,*) fbase,xdim,fname
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_x(:) 
       CLOSE(funit)
    ENDIF
    IF ( wrt(DY) .EQV. .TRUE. ) THEN   
       fname=TRIM(fbase)//TRIM(ydim)
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_y(:) 
       CLOSE(funit) 
    ENDIF
    IF ( wrt(DZ) .EQV. .TRUE. ) THEN   
       fname=TRIM(fbase)//TRIM(zdim)
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_z(:) 
       CLOSE(funit) 
    ENDIF
    IF ( wrt(DT) .EQV. .TRUE. .AND. nd .GT. 3) THEN   
       fname=TRIM(fbase)//TRIM(tdim)
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_t(:) 
       CLOSE(funit) 
    ENDIF
    !
    ! 2-DIMENSIONAL FIELDS 
    !
    IF ( wrt(DYX) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase)//TRIM(ydim)//'_'//TRIM(xdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yx(:,:) 
    ENDIF
    IF ( wrt(DYZ) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase)//TRIM(ydim)//'_'//TRIM(zdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yz(:,:) 
    ENDIF
    IF ( wrt(DYT) .EQV. .TRUE. .AND. nd .GT. 3) THEN 
       fname=TRIM(fbase)//TRIM(ydim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yt(:,:) 
    ENDIF
    IF ( wrt(DXZ) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase)//TRIM(xdim)//'_'//TRIM(zdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_xz(:,:) 
    ENDIF
    IF ( wrt(DXT) .EQV. .TRUE. .AND.nd.GT. 3) THEN 
       fname=TRIM(fbase)//TRIM(xdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_xt(:,:) 
    ENDIF
    IF ( wrt(DZT) .EQV. .TRUE. .AND.nd.GT.3) THEN 
       fname=TRIM(fbase)//TRIM(zdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_zt(:,:) 
    ENDIF
    !
    ! 3-DIMENSIONAL FIELDS 
    !
    IF ( wrt(DYXZ) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase)//TRIM(ydim)//'_'//TRIM(xdim)//'_'//TRIM(zdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yxz(:,:,:) 
    ENDIF
    IF ( wrt(DYXT) .EQV. .TRUE. .AND. nd.GT.3) THEN 
       fname=TRIM(fbase)//TRIM(ydim)//'_'//TRIM(xdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yxt(:,:,:) 
    ENDIF
    IF ( wrt(DYZT) .EQV. .TRUE. .AND. nd.GT.3) THEN 
       fname=TRIM(fbase)//TRIM(ydim)//'_'//TRIM(zdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yzt(:,:,:) 
    ENDIF
    IF ( wrt(DXZT) .EQV. .TRUE. .AND. nd.GT.3) THEN 
       fname=TRIM(fbase)//TRIM(xdim)//'_'//TRIM(zdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_xzt(:,:,:) 
    ENDIF
    
    
  END SUBROUTINE ANOVA_OUTPUT_BIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ANOVA ASCII OUTPUT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE ANOVA_OUTPUT_ASC(a,fname)    
    
    IMPLICIT NONE  

    ! PARAMETERS 
    TYPE(T_ANOVA), INTENT(IN):: a  

    ! LOCAL DEFINITIONS 
    CHARACTER(LEN=*)         :: fname  
    INTEGER,       PARAMETER :: funit=27 
    INTEGER                  :: ndim
    REAL(KIND=8)             :: si_max 

    ndim = a%ndim 
    
    IF ( ndim.LT.3 .OR. ndim .GT. 4 ) & 
       STOP 'ANOVA NOT IMPLEMENTED FOR LESS THAN 3 OR MORE THAN 4 DIMENSIONS' 

    si_max=MAXVAL(a%si(2:DLAST))

    OPEN(funit,FILE=fname,form='FORMATTED') 

    WRITE(funit,106) 'DIM', 'VARIANCE', 'SENSITIVITY','SENSITIVITY-2' 
    WRITE(funit,104) 'TOTAL VAR',   a%si(Dtot),1.0 
    WRITE(funit,104) 'RESIDUAL VAR',a%si_residual,a%si_residual/a%si(DTOT) 
    IF(ndim.EQ.4)WRITE(funit,105) a%ydim,a%xdim,a%zdim,a%tdim,a%si(Dyxzt),&
         100.*a%si(Dyxzt)/a%si(Dtot),100.*a%si(DYXZT)/(a%si(Dtot)-si_max)
    WRITE(funit,105)              a%ydim,a%xdim,a%zdim,'',    a%si(Dyxz), &
         100.*a%si(Dyxz)/a%si(Dtot), 100.*a%si(DYXZ)/(a%si(Dtot)-si_max)
    IF(ndim.EQ.4)WRITE(funit,105) a%ydim,a%xdim,a%tdim,'',    a%si(Dyxt), &
         100.*a%si(Dyxt)/a%si(Dtot), 100.*a%si(DYXT)/(a%si(Dtot)-si_max)
    IF(ndim.EQ.4)WRITE(funit,105) a%ydim,a%zdim,a%tdim,'',    a%si(Dyzt), &
         100.*a%si(Dyzt)/a%si(Dtot), 100.*a%si(DYZT)/(a%si(Dtot)-si_max)
    IF(ndim.EQ.4)WRITE(funit,105) a%xdim,a%zdim,a%tdim,'',    a%si(Dxzt), &
         100.*a%si(Dxzt)/a%si(Dtot), 100.*a%si(DXZT)/(a%si(Dtot)-si_max)
    WRITE(funit,105)              a%ydim,a%xdim,'','',        a%si(Dyx),  &
         100.*a%si(Dyx)/a%si(Dtot),  100.*a%si(DYX)/(a%si(Dtot)-si_max)
    WRITE(funit,105)              a%ydim,a%zdim,'','',        a%si(Dyz),  &
         100.*a%si(Dyz)/a%si(Dtot),  100.*a%si(DYZ)/(a%si(Dtot)-si_max)
    IF(ndim.EQ.4)WRITE(funit,105) a%ydim,a%tdim,'','',        a%si(Dyt),  &
         100.*a%si(Dyt)/a%si(Dtot),  100.*a%si(DYT)/(a%si(Dtot)-si_max)
    WRITE(funit,105)              a%xdim,a%zdim,'','',        a%si(Dxz),  &
         100.*a%si(Dxz)/a%si(Dtot),  100.*a%si(DXZ)/(a%si(Dtot)-si_max)
    IF(ndim.EQ.4)WRITE(funit,105) a%xdim,a%tdim,'','',        a%si(Dxt),  &
         100.*a%si(Dxt)/a%si(Dtot),  100.*a%si(DXT)/(a%si(Dtot)-si_max)
    IF(ndim.EQ.4)WRITE(funit,105) a%zdim,a%tdim,'','',        a%si(Dzt),  &
         100.*a%si(Dzt)/a%si(Dtot),  100.*a%si(DZT)/(a%si(Dtot)-si_max)
    WRITE(funit,105)              a%xdim,'',  '','',          a%si(Dx),   &
         100.*a%si(Dx)/a%si(Dtot),   100.*a%si(DX)/(a%si(Dtot)-si_max)
    WRITE(funit,105)              a%ydim,'',  '','',          a%si(Dy),   &
         100.*a%si(Dy)/a%si(Dtot),   100.*a%si(DY)/(a%si(Dtot)-si_max)
    WRITE(funit,105)              a%zdim,'',  '','',          a%si(Dz),   &
         100.*a%si(Dz)/a%si(Dtot),   100.*a%si(DZ)/(a%si(Dtot)-si_max)
    IF(ndim.EQ.4)WRITE(funit,105) a%tdim,'',  '','',          a%si(Dt),   &
         100.*a%si(Dt)/a%si(Dtot),   100.*a%si(DT)/(a%si(Dtot)-si_max)
  ! 
    CLOSE(funit) 
106 FORMAT('#',a3,';',a9,  ';',a11,';',a13)
105 FORMAT(a3,'-',a3,'-',a3,'-',a3,';',g10.3,';',g11.3'%;',g11.3,'%;') 
104 FORMAT(a15,';',g10.3,';',g11.3';') 
  END SUBROUTINE ANOVA_OUTPUT_ASC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ANOVA 3D - variance decomposition and sensitivity indices 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ANOVA_DECOMP3D(vname,xdim,ydim,zdim,fixed_dim,fixed_pos,an) 

    IMPLICIT NONE 

    INTEGER,          INTENT(IN) :: fixed_pos
    CHARACTER(LEN=*), INTENT(IN) :: vname,xdim,ydim,zdim,fixed_dim 
    TYPE(T_ANOVA)                :: an

    REAL(KIND=8)                              :: avg_yxz,rdum               !,f_empty
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: avg_yx, avg_yz,avg_xz      !,f_x,f_y,f_z  
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: avg_x,avg_y,avg_z,data,wrk !,f_yx,f_yz,f_xz   
    ! REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::                            !f_yxz
    INTEGER,      DIMENSION(2)                :: pos 
    CHARACTER,    DIMENSION(2)                :: yz_dims*32,xt_dims*32  

    INTEGER :: ix,iy,iz
    INTEGER :: nx,ny,nz,nt,nyz,nyx,nxz,nyxz 
    REAL :: chk_x_vs_y, chk_x_vs_z, chk_x_vs_yx, chk_x_vs_xz, chk_x_vs_yz
    REAL :: chk_y_vs_z,             chk_y_vs_yx, chk_y_vs_xz, chk_y_vs_yz
    REAL ::                         chk_z_vs_yx, chk_z_vs_xz, chk_z_vs_yz 
    REAL :: chk_yx_vs_yz, chk_yx_vs_xz, chk_yz_vs_xz 
    REAL :: max_err 

    CALL TIMER_START() 

    nx=ncrw_getdimlen(xdim);      an%nx=nx
    ny=ncrw_getdimlen(ydim);      an%ny=ny 
    nz=ncrw_getdimlen(zdim);      an%nz=nz 
    nt=ncrw_getdimlen(fixed_dim); an%nt=nt 

    IF ( ncrw_verbose ) THEN 
       WRITE(*,*) '============'
       WRITE(*,*) 'ANOVA_DECOMP3D: VARIABLE ',TRIM(vname),' NX=',nx,' NY=',ny,' NZ=',nz, ' FIXED_POS=',fixed_pos 
    ENDIF
    nyz=ny*nz 
    nxz=nx*nz 
    nyx=ny*nx
    nyxz=ny*nx*nz
    yz_dims(1) = ydim;  yz_dims(2) = zdim
    xt_dims(1) = xdim;  xt_dims(2) = fixed_dim


    ALLOCATE(avg_yx(nz),avg_yz(nx),avg_xz(ny)) 
    ALLOCATE(avg_x(ny,nz), avg_y(nx,nz), avg_z(ny,nx),data(ny,nz),wrk(ny,nz)) 
    
    CALL TIMER_FINISH('ANOVA_DECOMP3D: init time elapse')
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 1: INTEGRATION 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL TIMER_START() 

    IF ( ncrw_verbose ) & 
         WRITE(*,*) 'ANOVA_DECOMP3D: STEP 1 - INTEGRATION' 
 
    avg_yxz=0. 
    avg_yx(:)=0.; avg_yz(:)=0.; avg_xz(:)=0.
    avg_x(:,:)=0.;avg_y(:,:)=0.;avg_z(:,:)=0.
    pos(2) = fixed_pos
    DO ix=1,nx 
       pos(1)=ix  
       CALL ncrw_getvar_slice(vname,yz_dims,xt_dims,pos,data)    
       rdum=SUM(data)
       avg_yz(ix)=avg_yz(ix)+rdum/nyz
       avg_yxz   =avg_yxz   +rdum/nyz/nx   
       
       avg_x = avg_x+data/nx 

       DO iy=1,ny  
          rdum=SUM(data(iy,:))
          avg_xz(iy)=avg_xz(iy)+rdum/nxz
          avg_z(iy,ix)=avg_z(iy,ix)+rdum/nz 
       ENDDO

       DO iz=1,nz 
          rdum=SUM(data(:,iz))
          avg_yx(iz)=avg_yx(iz)+rdum/nyx
          avg_y(ix,iz)=avg_y(ix,iz)+rdum/ny
       ENDDO
    ENDDO

    ! Check if all fields average out to the same as the total data
    IF ( SUM(avg_yx)/nz-avg_yxz .GT. SMALL_DELTA .OR. &  
         SUM(avg_yz)/nx-avg_yxz .GT. SMALL_DELTA .OR. & 
         SUM(avg_xz)/ny-avg_yxz .GT. SMALL_DELTA .OR. & 
         SUM(avg_x)/nyz-avg_yxz .GT. SMALL_DELTA .OR. & 
         SUM(avg_y)/nxz-avg_yxz .GT. SMALL_DELTA .OR. & 
         SUM(avg_z)/nyx-avg_yxz .GT. SMALL_DELTA )  THEN  
       WRITE(*,*) 'yx:',SUM(avg_yx/nz) - avg_yxz
       WRITE(*,*) 'yz:',SUM(avg_yz/nx) - avg_yxz
       WRITE(*,*) 'xz:',SUM(avg_xz/ny) - avg_yxz
       WRITE(*,*) 'x: ',SUM(avg_x/nyz) - avg_yxz
       WRITE(*,*) 'y: ',SUM(avg_y/nxz) - avg_yxz
       WRITE(*,*) 'z: ',SUM(avg_z/nyx) - avg_yxz
       STOP 'CONSERVATION PROBLEM IN DECOMPOSITION'   
    ELSE IF ( ncrw_verbose ) THEN  
       WRITE(*,*) 'ANOVA_DECOMP3D: STEP 1 Finished, Integral Calculus consistent'  
       CALL TIMER_FINISH('ANOVA_DECOMP3D: STEP 1 Time elapse') 
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 2: DECOMPOSITION 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    CALL TIMER_START() 
    IF ( ncrw_verbose ) WRITE(*,*) 'ANOVA_DECOMP3D: STEP 2 - DECOMPOSITION AND SENSITIVITY INDICES' 
    an%f_empty=avg_yxz 
    ! 
    an%f_x=avg_yz-an%f_empty
    an%f_y=avg_xz-an%f_empty 
    an%f_z=avg_yx-an%f_empty  
    ! 
    DO iy=1,ny 
       DO ix=1,nx; an%f_yx(iy,ix)=avg_z(iy,ix)-an%f_x(ix)-an%f_y(iy)-an%f_empty;ENDDO 
       DO iz=1,nz; an%f_yz(iy,iz)=avg_x(iy,iz)-an%f_z(iz)-an%f_y(iy)-an%f_empty;ENDDO
    ENDDO
    DO ix=1,nx; 
       DO iz=1,nz; an%f_xz(ix,iz)=avg_y(ix,iz)-an%f_z(iz)-an%f_x(ix)-an%f_empty;ENDDO  
    ENDDO 

    an%si(:DLAST) = 0. 
    DO ix=1,nx  
       pos(1)=ix
       CALL ncrw_getvar_slice(vname,yz_dims,xt_dims,pos,data)   
       DO iy=1,ny; DO iz=1,nz  
          an%f_yxz(iy,ix,iz) = data(iy,iz) - an%f_yx(iy,ix) - an%f_xz(ix,iz) - an%f_yz(iy,iz) & 
               - an%f_x(ix) - an%f_y(iy) - an%f_z(iz) - an%f_empty
       ENDDO; ENDDO 
       data = data(:,:) - an%f_empty  
       an%si(DTOT) = an%si(DTOT) + ( SUM(data * data) ) / nyxz
       an%si(DYXZ) = an%si(DYXZ) + ( SUM(an%f_yxz(:,ix,:)*an%f_yxz(:,ix,:)) ) / nyxz
    ENDDO

    an%si(DX)  = SUM(an%f_x*an%f_x)/nx 
    an%si(DY)  = SUM(an%f_y*an%f_y)/ny
    an%si(DZ)  = SUM(an%f_z*an%f_z)/nz
    an%si(DYX) = SUM(an%f_yx*an%f_yx)/nyx
    an%si(DYZ) = SUM(an%f_yz*an%f_yz)/nyz 
    an%si(DXZ) = SUM(an%f_xz*an%f_xz)/nxz  
    an%si_residual = an%si(DTOT) - an%si(DYXZ)-an%si(DYX)-an%si(DYZ)-an%si(DXZ)-an%si(DX)-an%si(DY)-an%si(DZ) 

    IF ( ABS ( an%si_residual/an%si(DTOT)  ) .GT. 1e-6 ) THEN 
       WRITE(STDOUT,*) 'ANOVA_DECOMP3D: VARIANCE DECOMPOSITION FAILED, RESIDUAL:', an%si_residual
       STOP 'ANOVA_DECOMP3D: VARIANCE DECOMPOSITION FAILED' 
    ELSEIF ( ncrw_verbose ) THEN 
       WRITE(STDOUT,*) 'ANOVA_DECOMP3D: VARIANCE DECOMPOSED TO WITHIN', ABS(an%si_residual) / an%si(DTOT)  
    ENDIF  
104 FORMAT(A30,1x,G10.3) 
    !  
    ! Check orthogonality 
    ! xy 
    chk_x_vs_y=0.; chk_x_vs_z=0.; chk_x_vs_yx=0.; chk_x_vs_xz=0.; chk_x_vs_yz=0.
    chk_y_vs_z=0.;                chk_y_vs_yx=0.; chk_y_vs_xz=0.; chk_y_vs_yz=0.
    ;                             chk_z_vs_yx=0.; chk_z_vs_xz=0.; chk_z_vs_yz=0.
    chk_yx_vs_yz=0.; chk_yx_vs_xz=0.; chk_yz_vs_xz=0.; 
    DO ix=1,nx 
       DO iy=1,ny; 
          chk_x_vs_y = chk_x_vs_y  +an%f_x(ix)*an%f_y(iy); 
          chk_x_vs_yx = chk_x_vs_yx+an%f_x(ix)*an%f_yx(iy,ix)  
          chk_y_vs_yx = chk_y_vs_yx+an%f_y(iy)*an%f_yx(iy,ix)
          DO iz=1,nz  
             chk_z_vs_yx = chk_z_vs_yx + an%f_z(iz)*an%f_yx(iy,ix) 
             chk_x_vs_yz = chk_x_vs_yz + an%f_x(ix)*an%f_yz(iy,iz) 
             chk_y_vs_xz = chk_y_vs_xz + an%f_y(iy)*an%f_xz(ix,iz) 
             chk_yx_vs_yz= chk_yx_vs_yz+ an%f_yx(iy,ix)*an%f_yz(iy,iz)  
             chk_yx_vs_xz= chk_yx_vs_xz+ an%f_yx(iy,ix)*an%f_xz(ix,iz) 
             chk_yz_vs_xz= chk_yz_vs_xz+ an%f_yz(iy,iz)*an%f_xz(ix,iz) 
          ENDDO
       ENDDO
       DO iz=1,nz
          chk_x_vs_z = chk_x_vs_z + an%f_x(ix)*an%f_z(iz); 
          chk_x_vs_xz= chk_x_vs_xz+ an%f_x(ix)*an%f_xz(ix,iz); 
          chk_z_vs_xz= chk_z_vs_xz+ an%f_z(iz)*an%f_xz(ix,iz); 
       ENDDO
    ENDDO
    DO iy=1,ny 
       DO iz=1,nz 
          chk_y_vs_z = chk_y_vs_z + an%f_y(iy)*an%f_z(iz) 
          chk_y_vs_yz= chk_y_vs_yz+ an%f_y(iy)*an%f_yz(iy,iz) 
          chk_z_vs_yz= chk_z_vs_yz+ an%f_z(iz)*an%f_yz(iy,iz) 
       ENDDO 
    ENDDO

    chk_x_vs_y=chk_x_vs_y/nyx 
    chk_x_vs_yx=chk_x_vs_yx/nyx 
    chk_y_vs_yx=chk_y_vs_yx/nyx 
    chk_z_vs_yx=chk_z_vs_yx/nyxz
    chk_x_vs_yz=chk_x_vs_yz/nyxz
    chk_y_vs_xz=chk_y_vs_xz/nyxz
    chk_yx_vs_yz=chk_yx_vs_yz/nyxz 
    chk_yx_vs_xz=chk_yx_vs_xz/nyxz 
    chk_yz_vs_xz=chk_yz_vs_xz/nyxz 
    chk_x_vs_z=chk_x_vs_z/nxz
    chk_x_vs_xz=chk_x_vs_xz/nxz 
    chk_z_vs_xz=chk_z_vs_xz/nxz 
    chk_y_vs_z=chk_y_vs_z/nyz 
    chk_y_vs_yz=chk_y_vs_yz/nyz 
    chk_z_vs_yz=chk_z_vs_yz/nyz

    max_err = MAXVAL ( (/chk_x_vs_y, chk_x_vs_yx, chk_y_vs_yx, chk_z_vs_yx, chk_x_vs_yz,chk_y_vs_xz,&
         chk_yx_vs_yz, chk_yx_vs_xz, chk_yz_vs_xz,chk_x_vs_z,chk_x_vs_xz,chk_z_vs_xz,chk_y_vs_z,& 
         chk_y_vs_yz,chk_z_vs_yz/) )

    IF ( max_err /an%si(DTOT)  .GT. SQRT(SMALL_DELTA) ) THEN 
       WRITE(*,*) 'ANOVA_DECOMP3D: Orthogonality of ANOVA decomposition violated by', max_err /an%si(DTOT)
       STOP 'ANOVA_DECOMP3D: ORTHOGONALITY VIOLATED'
    ENDIF

    IF ( ncrw_verbose ) THEN 
       WRITE(*,*) 'ANOVA_DECOMP3D: GLOBAL AVERAGE:', an%f_empty 
!       WRITE(*,*) 'YZ:',SUM(an%f_yz*an%f_yz)/nyz, MINVAL(an%f_yz),MAXVAL(an%f_yz), SUM(an%f_yz)/nyz,nyz
!       WRITE(*,*) 'YX:',SUM(an%f_yx*an%f_yx)/nyx, MINVAL(an%f_yx),MAXVAL(an%f_yx), SUM(an%f_yx)/nyx,nyx
!       WRITE(*,*) 'XZ:',SUM(an%f_xz*an%f_xz)/nxz, MINVAL(an%f_xz),MAXVAL(an%f_xz), SUM(an%f_xz)/nxz,nxz
!       WRITE(*,*) 'X: ',SUM(an%f_x*an%f_x)/nx, MINVAL(an%f_x),MAXVAL(an%f_x),      SUM(an%f_x)/nx,  nx
!       WRITE(*,*) 'Y: ',SUM(an%f_y*an%f_y)/ny, MINVAL(an%f_y),MAXVAL(an%f_y),      SUM(an%f_y)/ny,  ny
!       WRITE(*,*) 'Z: ',SUM(an%f_z*an%f_z)/nz, MINVAL(an%f_z),MAXVAL(an%f_z),      SUM(an%f_z)/nz,  nz  
!       WRITE(*,*) ''
       WRITE(*,*) 'ANOVA_DECOMP3D: ORTHOGONALITY OF COMPONENTS ACHIEVED (MAXIMUM ERROR:',max_err,')'
!        WRITE(*,*) 'X VS Y  ', chk_x_vs_y, SQRT(d_arr(DX)*d_arr(DY))
!        WRITE(*,*) 'X VS Z  ', chk_x_vs_y, SQRT(d_arr(DX)*d_arr(DZ)) 
!        WRITE(*,*) 'X VS YX ', chk_x_vs_yx,SQRT(d_arr(DX)*d_arr(DYX)) 
!        WRITE(*,*) 'X VS YZ ', chk_x_vs_yz,SQRT(d_arr(DX)*d_arr(DYZ)) 
!        WRITE(*,*) 'X VS XZ ', chk_x_vs_xz,SQRT(d_arr(DX)*d_arr(DXZ)) 
!        WRITE(*,*) 'Y VS Z  ', chk_y_vs_z, SQRT(d_arr(DY)*d_arr(DZ)) 
!        WRITE(*,*) 'Y VS YX ', chk_y_vs_yx,SQRT(d_arr(DY)*d_arr(DYX))  
!        WRITE(*,*) 'Y VS YZ ', chk_y_vs_yz,SQRT(d_arr(DY)*d_arr(DYZ)) 
!        WRITE(*,*) 'Y VS XZ ', chk_y_vs_xz,SQRT(d_arr(DY)*d_arr(DXZ)) 
!        WRITE(*,*) 'Z VS YX ', chk_z_vs_yx,SQRT(d_arr(DZ)*d_arr(DYX)) 
!        WRITE(*,*) 'Z VS YZ ', chk_z_vs_yz,SQRT(d_arr(DZ)*d_arr(DYZ)) 
!        WRITE(*,*) 'Z VS XZ ', chk_z_vs_xz,SQRT(d_arr(DZ)*d_arr(DXZ)) 
!        WRITE(*,*) 'YX VS YZ ', chk_yx_vs_yz,SQRT(d_arr(DYX)*d_arr(DYZ)) 
!        WRITE(*,*) 'YX VS XZ ', chk_yx_vs_xz,SQRT(d_arr(DYX)*d_arr(DXZ)) 
!        WRITE(*,*) 'YZ VS XZ ', chk_yz_vs_xz,SQRT(d_arr(DYZ)*d_arr(DXZ)) 
    ENDIF 

    CALL TIMER_FINISH('ANOVA_DECOMP3D: STEP2 Time elapse') 

  END SUBROUTINE ANOVA_DECOMP3D   


  SUBROUTINE ANOVA_INIT(vname,xdim,ydim,zdim,tdim,a,ndim)   
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: vname,xdim,ydim,zdim,tdim   
    INTEGER,          INTENT(IN) :: ndim
    INTEGER :: nx,ny,nz,nt 
    TYPE(T_ANOVA) :: a 

    a%ndim=ndim 
    a%nx=ncrw_getdimlen(xdim);  a%xdim=xdim; nx=a%nx
    a%ny=ncrw_getdimlen(ydim);  a%ydim=ydim; ny=a%ny 
    a%nz=ncrw_getdimlen(zdim);  a%zdim=zdim; nz=a%nz 
    a%nt=ncrw_getdimlen(tdim);  a%tdim=tdim; nt=a%nt
    

    IF ( ndim .GE. 3 ) THEN 
       ALLOCATE(a%f_x(nx), a%f_y(ny), a%f_z(nz)) 
       ALLOCATE(a%f_xz(nx,nz), a%f_yz(ny,nz), a%f_yx(ny,nx))
       ALLOCATE(a%f_yxz(ny,nx,nz))
       ALLOCATE(a%si(DLAST))
       IF ( ndim .EQ. 4 ) THEN 
          ALLOCATE(a%f_t(nt)) 
          ALLOCATE(a%f_zt(nz,nt), a%f_xt(nx,nt), a%f_yt(ny,nt))
          ALLOCATE(a%f_yzt(ny,nz,nt), a%f_xzt(nx,nz,nt), a%f_yxt(ny,nx,nt))
       ELSE IF ( ndim .GT. 4 ) THEN 
          STOP 'ANOVA_INIT: MORE THAN 4D ANOVA NOT IMPLEMENTED'
       ENDIF
    ELSE 
       STOP 'ANOVA_INIT: LESS THAN 3D ANOVANOT IMPLEMENTED' 
    ENDIF
  END SUBROUTINE ANOVA_INIT

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ANOVA 4D - Variance decomposition and sensitivity indices 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  SUBROUTINE ANOVA_DECOMP4D(vname,si_res,xdim,ydim,zdim,tdim,a) 
    IMPLICIT NONE 

    CHARACTER(LEN=*), INTENT(IN) :: vname,xdim,ydim,zdim,tdim 
    REAL(KIND=8)                                :: si_res   
    TYPE(T_ANOVA) :: a
    ! 
    REAL(KIND=8)                                :: avg_yxzt,rdum,f_empty
    REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: avg_yxz,avg_yzt,avg_xzt,avg_yxt 
    REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: f_t,    f_x,    f_y,    f_z
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: avg_yx, avg_yz, avg_yt, avg_xz, avg_xt, avg_zt
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: f_zt,   f_xt,   f_xz,   f_yt,   f_yz,   f_yx, data 
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: avg_x,  avg_y,  avg_z,  avg_t
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: f_yzt,  f_xzt,  f_yxt,  f_yxz
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: f_yxzt
    INTEGER,      DIMENSION(2)                  :: pos 
    CHARACTER,    DIMENSION(2)                  :: yz_dims*32,xt_dims*32  

    INTEGER :: ix,iy,iz,it 
    INTEGER :: nx,ny,nz,nt
    INTEGER :: nyz,nyx,nyt,nxz,nxt,nzt
    INTEGER :: nyzt,nxzt,nyxt,nyxz 
    INTEGER :: nyxzt 

    CALL TIMER_START()

    nx=a%nx!ncrw_getdimlen(xdim);
    ny=a%ny!ncrw_getdimlen(ydim);
    nz=a%nz!ncrw_getdimlen(zdim); 
    nt=a%nt!ncrw_getdimlen(tdim); 

    
    IF ( ncrw_verbose ) THEN 
       WRITE(*,*) '=============='
       WRITE(*,*) 'ANOVA_DECOMP4D: VARIABLE ',TRIM(vname),' NX=',nx,' NY=',ny,' NZ=',nz, ' NT=',nt 
    ENDIF
    nyz=ny*nz;     nyx=ny*nx;     nyt=ny*nt;     nxz=nx*nz;     nxt=nx*nt; nzt=nz*nt  
    nyzt=ny*nz*nt; nxzt=nx*nz*nt; nyxt=ny*nx*nt; nyxz=ny*nx*nz
    nyxzt=ny*nx*nz*nt

    yz_dims(1) = ydim;  yz_dims(2) = zdim
    xt_dims(1) = xdim;  xt_dims(2) = tdim
    !
    ALLOCATE(avg_yxz(nt),avg_yzt(nx),avg_xzt(ny),avg_yxt(nz)) 
    ALLOCATE(avg_yx(nz,nt),avg_yz(nx,nt),avg_yt(nx,nz),avg_xz(ny,nt),avg_xt(ny,nz),avg_zt(ny,nx))  
    ALLOCATE(avg_x(ny,nz,nt), avg_y(nx,nz,nt), avg_z(ny,nx,nt),avg_t(ny,nx,nz)) 
    ALLOCATE(data(ny,nz),f_yxzt(ny,nz)) 
    
    CALL TIMER_FINISH('ANOVA_DECOMP4D: INIT time elapse') 

    CALL TIMER_START() 
    IF ( ncrw_verbose ) THEN
       WRITE(*,*) 'ANOVA_DECOMP4D: STEP 1 -- CALCULATING INTEGRALS' 
    ENDIF

    avg_yxzt=0. 
    avg_yxz(:)=0.;  avg_yzt(:)=0.;   avg_xzt(:)=0.;   avg_yxt(:)=0. 
    avg_yx(:,:)=0.; avg_yz(:,:)=0.;  avg_yt(:,:)=0.; avg_xz(:,:)=0.; avg_xt(:,:)=0.; avg_zt(:,:)=0.
    avg_x(:,:,:)=0.;avg_y(:,:,:)=0.; avg_z(:,:,:)=0.;avg_t(:,:,:)=0.; 
    a%f_empty=0.; a%f_x(:)=0.; a%f_y(:)=0.; a%f_z(:)=0.; a%f_t(:)=0.; 
    a%f_yx(:,:)=0.; a%f_yz(:,:)=0.; a%f_yt(:,:)=0.; a%f_xz(:,:)=0.; a%f_xt(:,:)=0.; a%f_zt(:,:)=0. 
    a%f_yxz(:,:,:)=0.; a%f_yzt(:,:,:)=0.; a%f_xzt(:,:,:)=0.; a%f_yxt(:,:,:)=0. 
    DO it=1,nt;  DO ix=1,nx 
       pos(1)=ix  
       pos(2)=it 
       CALL ncrw_getvar_slice(vname,yz_dims,xt_dims,pos,data)    
       rdum=SUM(data)
       avg_yz(ix,it)=avg_yz(ix,it)+rdum/nyz
       avg_yxz(it)=avg_yxz(it) +rdum/nyxz
       avg_yzt(ix)=avg_yzt(ix) +rdum/nyzt
       avg_yxzt=avg_yxzt+rdum/nyxzt
       
       avg_x(:,:,it) = avg_x(:,:,it)+data/nx  
       avg_t(:,ix,:) = avg_t(:,ix,:)+data/nt
       
       DO iy=1,ny  
          rdum=SUM(data(iy,:)) 
          avg_xzt(iy)  =avg_xzt(iy)+rdum/nxzt
          avg_xz(iy,it)=avg_xz(iy,it)+rdum/nxz  
          avg_zt(iy,ix)=avg_zt(iy,ix)+rdum/nzt
          avg_z(iy,ix,it)=avg_z(iy,ix,it)+rdum/nz 
          DO iz=1,nz
             avg_xt(iy,iz)=avg_xt(iy,iz)+data(iy,iz)/nxt 
          ENDDO
       ENDDO
          

       DO iz=1,nz 
          rdum=SUM(data(:,iz))
          avg_yxt(iz)=avg_yxt(iz)+rdum/nyxt
          avg_yx(iz,it)=avg_yx(iz,it)+rdum/nyx 
          avg_yt(ix,iz)=avg_yt(ix,iz)+rdum/nyt
          avg_y(ix,iz,it)=avg_y(ix,iz,it)+rdum/ny
       ENDDO
    ENDDO;ENDDO

    ! Check if all fields average out to the same as the total data
    IF ( SUM(avg_yxt/nz) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_yzt/nx) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_xzt/ny) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_yxz/nt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_yx/nzt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_yz/nxt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_y/nxzt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_xz/nyt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_xt/nyz) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_zt/nyx) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_x/nyzt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_y/nxzt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_z/nyxt) - avg_yxzt .GT. SMALL_DELTA .OR. &
         SUM(avg_t/nyxz) - avg_yxzt .GT. SMALL_DELTA ) THEN
       STOP 'ANOVA_DECOMP4D: CONSERVATION PROBLEM IN DECOMPOSITION'   
    ELSE 
       IF ( ncrw_verbose ) & 
            WRITE(*,*) 'ANOVA_DECOMP4D: STEP 1 FINISHED, INTEGRAL CALCULUS CONSISTENT' 
    ENDIF
    CALL TIMER_FINISH('ANOVA_DECOMP4D: STEP1 time elapse') 
    CALL TIMER_START() 


    IF ( ncrw_verbose ) & 
         WRITE(*,*) 'ANOVA_DECOMP4D: STEP 2 -- DECOMPOSITION' 
    
    a%f_empty = avg_yxzt 
    a%f_x = avg_yzt - a%f_empty 
    a%f_y = avg_xzt - a%f_empty 
    a%f_z = avg_yxt - a%f_empty 
    a%f_t = avg_yxz - a%f_empty
    
    DO iy=1,ny
       DO ix=1,nx; a%f_yx(iy,ix)=avg_zt(iy,ix)-a%f_y(iy)-a%f_x(ix)-a%f_empty; ENDDO 
       DO iz=1,nz; a%f_yz(iy,iz)=avg_xt(iy,iz)-a%f_y(iy)-a%f_z(iz)-a%f_empty; ENDDO 
       DO it=1,nt; a%f_yt(iy,it)=avg_xz(iy,it)-a%f_y(iy)-a%f_t(it)-a%f_empty; ENDDO 
    ENDDO 
    DO ix=1,nx
       DO iz=1,nz; a%f_xz(ix,iz)=avg_yt(ix,iz)-a%f_x(ix)-a%f_z(iz)-a%f_empty; ENDDO 
       DO it=1,nt; a%f_xt(ix,it)=avg_yz(ix,it)-a%f_x(ix)-a%f_t(it)-a%f_empty; ENDDO 
    ENDDO 

    DO iz=1,nz
       DO it=1,nt; a%f_zt(iz,it)=avg_yx(iz,it)-a%f_z(iz)-a%f_t(it)-a%f_empty; ENDDO 
    ENDDO

    DO iy=1,ny
       DO ix=1,nx; DO iz=1,nz
          a%f_yxz(iy,ix,iz)=avg_t(iy,ix,iz)-a%f_y(iy)-a%f_x(ix)-a%f_z(iz)-a%f_yx(iy,ix)-a%f_yz(iy,iz)-a%f_xz(ix,iz)-a%f_empty
       ENDDO; ENDDO 
       DO ix=1,nx; DO it=1,nt
          a%f_yxt(iy,ix,it)=avg_z(iy,ix,it)-a%f_y(iy)-a%f_x(ix)-a%f_t(it)-a%f_yx(iy,ix)-a%f_yt(iy,it)-a%f_xt(ix,it)-a%f_empty
       ENDDO; ENDDO 
       DO iz=1,nz; DO it=1,nt
          a%f_yzt(iy,iz,it)=avg_x(iy,iz,it)-a%f_y(iy)-a%f_z(iz)-a%f_t(it)-a%f_yz(iy,iz)-a%f_yt(iy,it)-a%f_zt(iz,it)-a%f_empty  
       ENDDO;ENDDO
    ENDDO

    DO ix=1,nx     
       DO iz=1,nz; DO it=1,nt
          a%f_xzt(ix,iz,it)=avg_y(ix,iz,it)-a%f_x(ix)-a%f_z(iz)-a%f_t(it)-a%f_xz(ix,iz)-a%f_xt(ix,it)-a%f_zt(iz,it)-a%f_empty  
       ENDDO;ENDDO 
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! STEP 3 Total variance and residual 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    a%si(:DLAST) = 0. 
    DO it=1,nt; DO ix=1,nx  
       pos(1)=ix 
       pos(2)=it 
       CALL ncrw_getvar_slice(vname,yz_dims,xt_dims,pos,data)   
       data = data(:,:) - a%f_empty 
       DO iy=1,ny; DO iz=1,nz  
          f_yxzt(iy,iz) = data(iy,iz) &
               - a%f_yxz(iy,ix,iz) - a%f_yxt(iy,ix,it) - a%f_yzt(iy,iz,it) - a%f_xzt(ix,iz,it) &
               - a%f_yx(iy,ix) - a%f_yz(iy,iz) - a%f_yt(iy,it) - a%f_xz(ix,iz) - a%f_xt(ix,it) - a%f_zt(iz,it) & 
               - a%f_x(ix) - a%f_y(iy) - a%f_z(iz) - a%f_t(it) ! -f_empty already above 
       ENDDO; ENDDO 
       a%si(DTOT)=a%si(DTOT) + ( SUM(data * data) ) 
       a%si(DYXZT)=a%si(DYXZT) + ( SUM(f_yxzt*f_yxzt) ) 
    ENDDO; ENDDO 

    a%si(DTOT) = a%si(DTOT) / nyxzt
    a%si(DYXZT)= a%si(DYXZT)/ nyxzt

    IF ( ncrw_verbose ) WRITE(*,*) 'ANOVA_DECOMP4D: Global average ', a%f_empty 

    a%si(DX)=SUM(a%f_x*a%f_x)/nx
    a%si(DY)=SUM(a%f_y*a%f_y)/ny
    a%si(DZ)=SUM(a%f_z*a%f_z)/nz
    a%si(DT)=SUM(a%f_t*a%f_t)/nt
    a%si(DYX)=SUM(a%f_yx*a%f_yx)/nyx
    a%si(DYZ)=SUM(a%f_yz*a%f_yz)/nyz
    a%si(DYT)=SUM(a%f_yt*a%f_yt)/nyt
    a%si(DXZ)=SUM(a%f_xz*a%f_xz)/nxz
    a%si(DXT)=SUM(a%f_xt*a%f_xt)/nxt
    a%si(DZT)=SUM(a%f_zt*a%f_zt)/nzt
    a%si(DYZT)=SUM(a%f_yzt*a%f_yzt)/nyzt
    a%si(DXZT)=SUM(a%f_xzt*a%f_xzt)/nxzt
    a%si(DYXT)=SUM(a%f_yxt*a%f_yxt)/nyxt
    a%si(DYXZ)=SUM(a%f_yxz*a%f_yxz)/nyxz

    a%si_residual = a%si(DTOT) & 
         - a%si(DX) - a%si(DY) -a%si(DZ) -a%si(DT) & 
         - a%si(DYX)- a%si(DYZ)-a%si(DYT)-a%si(DXZ)-a%si(DXT)-a%si(DZT) &
         - a%si(DYXZ)-a%si(DYXT)-a%si(DYZT)-a%si(DXZT) & 
         - a%si(DYXZT)

    IF ( a%si_residual .GT. SMALL_DELTA ) THEN 
       WRITE(*,*) 'ANOVA_DECOMP4D: RESIDUAL VARIANCE', a%si_residual, 'GLOBAL AVG:',a%f_empty 
       STOP 'ERROR: ANOVA_DECOMP4D VARIANCE NOT DECOMPOSED'
    ENDIF


    CALL TIMER_FINISH('ANOVA_DECOMP4D: STEP2 Time elapse') 

    RETURN 

  END SUBROUTINE ANOVA_DECOMP4D
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !
  SUBROUTINE CHECK_SLICES(vname,xdim,ydim,zdim,tdim) 

    IMPLICIT NONE 

    CHARACTER(LEN=*) :: vname 
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: data_yx,data_yz,data_yt,data_xz,data_xt,data_zt
    CHARACTER, DIMENSION(2) :: yx_dims*32,yz_dims*32,yt_dims*32,xz_dims*32,xt_dims*32,zt_dims*32 
    CHARACTER(LEN=*) :: xdim,ydim,zdim,tdim
    INTEGER,   DIMENSION(NCRW_MAXDIM) :: pos
    INTEGER :: nx,ny,nz,nt,t_start,t_finish
    INTEGER :: ix,iy,iz,it,max_dim 
    nx=ncrw_getdimlen(xdim) 
    ny=ncrw_getdimlen(ydim) 
    nz=ncrw_getdimlen(zdim) 
    nt=ncrw_getdimlen(tdim) 

    max_dim = MAXVAL( (/nx,ny,nz,nt/) )

    ALLOCATE(data_yx(ny,nx), data_yz(ny,nz),data_yt(ny,nt), data_xz(nx,nz), & 
         data_xt(nx,nt),data_zt(nz,nt))
    
    WRITE(*,*) 'NX=',nx,' NY=',ny,' NZ=',nz,' NT=',nt 

    yx_dims(1) = ydim;  yx_dims(2) = xdim 
    yz_dims(1) = ydim;  yz_dims(2) = zdim  
    yt_dims(1) = ydim;  yt_dims(2) = tdim 
    xz_dims(1) = xdim;  xz_dims(2) = zdim
    xt_dims(1) = xdim;  xt_dims(2) = tdim
    zt_dims(1) = zdim;  zt_dims(2) = tdim  
    
    pos(1) = 1 
    pos(2) = 1  

    CALL TIMER_START; DO iz=nz,nz-1,-1  
       pos(1) = iz 
       CALL ncrw_getvar_slice(vname,yx_dims,zt_dims,pos,data_yx)
    END DO; CALL TIMER_FINISH('yx_slice') 
    OPEN(17,FILE='yx_plane',ACCESS='STREAM',STATUS='UNKNOWN')
    WRITE(17) data_yx; CLOSE(17) 

    CALL TIMER_START; DO ix=nx,nx-1,-1
       pos(1) = ix
       CALL ncrw_getvar_slice(vname,yz_dims,xt_dims,pos,data_yz)
    ENDDO; CALL TIMER_FINISH('yz_slice')
    OPEN(17,FILE='yz_plane',ACCESS='STREAM',STATUS='UNKNOWN') 
    WRITE(17) data_yz; CLOSE(17)  
    
    CALL TIMER_START; DO iy=ny,ny-1,-1
       pos(1) = iy
       CALL ncrw_getvar_slice(vname,xz_dims,yt_dims,pos,data_xz)
    ENDDO; CALL TIMER_FINISH('xz_slice')
    OPEN(17,FILE='xz_plane',ACCESS='STREAM',STATUS='UNKNOWN')
    WRITE(17) data_xz; CLOSE(17) 

  END SUBROUTINE CHECK_SLICES 

END MODULE ANOVA
