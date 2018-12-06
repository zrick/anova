MODULE ANOVA  

  USE ANOVA_CONSTANTS
  USE NC_READWRITE 
  USE TIMER 

  PRIVATE :: ANOVA_GETDIMSTR

  TYPE T_ANOVA 
     INTEGER :: nx,ny,nz,nt,ndim
     INTEGER :: nt_sub,nx_sub 
     INTEGER :: ix_srt, ix_end, it_srt,it_end
     CHARACTER(LEN=64)                               :: xdim,ydim,zdim,tdim,vname,trange,xrange
     CHARACTER(LEN=64)                               :: fname 
     REAL(KIND=8)                                    :: f_empty,si_residual
     REAL(KIND=8),      DIMENSION(:),    ALLOCATABLE :: f_x,f_y,f_z,f_t 
     REAL(KIND=8),      DIMENSION(:,:),  ALLOCATABLE :: f_yx,f_yz,f_yt,f_xz,f_xt,f_zt
     REAL(KIND=8),      DIMENSION(:,:,:),ALLOCATABLE :: f_yxz,f_yxt,f_yzt,f_xzt 
     REAL(KIND=8),      DIMENSION(:),    ALLOCATABLE :: si 
     LOGICAL,           DIMENSION(:),    ALLOCATABLE :: f_save   
     INTEGER,           DIMENSION(:),    ALLOCATABLE :: si_order 
     CHARACTER(LEN=32), DIMENSION(:),    ALLOCATABLE :: si_tag
  END TYPE T_ANOVA

  CONTAINS

  SUBROUTINE GET_SORT_INDICES(n,arr,idx)
    
    IMPLICIT NONE 

    INTEGER,                      INTENT(IN) :: n 
    REAL(KIND=8),    DIMENSION(N),INTENT(IN) :: arr(n) 
    INTEGER,         DIMENSION(N),INTENT(OUT):: idx     
    !LOCAL 
    REAL a 
    INTEGER i,j 

    DO j=2, n
       a=arr(j)
       DO i=j-1,1,-1
          if (ARR(i)<=a) goto 10
          idx(i+1)=i
       ENDDO
       i=0
10     idx(i+1)=j
    ENDDO
    RETURN

  END SUBROUTINE GET_SORT_INDICES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET STRING ARRAY FOR DIMENSION NAMES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE ANOVA_SETDIMSTR(a) 
    
    IMPLICIT NONE

    INTEGER :: i  
    TYPE(T_ANOVA) :: a
    CHARACTER(LEN=20) :: x,y,z,t 
    x=a%xdim; y=a%ydim; z=a%zdim; t=a%tdim

    ALLOCATE(a%si_tag(DLAST)) 
    
    a%si_tag(:) =''

    WRITE(a%si_tag(DTOT),103) 'TOTAL VARIANCE' 
    WRITE(a%si_tag(DYXZ),105) y,x,z
    WRITE(a%si_tag(DYX),106)  y,x
    WRITE(a%si_tag(DYZ),106)  y,z
    WRITE(a%si_tag(DXZ),106)  x,z
    WRITE(a%si_tag(DY),107)   y
    WRITE(a%si_tag(DX),107)   x
    WRITE(a%si_tag(DZ),107)   z
    IF(a%ndim.GT.3) THEN 
       WRITE(a%si_tag(DYXZT),104)y,x,z,t 
       WRITE(a%si_tag(DYXT),105) y,x,t
       WRITE(a%si_tag(DYZT),105) y,z,t
       WRITE(a%si_tag(DXZT),105) x,z,t
       WRITE(a%si_tag(DYT),106)  y,t
       WRITE(a%si_tag(DXT),106)  x,t
       WRITE(a%si_tag(DZT),106)  z,t
       WRITE(a%si_tag(DT),107)   t 
    ENDIF 

103 FORMAT(A15)
104 FORMAT(A3,'-',A3,'-',A3,'-',A3,';')
105 FORMAT(A3,'-',A3,'-',A3,'    ;')
106 FORMAT(A3,'-',A3,'        ;')
107 FORMAT(A3,'            ;')


  END SUBROUTINE ANOVA_SETDIMSTR

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
    INTEGER :: nd,t0,t1,x0,x1
    CHARACTER(LEN=*) :: fbase
    CHARACTER(LEN=200) :: fname, fbase_loc

    xdim=a%xdim(:3); ydim=a%ydim(:3); zdim=a%zdim(:3); tdim=a%tdim(:3) 
    nd=a%ndim
    t0=a%it_srt
    t1=a%it_end
    x0=a%ix_srt 
    x1=a%ix_end
    
    WRITE(*,*) 'OUTPUT_BIN:', a%it_srt,a%it_end,a%nt_sub,a%nt,a%trange 

    fbase_loc = TRIM(fbase)  
    IF ( a%ix_srt .NE. 1 .OR. a%nx_sub .NE. a%nx) THEN 
       fbase_loc = TRIM(fbase_loc) // TRIM(a%xrange) 
    ENDIF
    IF ( a%it_srt .NE. 1 .OR. a%nt_sub .NE. a%nt) THEN 
       fbase_loc = TRIM(fbase_loc) // TRIM(a%trange) 
    ENDIF
    fbase_loc=TRIM(ADJUSTL(fbase_loc)) // '.' 

    WRITE(*,*) TRIM(fbase),'|',TRIM(fbase_loc),'|',TRIM(a%trange)

    IF ( wrt(DX) .EQV. .TRUE. ) THEN   
       fname=TRIM(fbase_loc)//TRIM(xdim)
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_x(x0:x1) 
       CLOSE(funit)
    ENDIF
    IF ( wrt(DY) .EQV. .TRUE. ) THEN   
       fname=TRIM(fbase_loc)//TRIM(ydim)
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_y(:) 
       CLOSE(funit) 
    ENDIF
    IF ( wrt(DZ) .EQV. .TRUE. ) THEN   
       fname=TRIM(fbase_loc)//TRIM(zdim)
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_z(:) 
       CLOSE(funit) 
    ENDIF
    IF ( wrt(DT) .EQV. .TRUE. .AND. nd .GT. 3) THEN   
       fname=TRIM(fbase_loc)//TRIM(tdim)
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_t(t0:t1) 
       CLOSE(funit) 
    ENDIF
    !
    ! 2-DIMENSIONAL FIELDS 
    !
    IF ( wrt(DYX) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase_loc)//TRIM(ydim)//'_'//TRIM(xdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yx(:,x0:x1) 
    ENDIF
    IF ( wrt(DYZ) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase_loc)//TRIM(ydim)//'_'//TRIM(zdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yz(:,:) 
    ENDIF
    IF ( wrt(DYT) .EQV. .TRUE. .AND. nd .GT. 3) THEN 
       fname=TRIM(fbase_loc)//TRIM(ydim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yt(:,t0:t1) 
    ENDIF
    IF ( wrt(DXZ) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase_loc)//TRIM(xdim)//'_'//TRIM(zdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_xz(x0:x1,:) 
    ENDIF
    IF ( wrt(DXT) .EQV. .TRUE. .AND.nd.GT. 3) THEN 
       fname=TRIM(fbase_loc)//TRIM(xdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_xt(x0:x1,t0:t1) 
    ENDIF
    IF ( wrt(DZT) .EQV. .TRUE. .AND.nd.GT.3) THEN 
       fname=TRIM(fbase_loc)//TRIM(zdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_zt(:,t0:t1) 
    ENDIF
    !
    ! 3-DIMENSIONAL FIELDS 
    !
    IF ( wrt(DYXZ) .EQV. .TRUE. ) THEN 
       fname=TRIM(fbase_loc)//TRIM(ydim)//'_'//TRIM(xdim)//'_'//TRIM(zdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yxz(:,x0:x1,:) 
    ENDIF
    IF ( wrt(DYXT) .EQV. .TRUE. .AND. nd.GT.3) THEN 
       fname=TRIM(fbase_loc)//TRIM(ydim)//'_'//TRIM(xdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yxt(:,x0:x1,t0:t1) 
    ENDIF
    IF ( wrt(DYZT) .EQV. .TRUE. .AND. nd.GT.3) THEN 
       fname=TRIM(fbase_loc)//TRIM(ydim)//'_'//TRIM(zdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_yzt(:,:,t0:t1) 
    ENDIF
    IF ( wrt(DXZT) .EQV. .TRUE. .AND. nd.GT.3) THEN 
       fname=TRIM(fbase_loc)//TRIM(xdim)//'_'//TRIM(zdim)//'_'//TRIM(tdim) 
       OPEN(funit,FILE=fname,ACCESS='STREAM',FORM='UNFORMATTED')  
       WRITE(funit) a%f_xzt(x0:x1,:,t0:t0) 
    ENDIF
    
    
  END SUBROUTINE ANOVA_OUTPUT_BIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ANOVA ASCII OUTPUT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE ANOVA_OUTPUT_ASC(a,fbase)    
    
    IMPLICIT NONE  

    ! PARAMETERS 
    TYPE(T_ANOVA), INTENT(IN):: a  
    CHARACTER(LEN=*)         :: fbase
    ! LOCAL DEFINITIONS 
    CHARACTER(LEN=200)       :: fname  
    INTEGER,       PARAMETER :: funit=27 
    INTEGER                  :: i
    REAL(KIND=8)             :: si_max 


    fname = TRIM(ADJUSTL(fbase))
    IF ( a%ix_srt .NE. 1 .OR. a%nx_sub .NE. a%nx ) THEN 
       fname = TRIM(ADJUSTL(fname)) // TRIM(ADJUSTL(a%xrange))
    ENDIF
    IF ( a%it_srt .NE. 1 .OR. a%nt_sub .NE. a%nt ) THEN 
       fname = TRIM(ADJUSTL(fname)) // TRIM(ADJUSTL(a%trange))
    ENDIF

    IF ( a%ndim.LT.3 .OR. a%ndim .GT. 4 ) & 
       STOP 'ANOVA NOT IMPLEMENTED FOR LESS THAN 3 OR MORE THAN 4 DIMENSIONS' 

    si_max=MAXVAL(a%si(2:DLAST))

    OPEN(funit,FILE=fname,form='FORMATTED') 

    WRITE(funit,104) 'TOTAL VAR',   a%si(Dtot),1.0 
    WRITE(funit,104) 'RESIDUAL VAR',a%si_residual,a%si_residual/a%si(DTOT)  
    WRITE(funit,106) 'DIM', 'VARIANCE', 'SENSITIVITY','SENSITIVITY-2'

    DO i=2,DLAST  
       IF(len(TRIM(a%si_tag(i))) .GT. 0) & 
            WRITE(funit,105) a%si_tag(i),a%si(i),&
            100.*a%si(i)/a%si(Dtot),100.*a%si(i)/(a%si(Dtot)-si_max)
    ENDDO
  ! 
    CLOSE(funit) 
106 FORMAT('#',a14,';',a10,  ';',a12,';',a14)
105 FORMAT(a15,';',g11.4,';',g12.4'%;',g12.4,'%;')
104 FORMAT('#',a14,';',g11.4,';',g12.4';')
  END SUBROUTINE ANOVA_OUTPUT_ASC




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ANOVA ASCII OUTPUT FOR PROFILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE ANOVA_OUTPUT_ASCP(a,fbase,p)    
    
    IMPLICIT NONE  

    ! PARAMETERS 
    TYPE(T_ANOVA),                    INTENT(IN) :: a  
    REAL(KIND=8),  DIMENSION(DLAST,*),INTENT(IN) :: p 
    CHARACTER(LEN=*),                 INTENT(IN) :: fbase

    ! LOCAL DEFINITIONS 
    INTEGER,       PARAMETER :: funit=27 
    INTEGER                  :: i
    REAL(KIND=8)             :: si_max 
    CHARACTER(LEN=200)       :: fname 

    fname = TRIM(ADJUSTL(fbase)) // TRIM(ADJUSTL(a%xrange)) //TRIM(ADJUSTL(a%trange))
 
    IF ( a%ndim.LT.3 .OR. a%ndim .GT. 4 ) & 
       STOP 'ANOVA NOT IMPLEMENTED FOR LESS THAN 3 OR MORE THAN 4 DIMENSIONS' 

    si_max=MAXVAL(a%si(2:DLAST))

    OPEN(funit,FILE=fname,form='FORMATTED') 

    WRITE(funit,108) a%si_tag((/1,3,7,8,10,13,14,15/))

    DO i=a%it_srt,a%it_end
       WRITE(funit,107) i,p( (/1,3,7,8,10,13,14,15 /), i)
    ENDDO
  ! 
    CLOSE(funit) 
108 Format('LEVEL;',8(A11,  ';'))
107 FORMAT(I5';',   8(G11.4,';')) 
  END SUBROUTINE ANOVA_OUTPUT_ASCP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ANOVA OUTPUT FOR DATA COMPRESSION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  SUBROUTINE ANOVA_OUTPUT_COMPRESS(a,fname)  
    
    IMPLICIT NONE 

    ! PARAMETERS 
    TYPE(T_ANOVA), INTENT(IN):: a  

    ! LOCAL DEFINITIONS 
    CHARACTER(LEN=*)         :: fname  
    INTEGER,       PARAMETER :: funit=27 
    INTEGER                  :: ndim
    REAL(KIND=8)             :: si_max 


    
    
  END SUBROUTINE ANOVA_OUTPUT_COMPRESS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ANOVA 3D - variance decomposition and sensitivity indices 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ANOVA_DECOMP3D(vname,fixed_pos,a)

    USE NC_READWRITE, ONLY : NCRW_MAXDIMLEN

    IMPLICIT NONE 

    INTEGER,          INTENT(IN) :: fixed_pos
    CHARACTER(LEN=*), INTENT(IN) :: vname
    TYPE(T_ANOVA)                :: a

    CHARACTER(LEN=NCRW_MAXDIMLEN)             :: xdim,ydim,zdim,fixed_dim

    REAL(KIND=8)                              :: avg_yxz,rdum               !,f_empty
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: avg_yx, avg_yz,avg_xz      !,f_x,f_y,f_z  
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: avg_x,avg_y,avg_z,data,wrk !,f_yx,f_yz,f_xz   
    ! REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::                            !f_yxz
    INTEGER,      DIMENSION(2)                :: pos 
    CHARACTER,    DIMENSION(2)                :: yz_dims*32,xt_dims*32  


    INTEGER :: ix,iy,iz
    INTEGER :: nx,ny,nz,nt,nyz,nyx,nxz,nyxz  
    INTEGER :: x0,x1
    REAL :: chk_x_vs_y, chk_x_vs_z, chk_x_vs_yx, chk_x_vs_xz, chk_x_vs_yz
    REAL :: chk_y_vs_z,             chk_y_vs_yx, chk_y_vs_xz, chk_y_vs_yz
    REAL ::                         chk_z_vs_yx, chk_z_vs_xz, chk_z_vs_yz 
    REAL :: chk_yx_vs_yz, chk_yx_vs_xz, chk_yz_vs_xz 
    REAL :: max_err 

    IF ( ncrw_verbose ) CALL TIMER_START()

    xdim=a%xdim
    ydim=a%ydim
    zdim=a%zdim
    fixed_dim=a%tdim

    nx=a%nx
    ny=a%ny
    nz=a%nz
    nt=a%nt

    IF ( ncrw_verbose ) THEN 
       WRITE(*,*) '============'
       WRITE(*,*) 'ANOVA_DECOMP3D: VARIABLE ',TRIM(vname),' NX=',nx,' NY=',ny,' NZ=',nz, ' FIXED_POS=',fixed_pos 
    ENDIF
    yz_dims(1) = ydim;  yz_dims(2) = zdim
    xt_dims(1) = xdim;  xt_dims(2) = fixed_dim

    ALLOCATE(avg_yx(nz),avg_yz(nx),avg_xz(ny)) 
    ALLOCATE(avg_x(ny,nz), avg_y(nx,nz), avg_z(ny,nx),data(ny,nz),wrk(ny,nz)) 

    nx=a%nx_sub
    x0=a%ix_srt 
    x1=a%ix_end  

    nyz=ny*nz 
    nxz=nx*nz 
    nyx=ny*nx
    nyxz=ny*nx*nz

    IF ( ncrw_verbose ) THEN 
       CALL TIMER_FINISH('ANOVA_DECOMP3D: init time elapse')
       CALL TIMER_START() 
       WRITE(*,*) 'ANOVA_DECOMP3D: STEP 1 - INTEGRATION' 
    ENDIF 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 1: INTEGRATION 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    avg_yxz=0. 
    avg_yx(:)=0.; avg_yz(:)=0.; avg_xz(:)=0.
    avg_x(:,:)=0.;avg_y(:,:)=0.;avg_z(:,:)=0.
    pos(2) = fixed_pos
    ! 
    DO ix=x0,x1
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
    ENDIF
    
    IF ( ncrw_verbose ) THEN  
       WRITE(*,*) 'ANOVA_DECOMP3D: STEP 1 Finished, Integral Calculus consistent'  
       CALL TIMER_FINISH('ANOVA_DECOMP3D: STEP 1 Time elapse') 
       CALL TIMER_START() 
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 2: DECOMPOSITION 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    IF ( ncrw_verbose ) WRITE(*,*) 'ANOVA_DECOMP3D: STEP 2 - DECOMPOSITION AND SENSITIVITY INDICES' 
    a%f_empty=avg_yxz 
    ! 
    a%f_x=avg_yz-a%f_empty
    a%f_y=avg_xz-a%f_empty 
    a%f_z=avg_yx-a%f_empty  
    ! 
    DO iy=1,ny 
       DO ix=x0,x1; a%f_yx(iy,ix)=avg_z(iy,ix)-a%f_x(ix)-a%f_y(iy)-a%f_empty;ENDDO 
       DO iz=1,nz;  a%f_yz(iy,iz)=avg_x(iy,iz)-a%f_z(iz)-a%f_y(iy)-a%f_empty;ENDDO
    ENDDO
    !
    DO ix=x0,x1
       DO iz=1,nz;  a%f_xz(ix,iz)=avg_y(ix,iz)-a%f_z(iz)-a%f_x(ix)-a%f_empty;ENDDO  
    ENDDO 

    a%si(:DLAST) = 0. 
    DO ix=x0,x1
       pos(1)=ix
       CALL ncrw_getvar_slice(vname,yz_dims,xt_dims,pos,data)   
       DO iy=1,ny; DO iz=1,nz  
          a%f_yxz(iy,ix,iz) = data(iy,iz) - a%f_yx(iy,ix) - a%f_xz(ix,iz) - a%f_yz(iy,iz) & 
               - a%f_x(ix) - a%f_y(iy) - a%f_z(iz) - a%f_empty
       ENDDO; ENDDO 
       data = data(:,:) - a%f_empty  
       a%si(DTOT) = a%si(DTOT) + ( SUM(data * data) ) / nyxz
       a%si(DYXZ) = a%si(DYXZ) + ( SUM(a%f_yxz(:,ix,:)*a%f_yxz(:,ix,:)) ) / nyxz
    ENDDO

    a%si(DX)  = SUM(a%f_x(x0:x1)*a%f_x(x0:x1))/nx 
    a%si(DY)  = SUM(a%f_y*a%f_y)/ny
    a%si(DZ)  = SUM(a%f_z*a%f_z)/nz
    a%si(DYX) = SUM(a%f_yx(:,x0:x1)*a%f_yx(:,x0:x1))/nyx
    a%si(DYZ) = SUM(a%f_yz*a%f_yz)/nyz 
    a%si(DXZ) = SUM(a%f_xz(x0:x1,:)*a%f_xz(x0:x1,:))/nxz  
    a%si_residual = a%si(DTOT) - a%si(DYXZ)-a%si(DYX)-a%si(DYZ)-a%si(DXZ)-a%si(DX)-a%si(DY)-a%si(DZ) 

    IF ( ABS ( a%si_residual/a%si(DTOT)  ) .GT. 1e-6 ) THEN 
       WRITE(STDOUT,*) 'ANOVA_DECOMP3D: VARIANCE DECOMPOSITION FAILED, RESIDUAL:', a%si_residual
       STOP 'ANOVA_DECOMP3D: VARIANCE DECOMPOSITION FAILED' 
    ELSEIF ( ncrw_verbose ) THEN 
       WRITE(STDOUT,*) 'ANOVA_DECOMP3D: VARIANCE DECOMPOSED TO WITHIN', ABS(a%si_residual) / a%si(DTOT)  
    ENDIF  
104 FORMAT(A30,1x,G10.3) 
    !  
    ! Check orthogonality 
    ! xy 
    chk_x_vs_y=0.; chk_x_vs_z=0.; chk_x_vs_yx=0.; chk_x_vs_xz=0.; chk_x_vs_yz=0.
    chk_y_vs_z=0.;                chk_y_vs_yx=0.; chk_y_vs_xz=0.; chk_y_vs_yz=0.
                                  chk_z_vs_yx=0.; chk_z_vs_xz=0.; chk_z_vs_yz=0.
    chk_yx_vs_yz=0.; chk_yx_vs_xz=0.; chk_yz_vs_xz=0.; 
    DO ix=x0,x1 
       DO iy=1,ny; 
          chk_x_vs_y = chk_x_vs_y  +a%f_x(ix)*a%f_y(iy); 
          chk_x_vs_yx = chk_x_vs_yx+a%f_x(ix)*a%f_yx(iy,ix)  
          chk_y_vs_yx = chk_y_vs_yx+a%f_y(iy)*a%f_yx(iy,ix)
          DO iz=1,nz  
             chk_z_vs_yx = chk_z_vs_yx + a%f_z(iz)*a%f_yx(iy,ix) 
             chk_x_vs_yz = chk_x_vs_yz + a%f_x(ix)*a%f_yz(iy,iz) 
             chk_y_vs_xz = chk_y_vs_xz + a%f_y(iy)*a%f_xz(ix,iz) 
             chk_yx_vs_yz= chk_yx_vs_yz+ a%f_yx(iy,ix)*a%f_yz(iy,iz)  
             chk_yx_vs_xz= chk_yx_vs_xz+ a%f_yx(iy,ix)*a%f_xz(ix,iz) 
             chk_yz_vs_xz= chk_yz_vs_xz+ a%f_yz(iy,iz)*a%f_xz(ix,iz) 
          ENDDO
       ENDDO
       DO iz=1,nz
          chk_x_vs_z = chk_x_vs_z + a%f_x(ix)*a%f_z(iz); 
          chk_x_vs_xz= chk_x_vs_xz+ a%f_x(ix)*a%f_xz(ix,iz); 
          chk_z_vs_xz= chk_z_vs_xz+ a%f_z(iz)*a%f_xz(ix,iz); 
       ENDDO
    ENDDO
    DO iy=1,ny 
       DO iz=1,nz 
          chk_y_vs_z = chk_y_vs_z + a%f_y(iy)*a%f_z(iz) 
          chk_y_vs_yz= chk_y_vs_yz+ a%f_y(iy)*a%f_yz(iy,iz) 
          chk_z_vs_yz= chk_z_vs_yz+ a%f_z(iz)*a%f_yz(iy,iz) 
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

    IF ( max_err /a%si(DTOT)  .GT. SQRT(SMALL_DELTA) ) THEN 
       WRITE(*,*) 'ANOVA_DECOMP3D: Orthogonality of ANOVA decomposition violated by', max_err /a%si(DTOT)
       STOP 'ANOVA_DECOMP3D: ORTHOGONALITY VIOLATED'
    ENDIF

    IF ( ncrw_verbose ) THEN 
       WRITE(*,*) 'ANOVA_DECOMP3D: GLOBAL AVERAGE:', a%f_empty 
!       WRITE(*,*) 'YZ:',SUM(a%f_yz*a%f_yz)/nyz, MINVAL(a%f_yz),MAXVAL(a%f_yz), SUM(a%f_yz)/nyz,nyz
!       WRITE(*,*) 'YX:',SUM(a%f_yx*a%f_yx)/nyx, MINVAL(a%f_yx),MAXVAL(a%f_yx), SUM(a%f_yx)/nyx,nyx
!       WRITE(*,*) 'XZ:',SUM(a%f_xz*a%f_xz)/nxz, MINVAL(a%f_xz),MAXVAL(a%f_xz), SUM(a%f_xz)/nxz,nxz
!       WRITE(*,*) 'X: ',SUM(a%f_x*a%f_x)/nx, MINVAL(a%f_x),MAXVAL(a%f_x),      SUM(a%f_x)/nx,  nx
!       WRITE(*,*) 'Y: ',SUM(a%f_y*a%f_y)/ny, MINVAL(a%f_y),MAXVAL(a%f_y),      SUM(a%f_y)/ny,  ny
!       WRITE(*,*) 'Z: ',SUM(a%f_z*a%f_z)/nz, MINVAL(a%f_z),MAXVAL(a%f_z),      SUM(a%f_z)/nz,  nz  
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

    CALL GET_SORT_INDICES(DLAST,a%si,a%si_order)

    IF ( ncrw_verbose ) & 
         CALL TIMER_FINISH('ANOVA_DECOMP3D: STEP2 Time elapse') 

  END SUBROUTINE ANOVA_DECOMP3D   


  SUBROUTINE ANOVA_INIT(xdim,ydim,zdim,tdim,a,ndim,ix_srt,ix_end,it_srt,it_end)   
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: xdim,ydim,zdim,tdim   
    INTEGER,          INTENT(IN) :: ndim
    INTEGER,          INTENT(IN), OPTIONAL :: ix_srt,ix_end,it_srt, it_end
    INTEGER :: nx,ny,nz,nt 
    TYPE(T_ANOVA) :: a 

    WRITE(*,*) 'INITIALIZING', ndim,'-dim ANOVA'
    a%ndim=ndim 
    a%nx=ncrw_getdimlen(xdim);  a%xdim=xdim; nx=a%nx
    a%ny=ncrw_getdimlen(ydim);  a%ydim=ydim; ny=a%ny 
    a%nz=ncrw_getdimlen(zdim);  a%zdim=zdim; nz=a%nz 
    a%nt=ncrw_getdimlen(tdim);  a%tdim=tdim; nt=a%nt  
    !
    ! Parse inpupt for subset in last index 
    a%ix_srt = 1 
    a%ix_end = nx 
    IF ( PRESENT(ix_srt) ) THEN 
       IF ( ix_srt .GT. 0 ) a%ix_srt=ix_srt 
    ENDIF
    IF ( PRESENT(ix_end) ) THEN 
       IF ( ix_end .GT. 0 ) a%ix_end=ix_end  
    ENDIF

    a%nx_sub = (a%ix_end - a%ix_srt) + 1 
    IF ( a%nx_sub .GT. nx .OR. a%ix_end .GT. nx ) THEN 
       WRITE(*,*) 'ERROR (ANOVA_INIT): index bounds given for exceed range of data'
       WRITE(*,*) '                        ix_srt',a%ix_srt,' a%ix_end',a%ix_end, 'a%nx_sub', a%nx_sub     
       STOP 'ERROR' 
    ENDIF
    !
    IF ( a%ix_srt .NE. 1 .OR. a%nx_sub .NE. a%nx ) THEN  
       WRITE(a%xrange,109) a%ix_srt,a%ix_end  ! Format defined below 
       a%xrange = '.'//TRIM(ADJUSTL(a%xdim))//TRIM(ADJUSTL(a%xrange))
    ELSE 
       a%trange = '' 
    ENDIF


    a%it_srt = 1 
    a%it_end = nt 
    IF ( PRESENT(it_srt) ) THEN 
       IF ( it_srt .GT. 0 ) a%it_srt=it_srt
    ENDIF
    ! 
    IF ( PRESENT(it_end) ) THEN 
       IF ( it_end .GT. 0 ) a%it_end=it_end
    ENDIF
    a%nt_sub = (a%it_end - a%it_srt) + 1 
    IF ( a%nt_sub .GT. nt .OR. a%it_end .GT. nt ) THEN 
       WRITE(*,*) 'ERROR (ANOVA_INIT): index bounds given for exceed range of data'
       WRITE(*,*) '                        it_srt',a%it_srt,' a%it_end',a%it_end, 'a%nt_sub', a%nt_sub     
       STOP 'ERROR' 
    ENDIF
    !
     IF ( a%it_srt .NE. 1 .OR. a%nt_sub .NE. a%nt ) THEN  
        WRITE(a%trange,109) a%it_srt,a%it_end 
        a%trange = '.'//TRIM(ADJUSTL(a%tdim))//TRIM(ADJUSTL(a%trange))
109     FORMAT(I4.4,'-',I4.4)
     ELSE 
        a%trange = '' 
     ENDIF

    WRITE(*,*) 'ANOVA_INIT: processing x from', a%ix_srt, ' to ', a%ix_end, '(', a%nx_sub, ')'
    WRITE(*,*) 'ANOVA_INIT: processing t from', a%it_srt, ' to ', a%it_end, '(', a%nt_sub, ')'
    !
    CALL ANOVA_SETDIMSTR(a) 
    ALLOCATE(a%si_order(DLAST)) 
    ! 
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

  SUBROUTINE ANOVA_DECOMP4D(vname,si_res,a)

    USE NC_READWRITE, ONLY : NCRW_MAXDIMLEN

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)           :: vname
    REAL(KIND=8)  :: si_res   
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

    CHARACTER(LEN=NCRW_MAXDIMLEN) :: xdim,ydim,zdim,tdim

    INTEGER :: ix,iy,iz,it   
    INTEGER :: nx,ny,nz,nt
    INTEGER :: nyz,nyx,nyt,nxz,nxt,nzt
    INTEGER :: nyzt,nxzt,nyxt,nyxz 
    INTEGER :: nyxzt 
    INTEGER :: t0,t1,x0,x1
    CALL TIMER_START()

    xdim=a%xdim
    ydim=a%ydim
    zdim=a%zdim
    tdim=a%tdim

    nx=a%nx !ncrw_getdimlen(xdim);
    ny=a%ny !ncrw_getdimlen(ydim);
    nz=a%nz !ncrw_getdimlen(zdim); 
    nt=a%nt !ncrw_getdimlen(tdim); 
    ! 
    t0=a%it_srt 
    t1=a%it_end 
    x0=a%ix_srt 
    x1=a%ix_end
    
    IF ( ncrw_verbose ) THEN 
       WRITE(*,*) '=============='
       WRITE(*,*) 'ANOVA_DECOMP4D: VARIABLE ',TRIM(vname),' NX=',a%nx,' NY=',a%ny,' NZ=',a%nz, ' NT=',a%nt  
       WRITE(*,*) '                PROCESSING RANGE IN LAST INDEX:', t0,t1,nt
       WRITE(*,*) '                PROCESSING RANGE IN FIRST INDEX:',x0,x1,nx
    ENDIF

    !
    yz_dims(1) = ydim;  yz_dims(2) = zdim
    xt_dims(1) = xdim;  xt_dims(2) = tdim
    !
    ALLOCATE(avg_yxz(nt),avg_yzt(nx),avg_xzt(ny),avg_yxt(nz)) 
    ALLOCATE(avg_yx(nz,nt),avg_yz(nx,nt),avg_yt(nx,nz),avg_xz(ny,nt),avg_xt(ny,nz),avg_zt(ny,nx))  
    ALLOCATE(avg_x(ny,nz,nt), avg_y(nx,nz,nt), avg_z(ny,nx,nt),avg_t(ny,nx,nz)) 
    ALLOCATE(data(ny,nz),f_yxzt(ny,nz)) 

    nt=a%nt_sub   ! This hack allocates the full arrays, 
    nx=a%nx_sub   ! }but normalizes only by the sub-sampled data

    nyz=ny*nz;     nyx=ny*nx;     nyt=ny*nt;     nxz=nx*nz;     nxt=nx*nt; nzt=nz*nt
    nyzt=ny*nz*nt; nxzt=nx*nz*nt; nyxt=ny*nx*nt; nyxz=ny*nx*nz
    nyxzt=ny*nx*nz*nt



    IF ( ncrw_verbose ) & 
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
    ! 
    DO it=t0,t1;  DO ix=x0,x1 
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
    
    IF ( ncrw_verbose ) & 
         CALL TIMER_FINISH('ANOVA_DECOMP4D: STEP 1 time elapse') 
    CALL TIMER_START() 


    IF ( ncrw_verbose ) & 
         WRITE(*,*) 'ANOVA_DECOMP4D: STEP 2 -- DECOMPOSITION' 
    
    a%f_empty = avg_yxzt 
    a%f_x = avg_yzt - a%f_empty 
    a%f_y = avg_xzt - a%f_empty 
    a%f_z = avg_yxt - a%f_empty 
    a%f_t = avg_yxz - a%f_empty
    
    DO iy=1,ny
       DO ix=x0,x1; a%f_yx(iy,ix)=avg_zt(iy,ix)-a%f_y(iy)-a%f_x(ix)-a%f_empty; ENDDO 
       DO iz=1,nz;  a%f_yz(iy,iz)=avg_xt(iy,iz)-a%f_y(iy)-a%f_z(iz)-a%f_empty; ENDDO 
       DO it=t0,t1; a%f_yt(iy,it)=avg_xz(iy,it)-a%f_y(iy)-a%f_t(it)-a%f_empty; ENDDO 
    ENDDO 
    DO ix=x0,x1
       DO iz=1,nz;  a%f_xz(ix,iz)=avg_yt(ix,iz)-a%f_x(ix)-a%f_z(iz)-a%f_empty; ENDDO 
       DO it=t0,t1; a%f_xt(ix,it)=avg_yz(ix,it)-a%f_x(ix)-a%f_t(it)-a%f_empty; ENDDO 
    ENDDO 

    DO iz=1,nz
       DO it=t0,t1; a%f_zt(iz,it)=avg_yx(iz,it)-a%f_z(iz)-a%f_t(it)-a%f_empty; ENDDO 
    ENDDO

    DO iy=1,ny
       DO ix=x0,x1; DO iz=1,nz
          a%f_yxz(iy,ix,iz)=avg_t(iy,ix,iz)-a%f_y(iy)-a%f_x(ix)-a%f_z(iz)-a%f_yx(iy,ix)-a%f_yz(iy,iz)-a%f_xz(ix,iz)-a%f_empty
       ENDDO; ENDDO 
       DO ix=x0,x1; DO it=t0,t1
          a%f_yxt(iy,ix,it)=avg_z(iy,ix,it)-a%f_y(iy)-a%f_x(ix)-a%f_t(it)-a%f_yx(iy,ix)-a%f_yt(iy,it)-a%f_xt(ix,it)-a%f_empty
       ENDDO; ENDDO 
       DO iz=1,nz; DO it=t0,t1
          a%f_yzt(iy,iz,it)=avg_x(iy,iz,it)-a%f_y(iy)-a%f_z(iz)-a%f_t(it)-a%f_yz(iy,iz)-a%f_yt(iy,it)-a%f_zt(iz,it)-a%f_empty  
       ENDDO;ENDDO
    ENDDO

    DO ix=x0,x1     
       DO iz=1,nz; DO it=t0,t1
          a%f_xzt(ix,iz,it)=avg_y(ix,iz,it)-a%f_x(ix)-a%f_z(iz)-a%f_t(it)-a%f_xz(ix,iz)-a%f_xt(ix,it)-a%f_zt(iz,it)-a%f_empty  
       ENDDO;ENDDO 
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! STEP 3 Total variance and residual 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    a%si(:DLAST) = 0. 
    DO it=t0,t1; DO ix=x0,x1  
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

    IF ( a%si_residual/a%si(DTOT) .GT. SMALL_DELTA ) THEN 
       WRITE(*,*) 'ANOVA_DECOMP4D: RESIDUAL VARIANCE', a%si_residual, 'GLOBAL AVG:',a%f_empty 
       STOP 'ERROR: ANOVA_DECOMP4D VARIANCE NOT DECOMPOSED'
    ENDIF

    CALL GET_SORT_INDICES(DLAST,a%si,a%si_order)

    IF ( ncrw_verbose ) & 
         CALL TIMER_FINISH('ANOVA_DECOMP4D: STEP 2 Time elapse') 

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
