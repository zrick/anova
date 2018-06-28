! Cedrick Ansorge, cedrick@posteo.de 15 Mar 2015 
! Module to parse a netCDF file 

MODULE nc_readwrite 

  USE netcdf 

  IMPLICIT NONE  
  SAVE 

  INTEGER, PARAMETER :: NCRW_MAXVAR=20, NCRW_MAXDIM=20,NCRW_MAXATT=20
  
  INTEGER :: ncrw_id  
  INTEGER :: ncrw_ndim, ncrw_nvar, ncrw_natt

  INTEGER,     DIMENSION(NCRW_MAXVAR)             :: ncrw_varids,ncrw_vartype,ncrw_varnatt,ncrw_varndim    
  INTEGER,     DIMENSION(NCRW_MAXVAR,NCRW_MAXDIM) :: ncrw_vardimids 
  LOGICAL,     DIMENSION(NCRW_MAXVAR)             :: ncrw_varisscaled,ncrw_varisdim  
  REAL(KIND=8),DIMENSION(NCRW_MAXVAR)             :: ncrw_varscale,ncrw_varoffset
  CHARACTER,   DIMENSION(NCRW_MAXVAR)             :: ncrw_varnames*32
  CHARACTER,   DIMENSION(NCRW_MAXVAR,NCRW_MAXDIM) :: ncrw_vardimnames*32  
  INTEGER,     DIMENSION(NCRW_MAXVAR,NCRW_MAXATT):: ncrw_varatttypes,ncrw_varattlen
  CHARACTER,   DIMENSION(NCRW_MAXVAR,NCRW_MAXATT):: ncrw_varattnames*32

  INTEGER,     DIMENSION(NCRW_MAXDIM)             :: ncrw_dimlen
  CHARACTER,   DIMENSION(NCRW_MAXDIM)             :: ncrw_dimnames*32

  LOGICAL :: ncrw_verbose

  REAL(KIND=8),DIMENSION(:,:), ALLOCATABLE        :: ncrw_grid

  PRIVATE check, initgrid, getvarid, getvar_slice_wrapper, getdimid

CONTAINS 
  
  SUBROUTINE ncrw_init(fname,verbose) 
    CHARACTER(LEN=*) fname
    INTEGER :: idim, ivar, iatt, idum
    INTEGER, OPTIONAL :: verbose 

    ncrw_verbose = .FALSE. 
    IF ( present(verbose) ) THEN  
       IF ( verbose .GT. 0 ) ncrw_verbose = .TRUE.  
    ENDIF
    
    
    CALL check ( nf90_open(fname,NF90_NOWRITE, ncrw_id) ) 
    CALL check ( nf90_inquire(ncrw_id,nDimensions=ncrw_ndim, &
         nVariables=ncrw_nvar,nAttributes=ncrw_natt) )

    IF ( ncrw_ndim .GT. NCRW_MAXDIM .OR. ncrw_nvar .GT. NCRW_MAXVAR ) THEN 
       STOP 'MODULE nc_read_write ERROR: number of dimensions or variables not supported'
    ENDIF

    IF ( ncrw_verbose ) WRITE(*,*) ncrw_ndim, 'Dimensions; ', ncrw_nvar, 'Variables; ', ncrw_natt, 'global Attributes;' 

    DO idim=1,ncrw_ndim  
       CALL check ( nf90_inquire_dimension(ncrw_id, idim, ncrw_dimnames(idim),ncrw_dimlen(idim) ) ) 
       IF ( ncrw_verbose ) THEN  
          WRITE(*,*) 'DIMENSION ', ncrw_dimnames(idim),' LEN=',ncrw_dimlen(idim) 
       ENDIF
    ENDDO
    IF ( ncrw_verbose ) WRITE(*,*) '================'
       

    DO ivar=1,ncrw_nvar 
       CALL check ( nf90_inquire_variable(ncrw_id, ivar, ncrw_varnames(ivar),ncrw_vartype(ivar),&
            ncrw_varndim(ivar),ncrw_vardimids(ivar,:), ncrw_varnatt(ivar) ) )
          
       IF ( ncrw_verbose ) THEN  
          WRITE(*,*) 'VARIABLE ', ncrw_varnames(ivar),'TYPE:',ncrw_vartype(ivar),'NDIM:',ncrw_varndim(ivar) ,&
            'NATT:', ncrw_varnatt(ivar) 
       ENDIF
       IF ( ncrw_varnatt(ivar) .GT. NCRW_MAXATT ) STOP 'MODULE nc_read_write ERROR: number of attributes too large' 
       
       DO idim=1,ncrw_varndim(ivar) 
          ncrw_vardimnames(ivar,idim)=ncrw_dimnames(ncrw_vardimids(ivar,idim)) 
       ENDDO
       IF ( ncrw_verbose ) WRITE(*,*) '   DIMs:',ncrw_vardimnames(ivar,:ncrw_varndim(ivar))
       ! IF ( ncrw_verbose ) WRITE(*,*) '        ',ncrw_vardimids(ivar,:ncrw_varndim(ivar))
       
       DO iatt=1,ncrw_varnatt(ivar) 
          CALL check ( nf90_inq_attname(ncrw_id,ivar,iatt,ncrw_varattnames(ivar,iatt) ) )
          CALL check ( nf90_inquire_attribute(ncrw_id,ivar,ncrw_varattnames(ivar,iatt), ncrw_varatttypes(ivar,iatt) ) ) 
          !CALL check ( nf90_get_att(ncrw_id, ivar, ncrw_varattnames(ivar,iatt) ) ) 
          IF (     ncrw_varattnames(ivar,iatt) .EQ. 'var_scale_factor' ) THEN  
             ncrw_varisscaled(ivar) = .TRUE.  
             CALL check ( nf90_get_att(ncrw_id,ivar,ncrw_varattnames(ivar,iatt),ncrw_varscale(ivar) ) )
          ELSEIF ( ncrw_varattnames(ivar,iatt) .EQ. 'var_add_offset' ) THEN 
             CALL check ( nf90_get_att(ncrw_id,ivar,ncrw_varattnames(ivar,iatt),ncrw_varoffset(ivar) ) )
          ELSEIF ( ncrw_varattnames(ivar,iatt) .EQ. 'axis' ) THEN 
             ncrw_varisdim(ivar) = .TRUE.  
          ENDIF
          ! IF ( ncrw_verbose ) WRITE(*,*) '   ATT:',iatt,ncrw_varattnames(ivar,iatt),ncrw_varatttypes(ivar,iatt) 
       ENDDO 
       
       IF ( ncrw_verbose ) THEN 
          IF ( ncrw_varisscaled(ivar)  ) THEN 
             WRITE(*,*) '   SCALE:',ncrw_varscale(ivar),ncrw_varoffset(ivar) 
          ELSE  
             WRITE(*,*) '   Variable is not scaled'
          ENDIF
          IF ( ncrw_varisdim(ivar) ) WRITE(*,*) '   variable is dimension data' 
       ENDIF
    ENDDO 
    
    CALL initgrid

    ! CHECK VARIABLES 
    DO ivar=1,ncrw_nvar 
       IF ( ivar .NE. getvarid(ncrw_varnames(ivar)) ) THEN 
          STOP 'ERROR INITIALIZING VARIABLES' 
       END IF
    ENDDO

    IF ( ncrw_verbose ) WRITE(*,*) '=========== FINISHED ncrw_init' 

  END SUBROUTINE ncrw_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GET LENGTH OF DIMENSION BY NAME 
  FUNCTION ncrw_getdimlen(n)  
    INTEGER :: ncrw_getdimlen 
    CHARACTER(len=*),INTENT(IN) :: n 
    ncrw_getdimlen=ncrw_dimlen(getdimid(n))
  END FUNCTION ncrw_getdimlen 

  FUNCTION getvardimpos(vid,n)
    INTEGER :: getvardimpos
    INTEGER,          INTENT(IN) :: vid
    CHARACTER(len=*), INTENT(IN) :: n 
    INTEGER :: idim,ndim,dum,dim_loc
    INTEGER,DIMENSION(NCRW_MAXDIM) :: dimids

    ndim=ncrw_varndim(vid)
    DO idim=1,ndim
       IF ( ncrw_vardimnames(vid,idim) .EQ. n ) THEN 
          getvardimpos=idim
          EXIT   
       ENDIF
    ENDDO
  END FUNCTION getvardimpos

  FUNCTION getdimid(n)   
    INTEGER :: getdimid
    CHARACTER(len=*), INTENT(IN) :: n 
    INTEGER :: idim
    
    DO idim=1,ncrw_ndim 
       IF ( ncrw_dimnames(idim) .EQ. n ) EXIT 
    ENDDO

    IF ( idim .GT. ncrw_ndim ) STOP 'MODULE nc_read_write ERROR: dimension not found'  

    getdimid = idim 
    RETURN 
  END FUNCTION getdimid
    
  FUNCTION ncrw_inq_dimtranspose(vname,d1,d2) 
    CHARACTER(LEN=*) :: d1, d2,vname 
    INTEGER :: dp1,dp2,vid  
    LOGICAL ncrw_inq_dimtranspose 

    vid=getvarid(vname) 
    dp1=getvardimpos(vid,d1) 
    dp2=getvardimpos(vid,d2) 
    
    IF ( dp1 .LT. dp2 ) THEN 
       ncrw_inq_dimtranspose=.FALSE. 
    ELSE 
       ncrw_inq_dimtranspose=.TRUE. 
    ENDIF

    RETURN
    
  END FUNCTION ncrw_inq_dimtranspose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ DATA SLICE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ncrw_getvar_slice(vname,slicedims,fixeddims,pos,values)  
    CHARACTER(LEN=*), INTENT(IN) :: vname   
    CHARACTER, DIMENSION(*)      :: slicedims*32 
    CHARACTER, DIMENSION(*)      :: fixeddims*32
    INTEGER, DIMENSION(*)        :: pos 
    REAL(KIND=8),DIMENSION(*)    :: values
    INTEGER                      :: vid  

    INTEGER :: vstart(ncrw_maxdim), vcount(ncrw_maxdim)
    INTEGER :: dummy_dimid,i
    INTEGER,DIMENSION(2) :: slice_size,slice_dimpos
    
    vid=getvarid(vname) 
    IF ( ncrw_varndim(vid) .LT. 2 ) STOP 'MODULE nc_read_write variables dimensions & 
         insufficient for read command'  

    DO i=1,2
       slice_dimpos(i) = getvardimpos(vid,slicedims(i)) 
       vstart(slice_dimpos(i)) = 1
       vcount(slice_dimpos(i)) = ncrw_dimlen(getdimid(slicedims(i)))
       slice_size(i)=ncrw_dimlen(slice_dimpos(i)) 
    ENDDO
 
    IF ( ncrw_inq_dimtranspose(vname,slicedims(1),slicedims(2)) ) THEN 
       STOP 'ERROR (MODULE_NC_READWRITE): TRANSPOSITION REQUIRED' 
    ENDIF

    DO i=1,ncrw_varndim(vid)-2  
       dummy_dimid = getvardimpos(vid,fixeddims(i))
       vstart(dummy_dimid) = pos(i) 
       vcount(dummy_dimid) = 1 
    ENDDO

    ! we need this wrapper as nf90_get_var requires the data array to have known dimensions 
    ! and we want to avoid having to pass the size of the slice to the public routine  
    CALL getvar_slice_wrapper(vid,vstart,vcount,ncrw_varndim(vid),slice_size(1),slice_size(2),values)

  END SUBROUTINE ncrw_getvar_slice 
  
  SUBROUTINE getvar_slice_wrapper(vid,vs,vc,ndim,n1,n2,v) 
    INTEGER,                       INTENT(IN)  :: vid,n1,n2,ndim
    INTEGER,      DIMENSION(ndim), INTENT(IN)  :: vs,vc 
    REAL(KIND=8), DIMENSION(n1,n2),INTENT(OUT) :: v
   
    !IF ( ncrw_verbose ) WRITE(*,*) 'IN getvar_slice_wrapper          ', ncrw_id,vid,' start:',vs,' count:',vc 
    !IF ( ncrw_verbose ) WRITE(*,*) 'IN getvar_slice_wrapper          ',vs,vc

    CALL check( nf90_get_var(ncrw_id,vid,v,start=vs,count=vc) )

    IF ( ncrw_varisscaled(vid) .EQV. .TRUE. ) THEN  
       v(:,:) = v(:,:)*ncrw_varscale(vid) - ncrw_varoffset(vid) 
    ENDIF
    
  END SUBROUTINE getvar_slice_wrapper


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ DATA LINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ncrw_getvar_line(vname)  
    CHARACTER(LEN=*), INTENT(IN) :: vname 
    INTEGER :: vid 

    vid = getvarid(vname) 

  END SUBROUTINE ncrw_getvar_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARID 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  FUNCTION getvarid(n)  
    INTEGER :: getvarid 
    CHARACTER(LEN=*), INTENT(IN) :: n  
    INTEGER :: i 
    DO i=1,ncrw_nvar 
       IF ( ncrw_varnames(i) .EQ. n ) EXIT 
    ENDDO
    

    IF ( i .GT. ncrw_nvar ) THEN  
       STOP 'MODULE nc_read_write ERROR: Variable not found'  
    ENDIF 
    getvarid=i 
    RETURN 
  END FUNCTION getvarid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GRID INITIALIZATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE initgrid 
    INTEGER :: ivar,dimlen,maxlen,dimid 
    INTEGER, DIMENSION(1) :: vstart,vcount
    maxlen = 0
    DO ivar=1,ncrw_nvar  
       IF ( ncrw_varisdim(ivar) ) THEN  
          dimid = ncrw_vardimids(ivar,1) 
          dimlen= ncrw_dimlen(dimid)
          IF ( dimlen .GT. maxlen ) maxlen = dimlen
       ENDIF
    ENDDO 
     
    ! The grid has zero fill values to acquire regular array shape
    ALLOCATE(ncrw_grid(ncrw_ndim,maxlen))  
    IF ( ncrw_verbose ) THEN 
       WRITE(*,*) '===============' 
       WRITE(*,*) 'INITIALIZE GRID'
    ENDIF
    DO ivar=1,ncrw_nvar  
       IF ( ncrw_varisdim(ivar) ) THEN   
          dimid = ncrw_vardimids(ivar,1) 
          dimlen= ncrw_dimlen(dimid)  
          IF ( ncrw_verbose ) WRITE(*,*) '   DATA FROM VAR ',ncrw_varnames(ivar), & 
               'FOR DIM ', ncrw_dimnames(dimid), 'N=',dimlen 
          vcount(1) = dimlen
          vstart(1) = 1 
          CALL check( nf90_get_var(ncrw_id,ivar,ncrw_grid(dimid,:),& 
               start=vstart,count=vcount ) )
       ENDIF
    ENDDO
  END SUBROUTINE initgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ERROR HANDLING 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE check(status)
    INTEGER, intent ( in) :: status
    
    IF(status /= nf90_noerr) then 
       PRINT *, trim(nf90_strerror(status))
       STOP "Stopped"
    ENDIF
  END SUBROUTINE check

END MODULE nc_readwrite
