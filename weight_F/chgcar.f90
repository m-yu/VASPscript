MODULE chgcar
  USE prec
  USE options
  USE matrix

  IMPLICIT NONE
!-----------------------------------------------------------------------------!
! Variable description:
! ntypes : number of types
! nitype : number of ions for each type
! nions : number of ions
! ionpos_dir : ionic positions in direct coords.
! ionpos_lat : ionic positions in lattice coords.
! dir2car : matrix converting direct corrds into cartesian coords
! lat2car : matrix converting lattice corrds into cartesian coords
! np(3) : 3D mesh grids 
! nrho : total grid points
! rho : charge density on grids
! Nneigh : number of Wigner-Seitz neighbors
! nvert : index of Wigner-Seitz neighbors 
! alpha : surface area between WS neighbors 
!-----------------------------------------------------------------------------!
  TYPE charge
    CHARACTER(40) :: SZNAM
    REAL(q) :: SCALE, A(3,3)
    CHARACTER :: mode
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: nitype 
    INTEGER :: ntypes,nions,np(3),nrho,Nneigh
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: ionpos_dir,ionpos_lat  
    REAL(q),ALLOCATABLE,DIMENSION(:) :: rho
    REAL(q),DIMENSION(3,3) :: dir2car,lat2car 
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: nvert  
    REAL(q),ALLOCATABLE,DIMENSION(:) :: alpha
    REAL(q) :: dminsq
  END TYPE

  PUBLIC :: charge
  PUBLIC :: read_chgcar,write_chgcar
  PUBLIC :: pbc1D,pbc3D,pbc_dcar

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE assign_charge
  END INTERFACE

  CONTAINS

  SUBROUTINE read_chgcar(chgcarfile,chg)

  TYPE(charge) :: chg

  CHARACTER(40) :: chgcarfile
  INTEGER :: IU,I,nt,k
  INTEGER,DIMENSION(20) :: iontype=0
 
  IU =1000

  WRITE(*,'(2X,1A10,1A20)') 'OPEN ... ',chgcarfile
  OPEN(unit=IU,FILE=chgcarfile,STATUS='old',ACTION='read')

! the 'name' of the system
    READ(IU,'(A40)')chg%SZNAM

! scaling factor
    READ(IU,*)chg%SCALE

! three lattice vectors
    READ(IU,*) (chg%A(I,1:3) , I=1,3)  

! the number of atoms per atomic species
    READ(IU,'(20I4)') iontype
    DO nt=1,20
      IF(iontype(nt).eq.0) GOTO 100
    ENDDO
100 continue
    chg%ntypes=nt-1
    ALLOCATE(chg%nitype(chg%ntypes))
    DO nt=1,chg%ntypes
      chg%nitype(nt)=iontype(nt)
    END DO
    chg%nions=SUM(chg%nitype)

! direct coordinates mode
    READ(IU,'(A)') chg%mode
    IF (chg%mode .ne. 'd' .AND. chg%mode .ne. 'D') THEN
        WRITE (*,*) 'ERROR: READ IN chargefile not Direct mode'
        GOTO 200
    ENDIF

! direct ionic positions
    ALLOCATE(chg%ionpos_dir(chg%nions,3),chg%ionpos_lat(chg%nions,3))
    READ(IU,'(3F10.6)')(chg%ionpos_dir(I,1:3), I=1,chg%nions)

! grid meshes
    READ(IU,*)
    READ(IU,*)chg%np(1:3)
    chg%nrho=chg%np(1)*chg%np(2)*chg%np(3)

! charge density
    ALLOCATE(chg%rho(chg%nrho))
    READ(IU,*) (chg%rho(k), k=1,chg%nrho)
    CLOSE(IU)
    WRITE(*,'(2X,1A10,1A20)') 'CLOSE ... ',chgcarfile

    CALL matrix_transpose(chg%SCALE*chg%A, chg%dir2car) 
    DO I =1,3
      chg%ionpos_lat(:, I) = chg%ionpos_dir(:,I)*chg%np(I)+1._q
      chg%lat2car(:,I) = chg%dir2car(:,I)/REAL(chg%np(I),q) 
    ENDDO
200 continue
  END SUBROUTINE read_chgcar

!-----------------------------------------------------------------------------!
! Write the charge density from a file in chgcar format
!-----------------------------------------------------------------------------!
    
  SUBROUTINE write_chgcar(writeoutfile,chg,dens)

    TYPE(charge) :: chg
    CHARACTER(40) :: writeoutfile
    REAL(q),INTENT(IN) :: dens(chg%nrho)
    INTEGER :: IU,I,nt,k
    INTEGER :: NWRITE, NWRITTEN
    CHARACTER(40) :: FORM

    IU =1000

   WRITE(*,'(1A22,1A20)') 'WRITE OUT DENSITY ... ',writeoutfile
   OPEN(unit=IU,FILE=writeoutfile)
      WRITE(IU,'(A40)')chg%SZNAM
      WRITE(IU,*)chg%SCALE
      WRITE(IU,'(1X,3F12.6)')(chg%A(I,1:3),I=1,3)
      WRITE(IU,'(20I4)')(chg%nitype(nt),nt=1,chg%ntypes)
      WRITE(IU,'(A6)')'Direct'
      WRITE(IU,'(3F10.6)')(chg%ionpos_dir(I,1:3),I=1,chg%nions)
      WRITE(IU,*)
!      FORM='(1(1X,E17.11))'
      FORM='(1(1X,G12.5))'
      NWRITE=5
      WRITE(IU,'(3I5)')chg%np(1),chg%np(2),chg%np(3)
      NWRITTEN=0
      DO k=1,chg%nrho
           NWRITTEN=NWRITTEN+1
           IF (MOD(NWRITTEN,NWRITE).eq.0) THEN
              WRITE(IU,FORM) dens(k)
           ELSE
              WRITE(IU,FORM,ADVANCE='NO') dens(k)
           ENDIF
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE).ne.0 ) WRITE(IU,*)''
    CLOSE(IU)
    WRITE(*,'(1A10,1A20)') 'CLOSE ... ',writeoutfile

    RETURN
  END SUBROUTINE write_chgcar

!-----------------------------------------------------------------------------!
!copy one charge object to another
!-----------------------------------------------------------------------------!

  SUBROUTINE assign_charge(chg1,chg2)

    TYPE(charge), INTENT(INOUT) :: chg1
    TYPE(charge), INTENT(IN) :: chg2

    chg1%lat2car=chg2%lat2car
    chg1%np=chg2%np
    chg1%nrho=chg2%nrho
    chg1%nions=chg2%nions

    ALLOCATE(chg1%rho(chg1%nrho))
    chg1%rho=chg2%rho

    END SUBROUTINE

!-----------------------------------------------------------------------------!
!wrap 3D point (p(1),p(2),p(3)) to [1,pmax]
!-----------------------------------------------------------------------------!

  SUBROUTINE pbc3D(p,pmax)

    INTEGER,INTENT(INOUT),DIMENSION(3) :: p
    INTEGER,INTENT(IN),DIMENSION(3) :: pmax
    INTEGER :: I

    DO I=1,3
      DO
        IF(p(I) > 0) EXIT
        p(I)=p(I)+pmax(I)
      END DO
      DO
        IF(p(I) <= pmax(I)) EXIT
        p(I)=p(I)-pmax(I)
      END DO
    END DO

  END SUBROUTINE pbc3D

!-----------------------------------------------------------------------------!
! wrap 1D point p to [1,pmax]
!-----------------------------------------------------------------------------!

  SUBROUTINE pbc1D(p,pmax)

    INTEGER,INTENT(INOUT) :: p
    INTEGER,INTENT(IN) :: pmax

      DO
        IF(p > 0) EXIT
        p=p+pmax
      END DO
      DO
        IF(p <= pmax) EXIT
        p=p-pmax
      END DO

  END SUBROUTINE pbc1D

!-----------------------------------------------------------------------------!
!Find the smallest distance vector in Cartesian coordinates
!-----------------------------------------------------------------------------!

  SUBROUTINE pbc_dcar(cell,dcar)
    
    REAL(q),INTENT(IN), DIMENSION(3,3) :: cell 
    REAL(q),INTENT(INOUT),DIMENSION(3) :: dcar
    REAL(q),DIMENSION(3) :: dcart,dcarmin,lat1,lat2,lat3
    REAL(q) :: dsq,dsqmin
    INTEGER :: d1,d2,d3

    dcarmin=dcar
    dsqmin=dot(dcar,dcar)
      DO d1=-1,1
         lat1=cell(1,:)*REAL(d1,q)
      DO d2=-1,1
         lat2=cell(2,:)*REAL(d2,q)
      DO d3=-1,1
         lat3=cell(3,:)*REAL(d3,q)

         dcart=dcar+lat1+lat2+lat3
         dsq=dot(dcart,dcart)
         IF(dsq<dsqmin) THEN
            dcarmin=dcart
            dsqmin=dsq
         END IF

      END DO
      END DO
      END DO
      dcar=dcarmin

  END SUBROUTINE pbc_dcar

!-----------------------------------------------------------------------------!
!Evaluate max. atomic distance
!-----------------------------------------------------------------------------!
  SUBROUTINE atomdist(chg)

    TYPE(charge) :: chg
    INTEGER,DIMENSION(3) :: p
    REAL(q),DIMENSION(3) :: dlat,dcar
    REAL(q) :: dsq
    INTEGER :: I,J

    chg%dminsq = 1.0e8_q
    IF (chg%nions == 1) THEN
        WRITE (*,'(2X,A)')'WARNING: Can not calculate atomic distance,  &
&                                   system contains only one atom.'
    chg%dminsq = 1.0e8_q
    ELSE 
      DO I=1,chg%nions-1
       DO J = I+1, chg%nions
        dlat = chg%ionpos_lat(I,:)-chg%ionpos_lat(J,:)
        CALL mult_matrix_vector(chg%lat2car,dlat,dcar)
        CALL pbc_dcar(chg%SCALE*chg%A,dcar)
        dsq = dot(dcar,dcar)
        IF (dsq < chg%dminsq) THEN
           chg%dminsq = dsq
        END IF
       END DO
      END DO
    END IF

  END SUBROUTINE atomdist


END MODULE chgcar
