
SUBROUTINE wsflux (chg,opts) 
  USE prec
  USE options
  USE matrix
  USE chgcar
  IMPLICIT NONE
  TYPE(opt) :: opts
  TYPE(charge) :: chg

  INTEGER,PARAMETER :: Nrange=3	! (generally) overkill; need a way to compute this        
  INTEGER :: Nneigh  
  REAL(q),ALLOCATABLE,DIMENSION(:,:) :: R
  REAL(q) :: Rd2(3)
  INTEGER :: MAXVERT, Nnonzero
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: vert
  REAL(q),ALLOCATABLE,DIMENSION(:) :: alpha
  REAL(q),ALLOCATABLE,DIMENSION(:,:) :: rvert 
  REAL(q) :: Rdot(3,3), detR, Rinv(3,3), R2(3)
  INTEGER :: IX,IY,IZ
  INTEGER :: I,J,d,n,nvert,nv,nvp,nA,nB
  REAL(q) :: rx(3), ry(3), rdRn
  LOGICAL :: incell,zeroarea
 
  Nneigh = (2*Nrange+1)*(2*Nrange+1)*(2*Nrange+1)-1
  ALLOCATE(R(Nneigh,4),vert(Nneigh,3))

!-----------------------------------------------------------------------------!  
!ANALYSIS
! 0. Generate (pruned) list of neighboring vectors that bound 
!the Wigner-Seitz cell.
!-----------------------------------------------------------------------------!   
  Nneigh=0

  DO IX=-Nrange, Nrange
  DO IY=-Nrange, Nrange
  DO IZ=-Nrange, Nrange
	IF (IX .ne. 0 .OR. IY .ne. 0 .OR. IZ .ne. 0) THEN
	  Nneigh = Nneigh+1
          vert(Nneigh,:)=(/IX,IY,IZ/)
	  CALL mult_matrix_vector(chg%lat2car, REAL(vert(Nneigh,:),q), R(Nneigh,1:3))
	  R(Nneigh,4) = 0.5*dot(R(Nneigh,1:3), R(Nneigh,1:3))
	END IF
  ENDDO
  ENDDO
  ENDDO

!-----------------------------------------------------------------------------!
! prune that list: basically, if R/2 isn't inside the WS cell, 
! then R isn't involved in bounding the WS cell.
!-----------------------------------------------------------------------------!
  n=Nneigh
  Nneigh =0
  DO I=1,n
    Rd2 = 0.5_q*R(I,1:3)
    incell =.TRUE.
    DO J=1, n
    IF (dot(Rd2, R(J,1:3)) > (R(J,4)+opts%WS_tol)) incell = .False.
    ENDDO
    IF (incell) THEN
       Nneigh = Nneigh+1
       vert(Nneigh,:) = vert(I,:)
       R(Nneigh,:) = R(I,:)
    ENDIF
  ENDDO
    DO I = Nneigh+1,n
       vert(I,:) =0 
       R(I,:)=0  
    ENDDO

!-----------------------------------------------------------------------------!
! 1. Run over each facet to find all of the vertex points
!-----------------------------------------------------------------------------!
  MAXVERT = (Nneigh-2)*(Nneigh-4)
  ALLOCATE(rvert(MAXVERT,4),alpha(Nneigh))
  Nnonzero = Nneigh
  DO n=1, Nneigh
     nvert = 1
     Rdot(1,1:3) = R(n,1:3)
     R2(1) = R(n,4)
  DO nA=1, Nneigh
     Rdot(2,1:3) = R(nA,1:3)
     R2(2) = R(nA,4)
  DO nB=(nA+1), Nneigh
     Rdot(3,1:3) = R(nB,1:3)
     R2(3) = R(nB,4)
     CALL matrix_inverse(Rdot, Rinv,detR)
     IF(ABS(detR)>opts%WS_tol) THEN
	  CALL mult_matrix_vector(Rinv, R2, rvert(nvert,1:3))
          incell =.TRUE.
          DO J=1, Nneigh
          IF (dot(rvert(nvert,1:3), R(J,1:3)) > (R(J,4)+opts%WS_tol)) incell = .False.
          ENDDO
	  IF (incell) nvert=nvert+1
     END IF
    ENDDO	    
    ENDDO
    nvert=nvert-1
!-----------------------------------------------------------------------------!
! check to make sure none of the vertices correspond to R/2:
!-----------------------------------------------------------------------------!
    zeroarea = .False.  
    DO nv=1, nvert
      IF (ABS(dot(rvert(nv,:), rvert(nv,:))-0.5*R(n,4))<opts%WS_tol) zeroarea=.True.
    ENDDO
    IF (zeroarea .OR. (nvert.eq.0)) THEN
      alpha(n)=0
      Nnonzero=Nnonzero-1
      GOTO 100
    ENDIF
!-----------------------------------------------------------------------------!
! Now we have a list of all the vertices for the polygon
! defining the facet along the direction R[n].
! Last step is to sort the list in terms of a winding angle around R[n].  
! To do that, we define rx and ry which are perpendicular
! to R[n], normalized, and right-handed: ry = R x rx, so that
! rx x ry points along R[n].
!-----------------------------------------------------------------------------!
        rx(:) = rvert(1,:)
        rdRn = dot(rx, R(n,1:3))/dot(R(n,1:3),R(n,1:3))
        rx(:) = rx(:)-rdRn*R(n,1:3)
        rdRn = sqrt(dot(rx, rx))
        rx(:) = rx(:)/rdRn
        CALL crossprod(R(n,1:3), rx, ry)
        rdRn = sqrt(dot(ry, ry))
        ry(:) = ry(:)/rdRn
! now compute winding angle phi_n
      DO nv=1, nvert
         rvert(nv,4) = atan2(dot(rvert(nv,1:3), ry), dot(rvert(nv,1:3), rx))
      ENDDO
! sort that list
      CALL sort(rvert(1:nvert,:), nvert,opts%WS_tol)
      alpha(n) = 0
      DO nv=1, nvert
      nvp=nv+1
      CALL pbc1D(nvp,nvert)
      alpha(n) = alpha(n) + vector_tripleprod(rvert(nv,1:3), rvert(nvp,1:3), R(n,1:3))
      ENDDO
      alpha(n) = alpha(n) *0.25/R(n,4)
      IF (ABS(alpha(n))<opts%WS_tol) THEN
      alpha(n)=0
      Nnonzero=Nnonzero-1
      ENDIF 
100 continue
  ENDDO

  DEALLOCATE(rvert,R)
!-----------------------------------------------------------------------------!
!                                         OUTPUT
!-----------------------------------------------------------------------------! 
  chg%Nneigh = Nnonzero
  ALLOCATE(chg%nvert(chg%Nneigh,3),chg%alpha(chg%Nneigh))
  WRITE(*,*) ''
  WRITE(*,'(2X,A,I3)')'Num of Wigner-Seitz neighbors : ',chg%Nneigh
  WRITE(*,*) '          X           Y           Z          alpha'
  chg%Nneigh =0
  DO n=1, Nneigh
    IF (alpha(n).ne.0) THEN
    chg%Nneigh = chg%Nneigh+1 
    chg%nvert(chg%Nneigh,:)=vert(n,:)   
    chg%alpha(chg%Nneigh)=alpha(n)
     WRITE(*,*)chg%nvert(chg%Nneigh,1:3),chg%alpha(chg%Nneigh)
    ENDIF
   ENDDO

END SUBROUTINE


!-----------------------------------------------------------------------------!
!arr(n,1:4) record the positions of arr(n,1:3)
!sort arr(n,1:4) into ascending order of arr(n,4) 
!-----------------------------------------------------------------------------!
  SUBROUTINE sort(arr,n,tol)
  USE prec
  INTEGER,INTENT(IN) :: n
  REAL(q),INTENT(INOUT),DIMENSION(n,4)  :: arr
  REAL(q),INTENT(IN) :: tol  
  INTEGER :: I,J,K
  REAL(q), DIMENSION(4) :: a

  DO I=1,n
   a = arr(I,:)
   K = I
   DO J=I+1,n
     IF (arr(J,4) < a(4)) THEN
       a = arr(J,:)
       K = J 
     END IF 
   ENDDO
       arr(K,:) = arr(I,:)
       arr(I,:) = a
  ENDDO

  END SUBROUTINE sort

