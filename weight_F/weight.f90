!-----------------------------------------------------------------------------!
! Atomic weight evaluation
!-----------------------------------------------------------------------------!

MODULE weights

  USE prec
  USE options
  USE matrix
  USE chgcar

  IMPLICIT NONE
!-----------------------------------------------------------------------------!
! Variable description:
! w(num_edge) : atomic weight on every mesh point
! ionden(nnions) : integrated density for each atom
! ionvol(nnions) : atomic volume
! num_edge : number of edge pts
! ionden : density summation for every ion 
! assign_2_ion : 
!   positive: the n.n. atomic index that a charge density maximal point is assigned to.
!   zero: edge points 
!-----------------------------------------------------------------------------!
  TYPE weight
    INTEGER :: num_edge
    REAL(q),ALLOCATABLE,DIMENSION(:) :: w,ionden,ionvol
    INTEGER,ALLOCATABLE,DIMENSION(:) :: assign_2_ion
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: neigh 
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: prob
  END TYPE

  PUBLIC :: weight
  PUBLIC :: weight_calc

  CONTAINS

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
  SUBROUTINE weight_calc(chgval,chgref,opts)

    TYPE(charge) :: chgref   !reference (charge/potential) density
    TYPE(charge) :: chgval   !energy density
    TYPE(opt) :: opts
    TYPE(weight) :: wt
!-----------------------------------------------------------------------------!
! Variable description:
! atmvol : atomic assignment according to max. weight
! wt_max : max. weight on each grid pt
!-----------------------------------------------------------------------------!
    INTEGER,DIMENSION(3) :: p,pt
    INTEGER :: d,I,J,k,kt,ion
    INTEGER :: num_assign,step
    INTEGER,ALLOCATABLE,DIMENSION(:) ::  neigh
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: R
    REAL(q),ALLOCATABLE,DIMENSION(:) ::  rhoR, grad, prob, arr_rho
    INTEGER,ALLOCATABLE,DIMENSION(:) :: arr_N,arr_N_back    
    LOGICAL :: Ledge_pt
    INTEGER :: t0,t1,t2,t3,count_rate,count_max
    INTEGER,ALLOCATABLE,DIMENSION(:) :: atmvol
    REAL(q),ALLOCATABLE,DIMENSION(:) :: wt_max
    REAL(q),ALLOCATABLE,DIMENSION(:) :: wden
    CHARACTER(LEN=40) :: writeoutfile    
    REAL(q) :: prob_norm,w_t

    CALL SYSTEM_CLOCK(t0,count_rate,count_max)

    ALLOCATE(wt%neigh(chgref%nrho,0:chgref%Nneigh),wt%prob(chgref%nrho,chgref%Nneigh))
    ALLOCATE(wt%assign_2_ion(chgref%nrho))
    ALLOCATE(arr_rho(chgref%nrho),arr_N(chgref%nrho),arr_N_back(chgref%nrho))
    ALLOCATE(neigh(0:chgref%Nneigh),R(chgref%Nneigh,1:3),rhoR(0:chgref%Nneigh))
    ALLOCATE(grad(chgref%Nneigh),prob(chgref%Nneigh))
    wt%assign_2_ion = 0
    arr_rho(:) = -1._q*chgref%rho(:)
    arr_N = 0
    arr_N_back = 0

! SORTRX return arr_N: indices of elements of rho arranged in descending order. 
! rho(arr_N(1)) will be the largest number in rho;
! rho(arr_N(N)) will be the smallest number in rho;
    WRITE(*,*) ''
    WRITE(*,'(2X,A)')'QUICK SORTING: reference charge density'
    CALL SYSTEM_CLOCK(t1,count_rate,count_max)
    CALL SORTRX(chgref%nrho,arr_rho,arr_N)
    CALL SYSTEM_CLOCK(t2,count_rate,count_max)
    WRITE(*,'(2X,A,1F7.2,1A8)') 'RUN TIME (SORTING): ',(t2-t1)/REAL(count_rate,q),'SECONDS'
! map arr_N back
    DO k = 1, chgref%nrho
      arr_N_back(arr_N(k)) = k
    END DO

    wt%neigh = 0
    wt%prob = 0
    wt%num_edge = 0
    CALL SYSTEM_CLOCK(t1,count_rate,count_max)
    WRITE(*,*) ''
    WRITE(*,'(2X,A)') 'CALCULATING JUMPING PROBABILITY'
    arrN: DO k = 1, chgref%nrho
      Ledge_pt = .FALSE.
      I = arr_N(k)
      CALL index1D_2_3D(I,chgref%np,p) 
! store indices of n.n. grid points; calculate charge densties.
! R(1:chg%Nneigh,1:3) : 3D coordinate indices of n.n. grid points 
! neigh(1:chg%Nneigh) : 1D indices of n.n. grid points
! rhoR(0:chg%Nneigh) : charge density of one grid point with its n.n.    
      rhoR(0) = chgref%rho(I)
      prob_norm = 0._q
      neigh = 0 
      grad = 0._q
      DO d = 1, chgref%Nneigh
        R(d,:) = p + chgref%nvert(d,:)
        CALL pbc3D(R(d,:),chgref%np)              
        J = chgref%np(1)*chgref%np(2)*(R(d,3)-1)+chgref%np(1)*(R(d,2)-1)+R(d,1)        
        rhoR(d) = chgref%rho(J) 
! calcualte jumping prob to the neighboring grid points with higher magnitude
        kt = arr_N_back(J)
        IF (kt < k) THEN  ! equal to (rho(J) > rho(I))
          neigh(0) = neigh(0) + 1 
          neigh(neigh(0)) = kt
          grad(neigh(0)) = (rhoR(d)-rhoR(0))*chgref%alpha(d)
          prob_norm = prob_norm + grad(neigh(0))   
        END IF
      END DO

! a grid point with density max.      
      IF (neigh(0) == 0) THEN
        wt%assign_2_ion(k) = assign_chg2atom(chgval,p,wt)       
      ELSE 
! check whether a grid point is an interior or a boundary point        
        ion = wt%assign_2_ion(neigh(1))
        Nneigh: DO d = 1, neigh(0)                 
! one grid point is defined as boundary point if        
! either neighbor with higher magnitude is assigned to different ion 
! or is a boundary point
          IF ( ion == 0 .or. wt%assign_2_ion(neigh(d)) /= ion ) THEN
            Ledge_pt = .TRUE.
            EXIT Nneigh  
          END IF
        END DO Nneigh
      END IF

      IF (Ledge_pt) THEN        ! a boundary point
        wt%num_edge = wt%num_edge + 1
        wt%neigh(k,:) = neigh(:) 
        wt%prob(k,:) = grad(:)/prob_norm 
      ELSE IF (neigh(0) /= 0) THEN !an interior point but not max. 
        wt%assign_2_ion(k) = ion
      END IF   

    END DO arrN
    DEALLOCATE(arr_rho,grad,neigh)

    WRITE(*,'(2X,A,I20)') 'num_edge=',wt%num_edge
    CALL SYSTEM_CLOCK(t2,count_rate,count_max)
    WRITE(*,'(2X,A,1F7.2,1A8)') 'RUN TIME (CALC. INI. PROB): ',(t2-t1)/REAL(count_rate,q),'SECONDS'

    ALLOCATE(wt%w(chgref%nrho)) 
    ALLOCATE(wt%ionden(chgref%nions),wt%ionvol(chgref%nions))
    wt%ionden = 0._q
    wt%ionvol = 0._q

    IF (opts%write_atmvol) THEN
    ALLOCATE(atmvol(chgref%nrho),wt_max(wt%num_edge))
    atmvol = 0
    wt_max = 0._q
    END IF
    IF (opts%write_weight) THEN
    ALLOCATE(wden(chgref%nrho))
    END IF      

    WRITE(*,*) ''
    WRITE(*,'(2X,A,I8,2X,A)') 'CALCULATING ATOMIC WEIGHT ON',chgref%nions,'IONS'
    WRITE(*,'(2X,A,$)') 'ION DONE:'
    nions: DO ion = 1, chgref%nions
      WRITE(*,'(I3,A,$)') ion,' '
      wt%w = 0._q
      IF (opts%write_weight .and. &
&         (opts%weight_atm == ion .or.  opts%weight_atm == 0)) wden = 0._q
      arrN2: DO k = 1, chgref%nrho
        IF (wt%assign_2_ion(k) == ion) THEN 
          wt%w(k) = 1._q
        ELSE IF (wt%assign_2_ion(k) == 0) THEN 
          DO d = 1, wt%neigh(k,0) 
            wt%w(k) = wt%w(k) + wt%w(wt%neigh(k,d))*wt%prob(k,d)
          END DO    
        END IF
! calculate integral
        wt%ionden(ion) = wt%ionden(ion) + chgval%rho(arr_N(k)) * wt%w(k)
        wt%ionvol(ion) = wt%ionvol(ion) + wt%w(k)
        IF (opts%write_atmvol .and. wt%w(k) > wt_max(k)) THEN
          atmvol(arr_N(k)) = ion
          wt_max(k) = wt%w(k)
        END IF
        IF (opts%write_weight .and. &
&         (opts%weight_atm == ion .or.  opts%weight_atm ==0)) &
&         wden(arr_N(k)) = wt%w(k)         
      END DO arrN2
       wt%ionden(ion) = wt%ionden(ion) / REAL(chgval%nrho,q)
       wt%ionvol(ion) = wt%ionvol(ion) / REAL(chgval%nrho,q)

! write out atomic weight function
      IF (opts%write_weight .and. &
&         (opts%weight_atm == ion .or.  opts%weight_atm ==0)) THEN
        WRITE(writeoutfile,'(A4,I3.3,A4)') "wion",ion,".dat"
        CALL write_chgcar(writeoutfile,chgval,wden) 
      END IF 

    END DO nions
    
 
     
!-----------------------------------------------------------------------------!
!                            OUTPUT weight
!-----------------------------------------------------------------------------! 
    WRITE(*,*) ''
    WRITE(*,'(2X,A)')'WRITEOUT:  ion      density integration       ion vol      '
    DO ion = 1, chgref%nions 
      WRITE(*,'(12X,I3,2(10X,1F10.5))')ion,wt%ionden(ion),wt%ionvol(ion)
    END DO  

! write out atomic assignment 
    IF (opts%write_atmvol) THEN
      writeoutfile = 'atmvol'
      CALL write_chgcar(writeoutfile,chgval,REAL(atmvol,q)) 
    END IF 

    CALL SYSTEM_CLOCK(t3,count_rate,count_max)
    WRITE(*,'(2X,A,1F7.2,1A8)') 'RUN TIME (AVERAGE ATOMIC WEIGHT): ',(t3-t2)/REAL(count_rate,q),'SECONDS'
    WRITE(*,*) ''
    WRITE(*,'(2X,A,1F7.2,1A8)') 'RUN TIME (TOTAL WEIGHT ANALYSIS): ',(t3-t0)/REAL(count_rate,q),'SECONDS'
  END SUBROUTINE weight_calc

!-----------------------------------------------------------------------------!
!Convert a grid point index from 1D to 3D.
!-----------------------------------------------------------------------------!
  SUBROUTINE index1D_2_3D(ind,np,p) 
    INTEGER,INTENT(IN) ::  ind
    INTEGER,DIMENSION(3),INTENT(IN) ::  np
    INTEGER,DIMENSION(3),INTENT(OUT) ::  p
    INTEGER :: n,n1,n2,n3
        
    n = ind
    n1 = mod(n,np(1)) 
    IF (n1 .eq. 0)  n1 = np(1)
    n = n - n1
    n2 = mod(n/np(1),np(2)) + 1
    IF (n2 .eq. 0)  n2 = np(2)
    n = n - np(1)*(n2-1) 
    n3 = n/(np(1)*np(2)) + 1
    p = (/n1,n2,n3/)

  END SUBROUTINE index1D_2_3D

!-----------------------------------------------------------------------------!
! Assign a charge density maximal point to n.n. atom
!-----------------------------------------------------------------------------!
  FUNCTION assign_chg2atom(chg,p,wt)

    TYPE(charge) :: chg
    TYPE(weight) :: wt
    INTEGER,DIMENSION(3) :: p
    REAL(q),DIMENSION(3) :: dlat,dcar
    REAL(q) :: dsq,dminsq
    INTEGER :: I,J,dindex,assign_chg2atom

      dlat = REAL(p(:),q)-chg%ionpos_lat(1,:)
      CALL mult_matrix_vector(chg%lat2car,dlat,dcar)
      CALL pbc_dcar(chg%SCALE*chg%A,dcar)
      dminsq = dot(dcar,dcar)
      dindex = 1
      DO J=2,chg%nions
        dlat = REAL(p(:),q)-chg%ionpos_lat(J,:)
        CALL mult_matrix_vector(chg%lat2car,dlat,dcar)
        CALL pbc_dcar(chg%SCALE*chg%A,dcar)
        dsq = dot(dcar,dcar)
        IF (dsq < dminsq) THEN
          dminsq = dsq
          dindex = J
        END IF
      END DO
      assign_chg2atom = dindex

  RETURN
  END FUNCTION assign_chg2atom

  END MODULE weights 
