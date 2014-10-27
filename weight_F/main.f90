!-----------------------------------------------------------------------------------!
! Wigner-Seitz flux algorithm 
! Version 1/1/11
!
! Authors:
!  Min Yu, Dallas Trinkle
!
! Based on the following publication:  
!  Accurate and efficient algorithm for Bader charge integration
!  arXiv:1010.4916  J. Chem. Phys. 134, 064111 (2011)
!-----------------------------------------------------------------------------------!

  PROGRAM Main 
     USE prec
     USE options
     USE chgcar
     USE weights

     IMPLICIT NONE

     TYPE(opt) :: opts
     TYPE(charge) :: chgval             !charge/energy density
     TYPE(charge) :: chgref 
       
! Write the version number
     WRITE(*,'(2X,A)') 'Wigner-Seitz flux algorithm  (Version 1/1/11)'

! Get arguments
     CALL get_opts(opts)
     WRITE(*,*) ''    
     WRITE(*,'(2X,A)') 'READ IN charge density'
     CALL read_chgcar(opts%chgcarfile,chgval)

   IF (opts%ref_flag) THEN
     WRITE(*,*) ''
     WRITE(*,'(2X,A)') 'READ IN reference charge density'
     CALL read_chgcar(opts%refchgcarfile,chgref)
   ELSE
     chgref = chgval
   END IF

! Construct the Wigner-Seitz Cell
     CALL wsflux(chgref,opts) 

! Calculate atomic weight
     CALL weight_calc(chgval,chgref,opts)

  END PROGRAM Main
