  MODULE options
    USE prec
    IMPLICIT NONE

!-----------------------------------------------------------------------------!
! Variable description:
! write_atmvol : .True. write out atomic assignment  
! write_weight : .True. write out atomic weight
!-----------------------------------------------------------------------------!
    TYPE :: opt
      CHARACTER(128) :: chgcarfile,refchgcarfile  
      LOGICAL :: ref_flag,write_atmvol,write_weight
      REAL(q) :: weight_tol, WS_tol
      INTEGER :: weight_atm
    END TYPE opt

    PRIVATE 
    PUBLIC :: get_opts,opt

    CONTAINS

!-------------------------------------------------------------------------------------!
! get_options: Read any input flags and the charge density file name
!-------------------------------------------------------------------------------------!

    SUBROUTINE get_opts(opts)

      TYPE(opt) :: opts
      LOGICAL :: existflag
      LOGICAL :: readchgflag
      INTEGER :: n,i,ip,m,it,ini,sel,istart,iend
      CHARACTER(128) :: p
      CHARACTER(128) :: inc
      INTEGER :: COMMAND_ARGUMENT_COUNT

! Default values
      opts%ref_flag = .FALSE.
      opts%write_atmvol = .FALSE.
      opts%write_weight = .FALSE.
      opts%WS_tol = 1.0e-8_q
      opts%weight_atm = 0
      n=COMMAND_ARGUMENT_COUNT()

! Loop over all arguments
      m=0
      readchgflag = .FALSE.
      readopts: DO WHILE(m<n)
100       m=m+1
        CALL GET_COMMAND_ARGUMENT(m,p)
        p=ADJUSTL(p)
        ip=LEN_TRIM(p)
        i=INDEX(p,'-')
!Charge density file name      
        IF (i /= 1) THEN
          IF (readchgflag) THEN
            WRITE(*,'(3A)') ' Option "',p(1:ip),'" is not valid'
            STOP
          END IF
          opts%chgcarfile=p
          INQUIRE(FILE=opts%chgcarfile,EXIST=existflag)
          IF (.NOT. existflag) THEN
            WRITE(*,'(2X,A)') opts%chgcarfile(1:ip),' does not exist'
            STOP
          END IF
          readchgflag=.TRUE.
!Help
        ELSEIF (p(1:ip) == '-h') THEN
          CALL write_expl()
          STOP
!Write out atomic assignment    
        ELSEIF (p(1:ip) == '-pv') THEN
          m=m+1
          opts%write_atmvol = .TRUE.
!Write out atomic weight 
        ELSEIF (p(1:ip) == '-pw') THEN  
          m=m+1
          opts%write_weight = .TRUE.
          CALL GET_COMMAND_ARGUMENT(m,inc)
          read(inc,*) opts%weight_atm
write(*,*) 'opts%weight_atm',opts%weight_atm
!Doing analysis with reference charge file
        ELSEIF (p(1:ip) == '-ref') THEN
          m=m+1
          CALL GET_COMMAND_ARGUMENT(m,inc)
          inc=ADJUSTL(inc)
          it=LEN_TRIM(inc)
          IF (inc(1:it) == 'NONE' .OR. inc(1:it) == 'none') THEN
            opts%ref_flag = .FALSE.
          ELSE
            opts%ref_flag = .TRUE.
            opts%refchgcarfile = inc(1:it)
          END IF
        ENDIF

      END DO readopts

    IF (.NOT. readchgflag) THEN
      WRITE(*,*) ' ERROR: expect to read a charge file name in the arguments'
      STOP
    ENDIF

    RETURN
    END SUBROUTINE get_opts


!-----------------------------------------------------------------------------------!
! write_expl: write help
!-----------------------------------------------------------------------------------!

    SUBROUTINE write_expl()

      WRITE(*,*) ''
      WRITE(*,*) 'Description of flags'
      WRITE(*,*) ''  
      WRITE(*,*) '   -h'
      WRITE(*,*) '        Help'
      WRITE(*,*) ''
      WRITE(*,*) '   -pv'
      WRITE(*,*) '        Write out atomic assignment: atmvol'
      WRITE(*,*) ''
      WRITE(*,*) '   -pw'
      WRITE(*,*) '        Write out atomic weight for selected atom:'
      WRITE(*,*) '        0 for all atoms; 1...n '
      WRITE(*,*) ''
      WRITE(*,*) '   -ref'
      WRITE(*,*) '        Analyze over reference charge file '
      WRITE(*,*) ''

    END SUBROUTINE write_expl

!-----------------------------------------------------------------------------------!

  END module options

