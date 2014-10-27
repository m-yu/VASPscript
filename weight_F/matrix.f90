MODULE matrix
  USE prec
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mult_matrix_vector
  PUBLIC :: matrix_transpose,matrix_inverse
  PUBLIC :: vector_tripleprod, dot, crossprod

  CONTAINS

!-----------------------------------------------------------------------------!
! mult_matrix_vector:  Vector C = Matrix A * Vector B 
!-----------------------------------------------------------------------------!

  SUBROUTINE mult_matrix_vector(A,B,C)

    REAL(q),INTENT(IN),DIMENSION(3,3) :: A
    REAL(q),INTENT(IN),DIMENSION(3) :: B
    REAL(q),INTENT(OUT),DIMENSION(3) :: C
    INTEGER :: I

    C=0._q
    DO I=1,3
      C=C+A(:,I)*B(I)
    END DO

  END SUBROUTINE mult_matrix_vector

!-----------------------------------------------------------------------------!
! matrix_transpose:  Matrix B = Matrix A ^{T}
!-----------------------------------------------------------------------------!
  SUBROUTINE matrix_transpose(A,B)

    REAL(q),INTENT(IN),DIMENSION(3,3) :: A
    REAL(q),INTENT(OUT),DIMENSION(3,3) :: B
    INTEGER :: I,J

    DO I=1,3
      DO J=1,3
        B(J,I)=A(I,J)
      END DO
    END DO

  END SUBROUTINE matrix_transpose

!-----------------------------------------------------------------------------!
! matrix_inverse: Matrix B = Matrix A ^{-1}
!-----------------------------------------------------------------------------!
  SUBROUTINE matrix_inverse(A,B,det)

    REAL(q),INTENT(IN),DIMENSION(3,3) :: A
    REAL(q),INTENT(OUT),DIMENSION(3,3) :: B
    REAL(q) :: det

    det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
    &  + A(1,2)*(A(2,3)*A(3,1) - A(2,1)*A(3,3)) &
    &  + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

    B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
    B(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
    B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)    
    B(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    B(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)    
    B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
    B(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
    B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)   

    B(:,:)=B(:,:)/det

  END SUBROUTINE matrix_inverse

!-----------------------------------------------------------------------------!
! vector_tripleprod: Return the triple product of three vectors.
!-----------------------------------------------------------------------------!
  FUNCTION vector_tripleprod(A,B,C)

    REAL(q) :: A(3),B(3),C(3)   
    REAL(q) :: vector_tripleprod
    
    vector_tripleprod = C(1)*(A(2)*B(3) - A(3)*B(2)) &
    &  + C(2)*(A(3)*B(1) - A(1)*B(3)) &
    &  + C(3)*(A(1)*B(2) - A(2)*B(1))

    RETURN
  END FUNCTION vector_tripleprod


!-----------------------------------------------------------------------------!
!dot: Return dot product of two vectors
!-----------------------------------------------------------------------------!
  FUNCTION dot (A, B)

    REAL(q),DIMENSION(3) :: A,B 
    REAL(q) :: dot 

    dot = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)

    RETURN
  END FUNCTION dot 

!-----------------------------------------------------------------------------------!
!crossprod: Vectors C = Vector A X Vector B 
!-----------------------------------------------------------------------------------!
  SUBROUTINE crossprod(A,B,C)

    REAL(q),INTENT(IN),DIMENSION(3) :: A,B
    REAL(q),INTENT(OUT),DIMENSION(3) :: C

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

  END SUBROUTINE  crossprod

END MODULE matrix

