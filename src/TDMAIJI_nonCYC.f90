!***********************************************************************
!> @author 
!> Wei Wang, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!>  SOLVER FOR REAL-DEFINITE, CYCLIC-TRIDIAGONAL, CONSTANT COEFFICIENTS
!>  LINEAR SYSTEMS OF EQUATIONS.
!>  THE SOLVER IS WRITTEN IN FORTRAN 90 LANGUAGE AND ADOPTS ONLY
!>  DOUBLE-PRECISION (8 BYTES) REAL VARIABLES.
!>********************************************************************************
!>  @par Methodology                                                             *
!>  THE METHOD SOLVES FIRSTLY FOR THE X(1) UNKNOWN, WHICH MAKES THE              *
!>  COEFFICIENT MATRIX CYCLIC: THIS IS THE INITIALIZATION PHASE, PERFORMED       *
!>  IN THE PRESENT SUBROUTINE. NOW THE SYSTEM IS A (N-1) BY (N-1) TRIDIAGONAL    *
!>  (NOT CYCLIC) SYSTEM, AND IS SOLVED WITH THE STANDARD (AND FAST) LU           *
!>  FACTORIZATION TECHNIQUE (ROUTINES LU_PRE AND LU_SOLV).                       *
!>                                                                               *
!>  COEFFICIENTS ON MAIN AND SECONDARY DIAGONALS MUST BE CONSTANT, I.E. WE ARE   *
!>  DEALING WITH SYSTEMS OF THE FOLLOWING KIND:                                  *
!>                                                                               *
!>                      | A  B  0  0  B | |X1|   |D1|                            *
!>                      | B  A  B  0  0 | |X2|   |D2|                            *
!>                      | 0  B  A  B  0 | |X3| = |D3|                            *
!>                      | 0  0  B  A  B | |X4|   |D4|                            *
!>                      | B  0  0  B  A | |X5|   |D5|                            *
!>                                                                               *
!>  OF COURSE, THE IMPROUVEMENT IN RESPECT OF USING THE CLASSIC AND EFFICIENT    *
!>  CTDMA SOLVER OF THOMAS ARISES WHEN ONE NEEDS TO SOLVE MANY CYCLIC TRIDIAGONAL*
!>  SYSTEMS WITH THE SAME COEFFICIENT MATRIX.                                    *
!>                                                                               *
!>  NOTE! IN THE INITIALIZATION PHASE WE PERFORM ALSO THE LU FACTORIZATION       *
!>  OF THE N-1 BY N-1 TRIDIAGONAL COEFFICIENT MATRIX.                            *
!>                                                                               *
!>  THE NON-UNITY DIAGONALS OF MATRICES L & U ARE STORED IN THE ARRAYS            * 
!>  DIAG_PRIN(1:N-1) & DIAG_SUPE(1:N-2) (IN THIS ORDER)                           *
!>                                                                                *
!>  THE RIGHT-HAND SIDE IS STORED IN THE ARRAY RHS (N ELEMENTS). IN               *
!>  OUTPUT FROM SUB. ALG4_SOL RHS CONTAINS THE SOLUTION.                          *
!>                                                                                *
!>  Z IS A REAL*8 ARRAY OF N ELEMENTS, AND DOESN'T NEED TO BE INITIALIZED.        *
!>                                                                                *
!> @par USAGE:                                                                    *
!>      CALL TDMAIJI_CYC(A,B,C,R,IS,ISZ,JS,JSZ)                                   *
!>      TRIDIAGNAL MARCH in I direction in the plane of I*J
!>
!> @par INPUT                                                                     *
!>      IS:  INTEGER, The starting index of I direction in the I*J plane          *
!>      ISZ: INTEGER, The total size of points in I direction                     *
!>      JS:  INTEGER, The starting index of J direction in the I*J plane          *
!>      JSZ: INTEGER, The total size of points in J direction                     *
!>      A(IS:1:IS+ISZ-1, JS:1:JS+JSZ-1) : Matrix ISZ*JSZ, Coef of u_{i-1,j}       *
!>      B(IS:1:IS+ISZ-1, JS:1:JS+JSZ-1) : Matrix ISZ*JSZ, Coef of u_{i,j}         *
!>      C(IS:1:IS+ISZ-1, JS:1:JS+JSZ-1) : Matrix ISZ*JSZ, Coef of u_{i+1,j}       *
!>      R(IS:1:IS+ISZ-1, JS:1:JS+JSZ-1) : Matrix ISZ*JSZ, RHS of the original     *
!>                                        tridiagnal sys                          *
!>
!>par   OUTPUT
!>      R(IS:1:IS+ISZ-1, JS:1:JS+JSZ-1) : Matrix ISZ*JSZ, Results u_{i,j}         *
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! 30/01/2013- Initial Version, by Wei Wang
!**********************************************************************
      SUBROUTINE TDMAIJI_nonCYC(A,B,C,R,BCI,IS,ISZ,JS,JSZ)  !(A,B,C,R,N,UUU,M)
      USE WPRECISION
      IMPLICIT NONE
      
      INTEGER(4),INTENT(IN)  :: JS, JSZ
      INTEGER(4),INTENT(IN)  :: IS, ISZ
      REAL(WP),INTENT(IN)      :: A(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP),INTENT(IN)      :: B(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP),INTENT(IN)      :: C(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP),INTENT(INOUT)   :: R(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP),INTENT(IN)      :: BCI(JS:JSZ+JS-1,2)
      
      REAL(WP)     :: H(IS:ISZ+IS-1,JS:JSZ+JS-1)
      REAL(WP)     :: G(IS:ISZ+IS-1,JS:JSZ+JS-1)
      
      REAL(WP)     :: PP, QQ
      INTEGER(4)  :: I, IEND
      INTEGER(4)  :: J, JEND
      
      H = 0.0_WP
      G = 0.0_WP
      PP= 0.0_WP
      QQ= 0.0_WP
      IEND = IS+ISZ-1
      JEND = JS+JSZ-1

      I = IS    
      DO J=JS,JEND
         H(I,J)=C(I,J)/B(I,J)
         G(I,J)=( R(I,J) - A(I,J)*BCI(J,1) )/B(I,J)
      END DO
   
       
      DO I=IS+1,IEND-1
         DO J=JS,JEND
            PP = 1.0_WP/( B(I,J)-H(I-1,J)*A(I,J) )
            H(I,J) = C(I,J)*PP      
            G(I,J) = ( R(I,J) - A(I,J) * G(I-1,J) ) * PP   
         END DO
      END DO
      
      I = IEND
      DO J=JS, JEND
         PP = 1.0_WP/( B(I,J)-H(I-1,J)*A(I,J) )
         G(I,J)=( R(I,J) - C(I,J)*BCI(J,2) - A(I,J) * G(I-1,J) ) * PP 
         R(I,J)=G(I,J)         
      END DO
      
   
      DO I=IEND-1,IS,-1
         DO J=JS,JEND
            R(I,J)=G(I,J)-H(I,J)*R(I+1,J)
         END DO
      END DO
   
      RETURN
      END SUBROUTINE
