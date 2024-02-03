      SUBROUTINE DORHR_COL( M, N, NB, A, LDA, T, LDT, D, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDT, M, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), D( * ), T( LDT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, JBTEMP1, JBTEMP2, JNB, NPLUSONE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAORHR_COL_GETRFNP, DSCAL, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( NB.LT.1 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDT.LT.MAX( 1, MIN( NB, N ) ) ) THEN
         INFO = -7
      END IF
*
*     Handle error in the input parameters.
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORHR_COL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N ).EQ.0 ) THEN
         RETURN
      END IF
*
*     On input, the M-by-N matrix A contains the orthogonal
*     M-by-N matrix Q_in.
*
*     (1) Compute the unit lower-trapezoidal V (ones on the diagonal
*     are not stored) by performing the "modified" LU-decomposition.
*
*     Q_in - ( S ) = V * U = ( V1 ) * U,
*            ( 0 )           ( V2 )
*
*     where 0 is an (M-N)-by-N zero matrix.
*
*     (1-1) Factor V1 and U.

      CALL DLAORHR_COL_GETRFNP( N, N, A, LDA, D, IINFO )
*
*     (1-2) Solve for V2.
*
      IF( M.GT.N ) THEN
         CALL DTRSM( 'R', 'U', 'N', 'N', M-N, N, ONE, A, LDA, A( N+1, 1 ), LDA )
      END IF
*
*     (2) Reconstruct the block reflector T stored in T(1:NB, 1:N)
*     as a sequence of upper-triangular blocks with NB-size column
*     blocking.
*
*     Loop over the column blocks of size NB of the array A(1:M,1:N)
*     and the array T(1:NB,1:N), JB is the column index of a column
*     block, JNB is the column block size at each step JB.
*
      NPLUSONE = N + 1
      DO JB = 1, N, NB
*
*        (2-0) Determine the column block size JNB.
*
         JNB = MIN( NPLUSONE-JB, NB )
*
*        (2-1) Copy the upper-triangular part of the current JNB-by-JNB
*        diagonal block U(JB) (of the N-by-N matrix U) stored
*        in A(JB:JB+JNB-1,JB:JB+JNB-1) into the upper-triangular part
*        of the current JNB-by-JNB block T(1:JNB,JB:JB+JNB-1)
*        column-by-column, total JNB*(JNB+1)/2 elements.
*
         JBTEMP1 = JB - 1
         DO J = JB, JB+JNB-1
            CALL DCOPY( J-JBTEMP1, A( JB, J ), 1, T( 1, J ), 1 )
         END DO
*
*        (2-2) Perform on the upper-triangular part of the current
*        JNB-by-JNB diagonal block U(JB) (of the N-by-N matrix U) stored
*        in T(1:JNB,JB:JB+JNB-1) the following operation in place:
*        (-1)*U(JB)*S(JB), i.e the result will be stored in the upper-
*        triangular part of T(1:JNB,JB:JB+JNB-1). This multiplication
*        of the JNB-by-JNB diagonal block U(JB) by the JNB-by-JNB
*        diagonal block S(JB) of the N-by-N sign matrix S from the
*        right means changing the sign of each J-th column of the block
*        U(JB) according to the sign of the diagonal element of the block
*        S(JB), i.e. S(J,J) that is stored in the array element D(J).
*
         DO J = JB, JB+JNB-1
            IF( D( J ).EQ.ONE ) THEN
               CALL DSCAL( J-JBTEMP1, -ONE, T( 1, J ), 1 )
            END IF
         END DO
*
*        (2-3) Perform the triangular solve for the current block
*        matrix X(JB):
*
*               X(JB) * (A(JB)**T) = B(JB), where:
*
*               A(JB)**T  is a JNB-by-JNB unit upper-triangular
*                         coefficient block, and A(JB)=V1(JB), which
*                         is a JNB-by-JNB unit lower-triangular block
*                         stored in A(JB:JB+JNB-1,JB:JB+JNB-1).
*                         The N-by-N matrix V1 is the upper part
*                         of the M-by-N lower-trapezoidal matrix V
*                         stored in A(1:M,1:N);
*
*               B(JB)     is a JNB-by-JNB  upper-triangular right-hand
*                         side block, B(JB) = (-1)*U(JB)*S(JB), and
*                         B(JB) is stored in T(1:JNB,JB:JB+JNB-1);
*
*               X(JB)     is a JNB-by-JNB upper-triangular solution
*                         block, X(JB) is the upper-triangular block
*                         reflector T(JB), and X(JB) is stored
*                         in T(1:JNB,JB:JB+JNB-1).
*
*             In other words, we perform the triangular solve for the
*             upper-triangular block T(JB):
*
*               T(JB) * (V1(JB)**T) = (-1)*U(JB)*S(JB).
*
*             Even though the blocks X(JB) and B(JB) are upper-
*             triangular, the routine DTRSM will access all JNB**2
*             elements of the square T(1:JNB,JB:JB+JNB-1). Therefore,
*             we need to set to zero the elements of the block
*             T(1:JNB,JB:JB+JNB-1) below the diagonal before the call
*             to DTRSM.
*
*        (2-3a) Set the elements to zero.
*
         JBTEMP2 = JB - 2
         DO J = JB, JB+JNB-2
            DO I = J-JBTEMP2, NB
               T( I, J ) = ZERO
            END DO
         END DO
*
*        (2-3b) Perform the triangular solve.
*
         CALL DTRSM( 'R', 'L', 'T', 'U', JNB, JNB, ONE, A( JB, JB ), LDA, T( 1, JB ), LDT )
*
      END DO
*
      RETURN
*
*     End of DORHR_COL
*
      END
