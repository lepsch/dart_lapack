      SUBROUTINE DLAKF2( M, N, A, LDA, B, D, E, Z, LDZ )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDZ, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDA, * ), D( LDA, * ), E( LDA, * ), Z( LDZ, * )
*     ..
*
*  ====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, IK, J, JK, L, MN, MN2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASET
*     ..
*     .. Executable Statements ..
*
*     Initialize Z
*
      MN = M*N
      MN2 = 2*MN
      CALL DLASET( 'Full', MN2, MN2, ZERO, ZERO, Z, LDZ )
*
      IK = 1
      DO 50 L = 1, N
*
*        form kron(In, A)
*
         DO 20 I = 1, M
            DO 10 J = 1, M
               Z( IK+I-1, IK+J-1 ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
*
*        form kron(In, D)
*
         DO 40 I = 1, M
            DO 30 J = 1, M
               Z( IK+MN+I-1, IK+J-1 ) = D( I, J )
   30       CONTINUE
   40    CONTINUE
*
         IK = IK + M
   50 CONTINUE
*
      IK = 1
      DO 90 L = 1, N
         JK = MN + 1
*
         DO 80 J = 1, N
*
*           form -kron(B', Im)
*
            DO 60 I = 1, M
               Z( IK+I-1, JK+I-1 ) = -B( J, L )
   60       CONTINUE
*
*           form -kron(E', Im)
*
            DO 70 I = 1, M
               Z( IK+MN+I-1, JK+I-1 ) = -E( J, L )
   70       CONTINUE
*
            JK = JK + M
   80    CONTINUE
*
         IK = IK + M
   90 CONTINUE
*
      RETURN
*
*     End of DLAKF2
*
      END
