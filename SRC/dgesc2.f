      SUBROUTINE DGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, N
      double             SCALE;
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), JPIV( * )
      double             A( LDA, * ), RHS( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE, TWO;
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, J
      double             BIGNUM, EPS, SMLNUM, TEMP;
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DSCAL
*     ..
*     .. External Functions ..
      int                IDAMAX
      double             DLAMCH;
      EXTERNAL           IDAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*      Set constant to control overflow
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Apply permutations IPIV to RHS
*
      CALL DLASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )
*
*     Solve for L part
*
      DO 20 I = 1, N - 1
         DO 10 J = I + 1, N
            RHS( J ) = RHS( J ) - A( J, I )*RHS( I )
   10    CONTINUE
   20 CONTINUE
*
*     Solve for U part
*
      SCALE = ONE
*
*     Check for scaling
*
      I = IDAMAX( N, RHS, 1 )
      IF( TWO*SMLNUM*ABS( RHS( I ) ).GT.ABS( A( N, N ) ) ) THEN
         TEMP = ( ONE / TWO ) / ABS( RHS( I ) )
         CALL DSCAL( N, TEMP, RHS( 1 ), 1 )
         SCALE = SCALE*TEMP
      END IF
*
      DO 40 I = N, 1, -1
         TEMP = ONE / A( I, I )
         RHS( I ) = RHS( I )*TEMP
         DO 30 J = I + 1, N
            RHS( I ) = RHS( I ) - RHS( J )*( A( I, J )*TEMP )
   30    CONTINUE
   40 CONTINUE
*
*     Apply permutations JPIV to the solution (RHS)
*
      CALL DLASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
      RETURN
*
*     End of DGESC2
*
      END
