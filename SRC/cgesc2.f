      SUBROUTINE CGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, N;
      REAL               SCALE
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      COMPLEX            A( LDA, * ), RHS( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               BIGNUM, EPS, SMLNUM
      COMPLEX            TEMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASWP, CSCAL
      // ..
      // .. External Functions ..
      int                ICAMAX;
      REAL               SLAMCH
      // EXTERNAL ICAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, REAL
      // ..
      // .. Executable Statements ..

      // Set constant to control overflow

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // Apply permutations IPIV to RHS

      CALL CLASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )

      // Solve for L part

      DO 20 I = 1, N - 1
         DO 10 J = I + 1, N
            RHS( J ) = RHS( J ) - A( J, I )*RHS( I )
   10    CONTINUE
   20 CONTINUE

      // Solve for U part

      SCALE = ONE

      // Check for scaling

      I = ICAMAX( N, RHS, 1 )
      IF( TWO*SMLNUM*ABS( RHS( I ) ).GT.ABS( A( N, N ) ) ) THEN
         TEMP = CMPLX( ONE / TWO, ZERO ) / ABS( RHS( I ) )
         CALL CSCAL( N, TEMP, RHS( 1 ), 1 )
         SCALE = SCALE*REAL( TEMP )
      END IF
      DO 40 I = N, 1, -1
         TEMP = CMPLX( ONE, ZERO ) / A( I, I )
         RHS( I ) = RHS( I )*TEMP
         DO 30 J = I + 1, N
            RHS( I ) = RHS( I ) - RHS( J )*( A( I, J )*TEMP )
   30    CONTINUE
   40 CONTINUE

      // Apply permutations JPIV to the solution (RHS)

      CALL CLASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
      RETURN

      // End of CGESC2

      END
