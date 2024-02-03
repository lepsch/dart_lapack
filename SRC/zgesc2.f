      SUBROUTINE ZGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      COMPLEX*16         A( LDA, * ), RHS( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             BIGNUM, EPS, SMLNUM;
      COMPLEX*16         TEMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASWP, ZSCAL
      // ..
      // .. External Functions ..
      int                IZAMAX;
      double             DLAMCH;
      // EXTERNAL IZAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      // Set constant to control overflow

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // Apply permutations IPIV to RHS

      zlaswp(1, RHS, LDA, 1, N-1, IPIV, 1 );

      // Solve for L part

      for (I = 1; I <= N - 1; I++) { // 20
         for (J = I + 1; J <= N; J++) { // 10
            RHS( J ) = RHS( J ) - A( J, I )*RHS( I )
         } // 10
      } // 20

      // Solve for U part

      SCALE = ONE

      // Check for scaling

      I = IZAMAX( N, RHS, 1 )
      if ( TWO*SMLNUM*ABS( RHS( I ) ) > ABS( A( N, N ) ) ) {
         TEMP = DCMPLX( ONE / TWO, ZERO ) / ABS( RHS( I ) )
         zscal(N, TEMP, RHS( 1 ), 1 );
         SCALE = SCALE*DBLE( TEMP )
      }
      DO 40 I = N, 1, -1
         TEMP = DCMPLX( ONE, ZERO ) / A( I, I )
         RHS( I ) = RHS( I )*TEMP
         for (J = I + 1; J <= N; J++) { // 30
            RHS( I ) = RHS( I ) - RHS( J )*( A( I, J )*TEMP )
         } // 30
      } // 40

      // Apply permutations JPIV to the solution (RHS)

      zlaswp(1, RHS, LDA, 1, N-1, JPIV, -1 );
      RETURN

      // End of ZGESC2

      }
