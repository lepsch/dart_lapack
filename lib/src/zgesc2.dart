      void zgesc2(N, A, LDA, RHS, IPIV, JPIV, SCALE ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      Complex         A( LDA, * ), RHS( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             BIGNUM, EPS, SMLNUM;
      Complex         TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASWP, ZSCAL
      // ..
      // .. External Functions ..
      //- int                IZAMAX;
      //- double             DLAMCH;
      // EXTERNAL IZAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      // Set constant to control overflow

      EPS = DLAMCH( 'P' );
      SMLNUM = DLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Apply permutations IPIV to RHS

      zlaswp(1, RHS, LDA, 1, N-1, IPIV, 1 );

      // Solve for L part

      for (I = 1; I <= N - 1; I++) { // 20
         for (J = I + 1; J <= N; J++) { // 10
            RHS[J] = RHS( J ) - A( J, I )*RHS( I );
         } // 10
      } // 20

      // Solve for U part

      SCALE = ONE;

      // Check for scaling

      I = IZAMAX( N, RHS, 1 );
      if ( TWO*SMLNUM*( RHS( I ) ).abs() > ( A( N, N ) ) ).abs() {
         TEMP = DCMPLX( ONE / TWO, ZERO ) / ( RHS( I ) ).abs();
         zscal(N, TEMP, RHS( 1 ), 1 );
         SCALE = SCALE*TEMP.toDouble();
      }
      for (I = N; I >= 1; I--) { // 40
         TEMP = DCMPLX( ONE, ZERO ) / A( I, I );
         RHS[I] = RHS( I )*TEMP;
         for (J = I + 1; J <= N; J++) { // 30
            RHS[I] = RHS( I ) - RHS( J )*( A( I, J )*TEMP );
         } // 30
      } // 40

      // Apply permutations JPIV to the solution (RHS)

      zlaswp(1, RHS, LDA, 1, N-1, JPIV, -1 );
      return;
      }