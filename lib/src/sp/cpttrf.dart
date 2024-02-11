      void cpttrf(final int N, final int D, final int E, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, N;
      double               D( * );
      Complex            E( * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int                I, I4;
      double               EII, EIR, F, G;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, MOD, REAL

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('CPTTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Compute the L*D*L**H (or U**H *D*U) factorization of A.

      I4 = (N-1 % 4);
      for (I = 1; I <= I4; I++) { // 10
         if ( D( I ) <= ZERO ) {
            INFO = I;
            GO TO 20;
         }
         EIR = double( E( I ) );
         EII = AIMAG( E( I ) );
         F = EIR / D( I );
         G = EII / D( I );
         E[I] = CMPLX( F, G );
         D[I+1] = D( I+1 ) - F*EIR - G*EII;
      } // 10

      for (I = I4+1; 4 < 0 ? I >= N - 4 : I <= N - 4; I += 4) { // 110

         // Drop out of the loop if d(i) <= 0: the matrix is not positive
         // definite.

         if ( D( I ) <= ZERO ) {
            INFO = I;
            GO TO 20;
         }

         // Solve for e(i) and d(i+1).

         EIR = double( E( I ) );
         EII = AIMAG( E( I ) );
         F = EIR / D( I );
         G = EII / D( I );
         E[I] = CMPLX( F, G );
         D[I+1] = D( I+1 ) - F*EIR - G*EII;

         if ( D( I+1 ) <= ZERO ) {
            INFO = I+1;
            GO TO 20;
         }

         // Solve for e(i+1) and d(i+2).

         EIR = double( E( I+1 ) );
         EII = AIMAG( E( I+1 ) );
         F = EIR / D( I+1 );
         G = EII / D( I+1 );
         E[I+1] = CMPLX( F, G );
         D[I+2] = D( I+2 ) - F*EIR - G*EII;

         if ( D( I+2 ) <= ZERO ) {
            INFO = I+2;
            GO TO 20;
         }

         // Solve for e(i+2) and d(i+3).

         EIR = double( E( I+2 ) );
         EII = AIMAG( E( I+2 ) );
         F = EIR / D( I+2 );
         G = EII / D( I+2 );
         E[I+2] = CMPLX( F, G );
         D[I+3] = D( I+3 ) - F*EIR - G*EII;

         if ( D( I+3 ) <= ZERO ) {
            INFO = I+3;
            GO TO 20;
         }

         // Solve for e(i+3) and d(i+4).

         EIR = double( E( I+3 ) );
         EII = AIMAG( E( I+3 ) );
         F = EIR / D( I+3 );
         G = EII / D( I+3 );
         E[I+3] = CMPLX( F, G );
         D[I+4] = D( I+4 ) - F*EIR - G*EII;
      } // 110

      // Check d(n) for positive definiteness.

      if( D( N ) <= ZERO ) INFO = N;

      } // 20
      }
