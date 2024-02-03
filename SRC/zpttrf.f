      void zpttrf(N, D, E, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             D( * );
      COMPLEX*16         E( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, I4;
      double             EII, EIR, F, G;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('ZPTTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Compute the L*D*L**H (or U**H *D*U) factorization of A.

      I4 = MOD( N-1, 4 );
      for (I = 1; I <= I4; I++) { // 10
         if ( D( I ) <= ZERO ) {
            INFO = I;
            GO TO 30;
         }
         EIR = DBLE( E( I ) );
         EII = DIMAG( E( I ) );
         F = EIR / D( I );
         G = EII / D( I );
         E( I ) = DCMPLX( F, G );
         D( I+1 ) = D( I+1 ) - F*EIR - G*EII;
      } // 10

      DO 20 I = I4 + 1, N - 4, 4;

         // Drop out of the loop if d(i) <= 0: the matrix is not positive
         // definite.

         if ( D( I ) <= ZERO ) {
            INFO = I;
            GO TO 30;
         }

         // Solve for e(i) and d(i+1).

         EIR = DBLE( E( I ) );
         EII = DIMAG( E( I ) );
         F = EIR / D( I );
         G = EII / D( I );
         E( I ) = DCMPLX( F, G );
         D( I+1 ) = D( I+1 ) - F*EIR - G*EII;

         if ( D( I+1 ) <= ZERO ) {
            INFO = I + 1;
            GO TO 30;
         }

         // Solve for e(i+1) and d(i+2).

         EIR = DBLE( E( I+1 ) );
         EII = DIMAG( E( I+1 ) );
         F = EIR / D( I+1 );
         G = EII / D( I+1 );
         E( I+1 ) = DCMPLX( F, G );
         D( I+2 ) = D( I+2 ) - F*EIR - G*EII;

         if ( D( I+2 ) <= ZERO ) {
            INFO = I + 2;
            GO TO 30;
         }

         // Solve for e(i+2) and d(i+3).

         EIR = DBLE( E( I+2 ) );
         EII = DIMAG( E( I+2 ) );
         F = EIR / D( I+2 );
         G = EII / D( I+2 );
         E( I+2 ) = DCMPLX( F, G );
         D( I+3 ) = D( I+3 ) - F*EIR - G*EII;

         if ( D( I+3 ) <= ZERO ) {
            INFO = I + 3;
            GO TO 30;
         }

         // Solve for e(i+3) and d(i+4).

         EIR = DBLE( E( I+3 ) );
         EII = DIMAG( E( I+3 ) );
         F = EIR / D( I+3 );
         G = EII / D( I+3 );
         E( I+3 ) = DCMPLX( F, G );
         D( I+4 ) = D( I+4 ) - F*EIR - G*EII;
      } // 20

      // Check d(n) for positive definiteness.

      if( D( N ) <= ZERO ) INFO = N;

      } // 30
      return;

      // End of ZPTTRF

      }
