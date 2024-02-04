      void spttrf(N, D, E, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double               D( * ), E( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, I4;
      double               EI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('SPTTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Compute the L*D*L**T (or U**T*D*U) factorization of A.

      I4 = (N-1 % 4);
      for (I = 1; I <= I4; I++) { // 10
         if ( D( I ) <= ZERO ) {
            INFO = I;
            GO TO 30;
         }
         EI = E( I );
         E[I] = EI / D( I );
         D[I+1] = D( I+1 ) - E( I )*EI;
      } // 10

      for (I = I4 + 1; 4 < 0 ? I >= N - 4 : I <= N - 4; I += 4) { // 20

         // Drop out of the loop if d(i) <= 0: the matrix is not positive
         // definite.

         if ( D( I ) <= ZERO ) {
            INFO = I;
            GO TO 30;
         }

         // Solve for e(i) and d(i+1).

         EI = E( I );
         E[I] = EI / D( I );
         D[I+1] = D( I+1 ) - E( I )*EI;

         if ( D( I+1 ) <= ZERO ) {
            INFO = I + 1;
            GO TO 30;
         }

         // Solve for e(i+1) and d(i+2).

         EI = E( I+1 );
         E[I+1] = EI / D( I+1 );
         D[I+2] = D( I+2 ) - E( I+1 )*EI;

         if ( D( I+2 ) <= ZERO ) {
            INFO = I + 2;
            GO TO 30;
         }

         // Solve for e(i+2) and d(i+3).

         EI = E( I+2 );
         E[I+2] = EI / D( I+2 );
         D[I+3] = D( I+3 ) - E( I+2 )*EI;

         if ( D( I+3 ) <= ZERO ) {
            INFO = I + 3;
            GO TO 30;
         }

         // Solve for e(i+3) and d(i+4).

         EI = E( I+3 );
         E[I+3] = EI / D( I+3 );
         D[I+4] = D( I+4 ) - E( I+3 )*EI;
      } // 20

      // Check d(n) for positive definiteness.

      if( D( N ) <= ZERO ) INFO = N;

      } // 30
      return;
      }