      void zptts2(final int IUPLO, final int N, final int NRHS, final int D, final int E, final int B, final int LDB) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IUPLO, LDB, N, NRHS;
      double             D( * );
      Complex         B( LDB, * ), E( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG

      // Quick return if possible

      if ( N <= 1 ) {
         if (N == 1) zdscal( NRHS, 1.0 / D( 1 ), B, LDB );
         return;
      }

      if ( IUPLO == 1 ) {

         // Solve A * X = B using the factorization A = U**H *D*U,
         // overwriting each right hand side vector with its solution.

         if ( NRHS <= 2 ) {
            J = 1;
            } // 10

            // Solve U**H * x = b.

            for (I = 2; I <= N; I++) { // 20
               B[I][J] = B( I, J ) - B( I-1, J )*DCONJG( E( I-1 ) );
            } // 20

            // Solve D * U * x = b.

            for (I = 1; I <= N; I++) { // 30
               B[I][J] = B( I, J ) / D( I );
            } // 30
            for (I = N - 1; I >= 1; I--) { // 40
               B[I][J] = B( I, J ) - B( I+1, J )*E( I );
            } // 40
            if ( J < NRHS ) {
               J = J + 1;
               GO TO 10;
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 70

               // Solve U**H * x = b.

               for (I = 2; I <= N; I++) { // 50
                  B[I][J] = B( I, J ) - B( I-1, J )*DCONJG( E( I-1 ) );
               } // 50

               // Solve D * U * x = b.

               B[N][J] = B( N, J ) / D( N );
               for (I = N - 1; I >= 1; I--) { // 60
                  B[I][J] = B( I, J ) / D( I ) - B( I+1, J )*E( I );
               } // 60
            } // 70
         }
      } else {

         // Solve A * X = B using the factorization A = L*D*L**H,
         // overwriting each right hand side vector with its solution.

         if ( NRHS <= 2 ) {
            J = 1;
            } // 80

            // Solve L * x = b.

            for (I = 2; I <= N; I++) { // 90
               B[I][J] = B( I, J ) - B( I-1, J )*E( I-1 );
            } // 90

            // Solve D * L**H * x = b.

            for (I = 1; I <= N; I++) { // 100
               B[I][J] = B( I, J ) / D( I );
            } // 100
            for (I = N - 1; I >= 1; I--) { // 110
               B[I][J] = B( I, J ) - B( I+1, J )*DCONJG( E( I ) );
            } // 110
            if ( J < NRHS ) {
               J = J + 1;
               GO TO 80;
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 140

               // Solve L * x = b.

               for (I = 2; I <= N; I++) { // 120
                  B[I][J] = B( I, J ) - B( I-1, J )*E( I-1 );
               } // 120

               // Solve D * L**H * x = b.

               B[N][J] = B( N, J ) / D( N );
               for (I = N - 1; I >= 1; I--) { // 130
                  B[I][J] = B( I, J ) / D( I ) - B( I+1, J )*DCONJG( E( I ) );
               } // 130
            } // 140
         }
      }

      }
