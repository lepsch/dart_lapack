      void cptts2(IUPLO, N, NRHS, D, E, B, LDB ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IUPLO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               D( * );
      COMPLEX            B( LDB, * ), E( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 1 ) {
         if (N == 1) csscal( NRHS, 1. / D( 1 ), B, LDB );
         return;
      }

      if ( IUPLO == 1 ) {

         // Solve A * X = B using the factorization A = U**H *D*U,
         // overwriting each right hand side vector with its solution.

         if ( NRHS <= 2 ) {
            J = 1;
            } // 5

            // Solve U**H * x = b.

            for (I = 2; I <= N; I++) { // 10
               B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) );
            } // 10

            // Solve D * U * x = b.

            for (I = 1; I <= N; I++) { // 20
               B( I, J ) = B( I, J ) / D( I );
            } // 20
            for (I = N - 1; I >= 1; I--) { // 30
               B( I, J ) = B( I, J ) - B( I+1, J )*E( I );
            } // 30
            if ( J < NRHS ) {
               J = J + 1;
               GO TO 5;
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 60

               // Solve U**H * x = b.

               for (I = 2; I <= N; I++) { // 40
                  B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) );
               } // 40

               // Solve D * U * x = b.

               B( N, J ) = B( N, J ) / D( N );
               for (I = N - 1; I >= 1; I--) { // 50
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I );
               } // 50
            } // 60
         }
      } else {

         // Solve A * X = B using the factorization A = L*D*L**H,
         // overwriting each right hand side vector with its solution.

         if ( NRHS <= 2 ) {
            J = 1;
            } // 65

            // Solve L * x = b.

            for (I = 2; I <= N; I++) { // 70
               B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 );
            } // 70

            // Solve D * L**H * x = b.

            for (I = 1; I <= N; I++) { // 80
               B( I, J ) = B( I, J ) / D( I );
            } // 80
            for (I = N - 1; I >= 1; I--) { // 90
               B( I, J ) = B( I, J ) - B( I+1, J )*CONJG( E( I ) );
            } // 90
            if ( J < NRHS ) {
               J = J + 1;
               GO TO 65;
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 120

               // Solve L * x = b.

               for (I = 2; I <= N; I++) { // 100
                  B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 );
               } // 100

               // Solve D * L**H * x = b.

               B( N, J ) = B( N, J ) / D( N );
               for (I = N - 1; I >= 1; I--) { // 110
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*CONJG( E( I ) );
               } // 110
            } // 120
         }
      }

      return;
      }
