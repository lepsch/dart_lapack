      void claptm(UPLO, N, NRHS, ALPHA, D, E, X, LDX, BETA, B, LDB ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDB, LDX, N, NRHS;
      double               ALPHA, BETA;
      double               D( * );
      Complex            B( LDB, * ), E( * ), X( LDX, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG

      if (N == 0) return;

      if ( BETA == ZERO ) {
         for (J = 1; J <= NRHS; J++) { // 20
            for (I = 1; I <= N; I++) { // 10
               B[I][J] = ZERO;
            } // 10
         } // 20
      } else if ( BETA == -ONE ) {
         for (J = 1; J <= NRHS; J++) { // 40
            for (I = 1; I <= N; I++) { // 30
               B[I][J] = -B( I, J );
            } // 30
         } // 40
      }

      if ( ALPHA == ONE ) {
         if ( lsame( UPLO, 'U' ) ) {

            // Compute B := B + A*X, where E is the superdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 60
               if ( N == 1 ) {
                  B[1][J] = B( 1, J ) + D( 1 )*X( 1, J );
               } else {
                  B[1][J] = B( 1, J ) + D( 1 )*X( 1, J ) + E( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) + CONJG( E( N-1 ) )* X( N-1, J ) + D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 50
                     B[I][J] = B( I, J ) + CONJG( E( I-1 ) )* X( I-1, J ) + D( I )*X( I, J ) + E( I )*X( I+1, J );
                  } // 50
               }
            } // 60
         } else {

            // Compute B := B + A*X, where E is the subdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 80
               if ( N == 1 ) {
                  B[1][J] = B( 1, J ) + D( 1 )*X( 1, J );
               } else {
                  B[1][J] = B( 1, J ) + D( 1 )*X( 1, J ) + CONJG( E( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) + E( N-1 )*X( N-1, J ) + D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 70
                     B[I][J] = B( I, J ) + E( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + CONJG( E( I ) )*X( I+1, J );
                  } // 70
               }
            } // 80
         }
      } else if ( ALPHA == -ONE ) {
         if ( lsame( UPLO, 'U' ) ) {

            // Compute B := B - A*X, where E is the superdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 100
               if ( N == 1 ) {
                  B[1][J] = B( 1, J ) - D( 1 )*X( 1, J );
               } else {
                  B[1][J] = B( 1, J ) - D( 1 )*X( 1, J ) - E( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) - CONJG( E( N-1 ) )* X( N-1, J ) - D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 90
                     B[I][J] = B( I, J ) - CONJG( E( I-1 ) )* X( I-1, J ) - D( I )*X( I, J ) - E( I )*X( I+1, J );
                  } // 90
               }
            } // 100
         } else {

            // Compute B := B - A*X, where E is the subdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 120
               if ( N == 1 ) {
                  B[1][J] = B( 1, J ) - D( 1 )*X( 1, J );
               } else {
                  B[1][J] = B( 1, J ) - D( 1 )*X( 1, J ) - CONJG( E( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) - E( N-1 )*X( N-1, J ) - D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 110
                     B[I][J] = B( I, J ) - E( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - CONJG( E( I ) )*X( I+1, J );
                  } // 110
               }
            } // 120
         }
      }
      }
