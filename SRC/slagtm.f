      SUBROUTINE SLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, B, LDB );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, N, NRHS;
      REAL               ALPHA, BETA;
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      if (N == 0) RETURN;

      // Multiply B by BETA if BETA != 1.

      if ( BETA == ZERO ) {
         for (J = 1; J <= NRHS; J++) { // 20
            for (I = 1; I <= N; I++) { // 10
               B( I, J ) = ZERO;
            } // 10
         } // 20
      } else if ( BETA == -ONE ) {
         for (J = 1; J <= NRHS; J++) { // 40
            for (I = 1; I <= N; I++) { // 30
               B( I, J ) = -B( I, J );
            } // 30
         } // 40
      }

      if ( ALPHA == ONE ) {
         if ( LSAME( TRANS, 'N' ) ) {

            // Compute B := B + A*X

            for (J = 1; J <= NRHS; J++) { // 60
               if ( N == 1 ) {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J );
               } else {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + DU( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) + DL( N-1 )*X( N-1, J ) + D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 50
                     B( I, J ) = B( I, J ) + DL( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + DU( I )*X( I+1, J );
                  } // 50
               }
            } // 60
         } else {

            // Compute B := B + A**T*X

            for (J = 1; J <= NRHS; J++) { // 80
               if ( N == 1 ) {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J );
               } else {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + DL( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) + DU( N-1 )*X( N-1, J ) + D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 70
                     B( I, J ) = B( I, J ) + DU( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + DL( I )*X( I+1, J );
                  } // 70
               }
            } // 80
         }
      } else if ( ALPHA == -ONE ) {
         if ( LSAME( TRANS, 'N' ) ) {

            // Compute B := B - A*X

            for (J = 1; J <= NRHS; J++) { // 100
               if ( N == 1 ) {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J );
               } else {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - DU( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) - DL( N-1 )*X( N-1, J ) - D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 90
                     B( I, J ) = B( I, J ) - DL( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - DU( I )*X( I+1, J );
                  } // 90
               }
            } // 100
         } else {

            // Compute B := B - A**T*X

            for (J = 1; J <= NRHS; J++) { // 120
               if ( N == 1 ) {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J );
               } else {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - DL( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) - DU( N-1 )*X( N-1, J ) - D( N )*X( N, J );
                  for (I = 2; I <= N - 1; I++) { // 110
                     B( I, J ) = B( I, J ) - DU( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - DL( I )*X( I+1, J );
                  } // 110
               }
            } // 120
         }
      }
      RETURN;

      // End of SLAGTM

      }
