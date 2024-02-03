      SUBROUTINE ZLAPTM( UPLO, N, NRHS, ALPHA, D, E, X, LDX, BETA, B, LDB )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDB, LDX, N, NRHS;
      double             ALPHA, BETA;
      // ..
      // .. Array Arguments ..
      double             D( * );
      COMPLEX*16         B( LDB, * ), E( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG
      // ..
      // .. Executable Statements ..

      IF( N.EQ.0 ) RETURN

      if ( BETA.EQ.ZERO ) {
         for (J = 1; J <= NRHS; J++) { // 20
            for (I = 1; I <= N; I++) { // 10
               B( I, J ) = ZERO
            } // 10
         } // 20
      } else if ( BETA.EQ.-ONE ) {
         for (J = 1; J <= NRHS; J++) { // 40
            for (I = 1; I <= N; I++) { // 30
               B( I, J ) = -B( I, J )
            } // 30
         } // 40
      }

      if ( ALPHA.EQ.ONE ) {
         if ( LSAME( UPLO, 'U' ) ) {

            // Compute B := B + A*X, where E is the superdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 60
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + E( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) + DCONJG( E( N-1 ) )* X( N-1, J ) + D( N )*X( N, J )
                  DO 50 I = 2, N - 1
                     B( I, J ) = B( I, J ) + DCONJG( E( I-1 ) )* X( I-1, J ) + D( I )*X( I, J ) + E( I )*X( I+1, J )
                  } // 50
               }
            } // 60
         } else {

            // Compute B := B + A*X, where E is the subdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 80
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + DCONJG( E( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) + E( N-1 )*X( N-1, J ) + D( N )*X( N, J )
                  DO 70 I = 2, N - 1
                     B( I, J ) = B( I, J ) + E( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + DCONJG( E( I ) )*X( I+1, J )
                  } // 70
               }
            } // 80
         }
      } else if ( ALPHA.EQ.-ONE ) {
         if ( LSAME( UPLO, 'U' ) ) {

            // Compute B := B - A*X, where E is the superdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 100
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - E( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) - DCONJG( E( N-1 ) )* X( N-1, J ) - D( N )*X( N, J )
                  DO 90 I = 2, N - 1
                     B( I, J ) = B( I, J ) - DCONJG( E( I-1 ) )* X( I-1, J ) - D( I )*X( I, J ) - E( I )*X( I+1, J )
                  } // 90
               }
            } // 100
         } else {

            // Compute B := B - A*X, where E is the subdiagonal of A.

            for (J = 1; J <= NRHS; J++) { // 120
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - DCONJG( E( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) - E( N-1 )*X( N-1, J ) - D( N )*X( N, J )
                  DO 110 I = 2, N - 1
                     B( I, J ) = B( I, J ) - E( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - DCONJG( E( I ) )*X( I+1, J )
                  } // 110
               }
            } // 120
         }
      }
      RETURN

      // End of ZLAPTM

      }
