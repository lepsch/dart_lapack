      void cstt22(final int N, final int M, final int KBAND, final int AD, final int AE, final int SD, final int SE, final Matrix<double> U, final int LDU, final Matrix<double> WORK, final int LDWORK, final Array<double> RWORK, final int RESULT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KBAND, LDU, LDWORK, M, N;
      double               AD( * ), AE( * ), RESULT( 2 ), RWORK( * ), SD( * ), SE( * );
      Complex            U( LDU, * ), WORK( LDWORK, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, J, K;
      double               ANORM, ULP, UNFL, WNORM;
      Complex            AUKJ;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, CLANSY, SLAMCH;
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0 || M <= 0) return;

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' );

      // Do Test 1

      // Compute the 1-norm of A.

      if ( N > 1 ) {
         ANORM = ( AD( 1 ) ).abs() + ( AE( 1 ) ).abs();
         for (J = 2; J <= N - 1; J++) { // 10
            ANORM = max( ANORM, ( AD( J ) ).abs()+( AE( J ) ).abs()+ ( AE( J-1 ) ).abs() );
         } // 10
         ANORM = max( ANORM, ( AD( N ) ).abs()+( AE( N-1 ) ).abs() );
      } else {
         ANORM = ( AD( 1 ) ).abs();
      }
      ANORM = max( ANORM, UNFL );

      // Norm of U*AU - S

      for (I = 1; I <= M; I++) { // 40
         for (J = 1; J <= M; J++) { // 30
            WORK[I][J] = CZERO;
            for (K = 1; K <= N; K++) { // 20
               AUKJ = AD( K )*U( K, J );
               if (K != N) AUKJ = AUKJ + AE( K )*U( K+1, J );
               IF( K != 1 ) AUKJ = AUKJ + AE( K-1 )*U( K-1, J );
               WORK[I][J] = WORK( I, J ) + U( K, I )*AUKJ;
            } // 20
         } // 30
         WORK[I][I] = WORK( I, I ) - SD( I );
         if ( KBAND == 1 ) {
            if (I != 1) WORK( I, I-1 ) = WORK( I, I-1 ) - SE( I-1 );
            IF( I != N ) WORK( I, I+1 ) = WORK( I, I+1 ) - SE( I );
         }
      } // 40

      WNORM = CLANSY( '1', 'L', M, WORK, M, RWORK );

      if ( ANORM > WNORM ) {
         RESULT[1] = ( WNORM / ANORM ) / ( M*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*ULP );
         } else {
            RESULT[1] = min( WNORM / ANORM, REAL( M ) ) / ( M*ULP );
         }
      }

      // Do Test 2

      // Compute  U*U - I

      cgemm('T', 'N', M, M, N, CONE, U, LDU, U, LDU, CZERO, WORK, M );

      for (J = 1; J <= M; J++) { // 50
         WORK[J][J] = WORK( J, J ) - ONE;
      } // 50

      RESULT[2] = min( double( M ), CLANGE( '1', M, M, WORK, M, RWORK ) ) / ( M*ULP );

      }
