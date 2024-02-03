      void zstt22(N, M, KBAND, AD, AE, SD, SE, U, LDU, WORK, LDWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KBAND, LDU, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             AD( * ), AE( * ), RESULT( 2 ), RWORK( * ), SD( * ), SE( * );
      COMPLEX*16         U( LDU, * ), WORK( LDWORK, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ANORM, ULP, UNFL, WNORM;
      COMPLEX*16         AUKJ;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANSY;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO;
      RESULT( 2 ) = ZERO;
      if (N <= 0 || M <= 0) return;

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Epsilon' );

      // Do Test 1

      // Compute the 1-norm of A.

      if ( N > 1 ) {
         ANORM = ABS( AD( 1 ) ) + ABS( AE( 1 ) );
         for (J = 2; J <= N - 1; J++) { // 10
            ANORM = max( ANORM, ABS( AD( J ) )+ABS( AE( J ) )+ ABS( AE( J-1 ) ) );
         } // 10
         ANORM = max( ANORM, ABS( AD( N ) )+ABS( AE( N-1 ) ) );
      } else {
         ANORM = ABS( AD( 1 ) );
      }
      ANORM = max( ANORM, UNFL );

      // Norm of U*AU - S

      for (I = 1; I <= M; I++) { // 40
         for (J = 1; J <= M; J++) { // 30
            WORK( I, J ) = CZERO;
            for (K = 1; K <= N; K++) { // 20
               AUKJ = AD( K )*U( K, J );
               if (K != N) AUKJ = AUKJ + AE( K )*U( K+1, J );
               IF( K != 1 ) AUKJ = AUKJ + AE( K-1 )*U( K-1, J );
               WORK( I, J ) = WORK( I, J ) + U( K, I )*AUKJ;
            } // 20
         } // 30
         WORK( I, I ) = WORK( I, I ) - SD( I );
         if ( KBAND == 1 ) {
            if (I != 1) WORK( I, I-1 ) = WORK( I, I-1 ) - SE( I-1 );
            IF( I != N ) WORK( I, I+1 ) = WORK( I, I+1 ) - SE( I );
         }
      } // 40

      WNORM = ZLANSY( '1', 'L', M, WORK, M, RWORK );

      if ( ANORM > WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( M*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT( 1 ) = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*ULP );
         } else {
            RESULT( 1 ) = min( WNORM / ANORM, DBLE( M ) ) / ( M*ULP );
         }
      }

      // Do Test 2

      // Compute  U*U - I

      zgemm('T', 'N', M, M, N, CONE, U, LDU, U, LDU, CZERO, WORK, M );

      for (J = 1; J <= M; J++) { // 50
         WORK( J, J ) = WORK( J, J ) - ONE;
      } // 50

      RESULT( 2 ) = min( DBLE( M ), ZLANGE( '1', M, M, WORK, M, RWORK ) ) / ( M*ULP );

      return;

      // End of ZSTT22

      }
