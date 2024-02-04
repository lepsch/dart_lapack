      void dstt22(N, M, KBAND, AD, AE, SD, SE, U, LDU, WORK, LDWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KBAND, LDU, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             AD( * ), AE( * ), RESULT( 2 ), SD( * ), SE( * ), U( LDU, * ), WORK( LDWORK, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ANORM, AUKJ, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0 || M <= 0) return;

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Epsilon' );

      // Do Test 1

      // Compute the 1-norm of A.

      if ( N > 1 ) {
         ANORM = ( AD( 1 ) ).abs() + ( AE( 1 ) ).abs();
         for (J = 2; J <= N - 1; J++) { // 10
            ANORM = max( ANORM, ( AD( J ) ).abs()+( AE( J ) ).abs()+ ( AE( J-1 ) ) ).abs();
         } // 10
         ANORM = max( ANORM, ( AD( N ) ).abs()+( AE( N-1 ) ) ).abs();
      } else {
         ANORM = ( AD( 1 ) ).abs();
      }
      ANORM = max( ANORM, UNFL );

      // Norm of U'AU - S

      for (I = 1; I <= M; I++) { // 40
         for (J = 1; J <= M; J++) { // 30
            WORK[I, J] = ZERO;
            for (K = 1; K <= N; K++) { // 20
               AUKJ = AD( K )*U( K, J );
               if (K != N) AUKJ = AUKJ + AE( K )*U( K+1, J );
               IF( K != 1 ) AUKJ = AUKJ + AE( K-1 )*U( K-1, J );
               WORK[I, J] = WORK( I, J ) + U( K, I )*AUKJ;
            } // 20
         } // 30
         WORK[I, I] = WORK( I, I ) - SD( I );
         if ( KBAND == 1 ) {
            if (I != 1) WORK( I, I-1 ) = WORK( I, I-1 ) - SE( I-1 );
            IF( I != N ) WORK( I, I+1 ) = WORK( I, I+1 ) - SE( I );
         }
      } // 40

      WNORM = DLANSY( '1', 'L', M, WORK, M, WORK( 1, M+1 ) );

      if ( ANORM > WNORM ) {
         RESULT[1] = ( WNORM / ANORM ) / ( M*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ( min( WNORM, M*ANORM ) / ANORM ) / ( M*ULP );
         } else {
            RESULT[1] = min( WNORM / ANORM, M.toDouble() ) / ( M*ULP );
         }
      }

      // Do Test 2

      // Compute  U'U - I

      dgemm('T', 'N', M, M, N, ONE, U, LDU, U, LDU, ZERO, WORK, M );

      for (J = 1; J <= M; J++) { // 50
         WORK[J, J] = WORK( J, J ) - ONE;
      } // 50

      RESULT[2] = min( M.toDouble(), DLANGE( '1', M, M, WORK, M, WORK( 1, M+1 ) ) ) / ( M*ULP );

      return;
      }
