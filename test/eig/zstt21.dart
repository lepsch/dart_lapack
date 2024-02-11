      void zstt21(final int N, final int KBAND, final int AD, final int AE, final int SD, final int SE, final Matrix<double> U, final int LDU, final Array<double> _WORK, final Array<double> RWORK, final int RESULT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KBAND, LDU, N;
      double             AD( * ), AE( * ), RESULT( 2 ), RWORK( * ), SD( * ), SE( * );
      Complex         U( LDU, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                J;
      double             ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHER, ZHER2, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN

      // 1)      Constants

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0) return;

      UNFL = dlamch( 'Safe minimum' );
      ULP = dlamch( 'Precision' );

      // Do Test 1

      // Copy A & Compute its 1-Norm:

      zlaset('Full', N, N, CZERO, CZERO, WORK, N );

      ANORM = ZERO;
      TEMP1 = ZERO;

      for (J = 1; J <= N - 1; J++) { // 10
         WORK[( N+1 )*( J-1 )+1] = AD( J );
         WORK[( N+1 )*( J-1 )+2] = AE( J );
         TEMP2 = ( AE( J ) ).abs();
         ANORM = max( ANORM, ( AD( J ) ).abs()+TEMP1+TEMP2 );
         TEMP1 = TEMP2;
      } // 10

      WORK[N**2] = AD( N );
      ANORM = max( ANORM, ( AD( N ) ).abs()+TEMP1, UNFL );

      // Norm of A - USU*

      for (J = 1; J <= N; J++) { // 20
         zher('L', N, -SD( J ), U( 1, J ), 1, WORK, N );
      } // 20

      if ( N > 1 && KBAND == 1 ) {
         for (J = 1; J <= N - 1; J++) { // 30
            zher2('L', N, -DCMPLX( SE( J ) ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N );
         } // 30
      }

      WNORM = ZLANHE( '1', 'L', N, WORK, N, RWORK );

      if ( ANORM > WNORM ) {
         RESULT[1] = ( WNORM / ANORM ) / ( N*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
         } else {
            RESULT[1] = min( WNORM / ANORM, N.toDouble() ) / ( N*ULP );
         }
      }

      // Do Test 2

      // Compute  U U**H - I

      zgemm('N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N );

      for (J = 1; J <= N; J++) { // 40
         WORK[( N+1 )*( J-1 )+1] = WORK( ( N+1 )*( J-1 )+1 ) - CONE;
      } // 40

      RESULT[2] = min( N.toDouble(), ZLANGE( '1', N, N, WORK, N, RWORK ) ) / ( N*ULP );

      }
