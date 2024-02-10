      void sstt21(N, KBAND, AD, AE, SD, SE, final Matrix<double> U, final int LDU, final Array<double> _WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KBAND, LDU, N;
      double               AD( * ), AE( * ), RESULT( 2 ), SD( * ), SE( * ), U( LDU, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLASET, SSYR, SSYR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL

      // 1)      Constants

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0) return;

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Precision' );

      // Do Test 1

      // Copy A & Compute its 1-Norm:

      slaset('Full', N, N, ZERO, ZERO, WORK, N );

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

      // Norm of A - USU'

      for (J = 1; J <= N; J++) { // 20
         ssyr('L', N, -SD( J ), U( 1, J ), 1, WORK, N );
      } // 20

      if ( N > 1 && KBAND == 1 ) {
         for (J = 1; J <= N - 1; J++) { // 30
            ssyr2('L', N, -SE( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N );
         } // 30
      }

      WNORM = SLANSY( '1', 'L', N, WORK, N, WORK( N**2+1 ) );

      if ( ANORM > WNORM ) {
         RESULT[1] = ( WNORM / ANORM ) / ( N*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
         } else {
            RESULT[1] = min( WNORM / ANORM, REAL( N ) ) / ( N*ULP );
         }
      }

      // Do Test 2

      // Compute  UU' - I

      sgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N );

      for (J = 1; J <= N; J++) { // 40
         WORK[( N+1 )*( J-1 )+1] = WORK( ( N+1 )*( J-1 )+1 ) - ONE;
      } // 40

      RESULT[2] = min( double( N ), SLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ) ) / ( N*ULP );

      }
