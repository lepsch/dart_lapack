      void dstt21(N, KBAND, AD, AE, SD, SE, U, LDU, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KBAND, LDU, N;
      // ..
      // .. Array Arguments ..
      double             AD( * ), AE( * ), RESULT( 2 ), SD( * ), SE( * ), U( LDU, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLASET, DSYR, DSYR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // 1)      Constants

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0) return;

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Precision' );

      // Do Test 1

      // Copy A & Compute its 1-Norm:

      dlaset('Full', N, N, ZERO, ZERO, WORK, N );

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
         dsyr('L', N, -SD( J ), U( 1, J ), 1, WORK, N );
      } // 20

      if ( N > 1 && KBAND == 1 ) {
         for (J = 1; J <= N - 1; J++) { // 30
            dsyr2('L', N, -SE( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N );
         } // 30
      }

      WNORM = DLANSY( '1', 'L', N, WORK, N, WORK( N**2+1 ) );

      if ( ANORM > WNORM ) {
         RESULT[1] = ( WNORM / ANORM ) / ( N*ULP );
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ( min( WNORM, N*ANORM ) / ANORM ) / ( N*ULP );
         } else {
            RESULT[1] = min( WNORM / ANORM, DBLE( N ) ) / ( N*ULP );
         }
      }

      // Do Test 2

      // Compute  UU' - I

      dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N );

      for (J = 1; J <= N; J++) { // 40
         WORK[( N+1 )*( J-1 )+1] = WORK( ( N+1 )*( J-1 )+1 ) - ONE;
      } // 40

      RESULT[2] = min( DBLE( N ), DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ) ) / ( N*ULP );

      return;
      }
