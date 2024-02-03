      SUBROUTINE DLARTGP( F, G, CS, SN, R );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             CS, F, G, R, SN;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      double             TWO;
      const              TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      // LOGICAL            FIRST
      int                COUNT, I;
      double             EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, SIGN, SQRT
      // ..
      // .. Save statement ..
      // SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
      // ..
      // .. Data statements ..
      // DATA               FIRST / true /
      // ..
      // .. Executable Statements ..

      // IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' );
         EPS = DLAMCH( 'E' );
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( DLAMCH( 'B' ) ) / TWO );
         SAFMX2 = ONE / SAFMN2;
         // FIRST = false;
      // END IF
      if ( G == ZERO ) {
         CS = SIGN( ONE, F );
         SN = ZERO;
         R = ABS( F );
      } else if ( F == ZERO ) {
         CS = ZERO;
         SN = SIGN( ONE, G );
         R = ABS( G );
      } else {
         F1 = F;
         G1 = G;
         SCALE = max( ABS( F1 ), ABS( G1 ) );
         if ( SCALE >= SAFMX2 ) {
            COUNT = 0;
            } // 10
            COUNT = COUNT + 1;
            F1 = F1*SAFMN2;
            G1 = G1*SAFMN2;
            SCALE = max( ABS( F1 ), ABS( G1 ) );
            if (SCALE >= SAFMX2 && COUNT < 20) GO TO 10;
            R = sqrt( F1**2+G1**2 );
            CS = F1 / R;
            SN = G1 / R;
            for (I = 1; I <= COUNT; I++) { // 20
               R = R*SAFMX2;
            } // 20
         } else if ( SCALE <= SAFMN2 ) {
            COUNT = 0;
            } // 30
            COUNT = COUNT + 1;
            F1 = F1*SAFMX2;
            G1 = G1*SAFMX2;
            SCALE = max( ABS( F1 ), ABS( G1 ) );
            if (SCALE <= SAFMN2) GO TO 30;
            R = sqrt( F1**2+G1**2 );
            CS = F1 / R;
            SN = G1 / R;
            for (I = 1; I <= COUNT; I++) { // 40
               R = R*SAFMN2;
            } // 40
         } else {
            R = sqrt( F1**2+G1**2 );
            CS = F1 / R;
            SN = G1 / R;
         }
         if ( R < ZERO ) {
            CS = -CS;
            SN = -SN;
            R = -R;
         }
      }
      return;

      // End of DLARTGP

      }
