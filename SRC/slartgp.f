      SUBROUTINE SLARTGP( F, G, CS, SN, R )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               CS, F, G, R, SN
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               ONE
      const              ONE = 1.0E0 ;
      REAL               TWO
      const              TWO = 2.0E0 ;
      // ..
      // .. Local Scalars ..
      // LOGICAL            FIRST
      int                COUNT, I;
      REAL               EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
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
         SAFMIN = SLAMCH( 'S' )
         EPS = SLAMCH( 'E' )
         SAFMN2 = SLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( SLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
         // FIRST = false;
      // END IF
      if ( G.EQ.ZERO ) {
         CS = SIGN( ONE, F )
         SN = ZERO
         R = ABS( F )
      } else if ( F.EQ.ZERO ) {
         CS = ZERO
         SN = SIGN( ONE, G )
         R = ABS( G )
      } else {
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         if ( SCALE.GE.SAFMX2 ) {
            COUNT = 0
            } // 10
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            if (SCALE.GE.SAFMX2 .AND. COUNT .LT. 20) GO TO 10;
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            for (I = 1; I <= COUNT; I++) { // 20
               R = R*SAFMX2
            } // 20
         } else if ( SCALE.LE.SAFMN2 ) {
            COUNT = 0
            } // 30
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            if (SCALE.LE.SAFMN2) GO TO 30;
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            for (I = 1; I <= COUNT; I++) { // 40
               R = R*SAFMN2
            } // 40
         } else {
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         }
         if ( R.LT.ZERO ) {
            CS = -CS
            SN = -SN
            R = -R
         }
      }
      RETURN

      // End of SLARTGP

      }
