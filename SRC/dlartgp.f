      SUBROUTINE DLARTGP( F, G, CS, SN, R )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             CS, F, G, R, SN;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      double             ONE;
      const              ONE = 1.0D0 ;
      double             TWO;
      const              TWO = 2.0D0 ;
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
      // DATA               FIRST / .TRUE. /
      // ..
      // .. Executable Statements ..

      // IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
         // FIRST = .FALSE.
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
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 .AND. COUNT .LT. 20 ) GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         } else if ( SCALE.LE.SAFMN2 ) {
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 ) GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
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

      // End of DLARTGP

      }
