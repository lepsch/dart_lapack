      SUBROUTINE SLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                J, JOB;
      REAL               C, GAMMA, S, SEST, SESTPR
      // ..
      // .. Array Arguments ..
      REAL               W( J ), X( J )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 ;
      REAL               HALF, FOUR
      const              HALF = 0.5E0, FOUR = 4.0E0 ;
      // ..
      // .. Local Scalars ..
      REAL               ABSALP, ABSEST, ABSGAM, ALPHA, B, COSINE, EPS, NORMA, S1, S2, SINE, T, TEST, TMP, ZETA1, ZETA2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN, SQRT
      // ..
      // .. External Functions ..
      REAL               SDOT, SLAMCH
      // EXTERNAL SDOT, SLAMCH
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' )
      ALPHA = SDOT( J, X, 1, W, 1 )

      ABSALP = ABS( ALPHA )
      ABSGAM = ABS( GAMMA )
      ABSEST = ABS( SEST )

      if ( JOB == 1 ) {

         // Estimating largest singular value

         // special cases

         if ( SEST == ZERO ) {
            S1 = MAX( ABSGAM, ABSALP )
            if ( S1 == ZERO ) {
               S = ZERO
               C = ONE
               SESTPR = ZERO
            } else {
               S = ALPHA / S1
               C = GAMMA / S1
               TMP = SQRT( S*S+C*C )
               S = S / TMP
               C = C / TMP
               SESTPR = S1*TMP
            }
            RETURN
         } else if ( ABSGAM.LE.EPS*ABSEST ) {
            S = ONE
            C = ZERO
            TMP = MAX( ABSEST, ABSALP )
            S1 = ABSEST / TMP
            S2 = ABSALP / TMP
            SESTPR = TMP*SQRT( S1*S1+S2*S2 )
            RETURN
         } else if ( ABSALP.LE.EPS*ABSEST ) {
            S1 = ABSGAM
            S2 = ABSEST
            if ( S1.LE.S2 ) {
               S = ONE
               C = ZERO
               SESTPR = S2
            } else {
               S = ZERO
               C = ONE
               SESTPR = S1
            }
            RETURN
         } else if ( ABSEST.LE.EPS*ABSALP || ABSEST.LE.EPS*ABSGAM ) {
            S1 = ABSGAM
            S2 = ABSALP
            if ( S1.LE.S2 ) {
               TMP = S1 / S2
               S = SQRT( ONE+TMP*TMP )
               SESTPR = S2*S
               C = ( GAMMA / S2 ) / S
               S = SIGN( ONE, ALPHA ) / S
            } else {
               TMP = S2 / S1
               C = SQRT( ONE+TMP*TMP )
               SESTPR = S1*C
               S = ( ALPHA / S1 ) / C
               C = SIGN( ONE, GAMMA ) / C
            }
            RETURN
         } else {

            // normal case

            ZETA1 = ALPHA / ABSEST
            ZETA2 = GAMMA / ABSEST

            B = ( ONE-ZETA1*ZETA1-ZETA2*ZETA2 )*HALF
            C = ZETA1*ZETA1
            if ( B > ZERO ) {
               T = C / ( B+SQRT( B*B+C ) )
            } else {
               T = SQRT( B*B+C ) - B
            }

            SINE = -ZETA1 / T
            COSINE = -ZETA2 / ( ONE+T )
            TMP = SQRT( SINE*SINE+COSINE*COSINE )
            S = SINE / TMP
            C = COSINE / TMP
            SESTPR = SQRT( T+ONE )*ABSEST
            RETURN
         }

      } else if ( JOB == 2 ) {

         // Estimating smallest singular value

         // special cases

         if ( SEST == ZERO ) {
            SESTPR = ZERO
            if ( MAX( ABSGAM, ABSALP ) == ZERO ) {
               SINE = ONE
               COSINE = ZERO
            } else {
               SINE = -GAMMA
               COSINE = ALPHA
            }
            S1 = MAX( ABS( SINE ), ABS( COSINE ) )
            S = SINE / S1
            C = COSINE / S1
            TMP = SQRT( S*S+C*C )
            S = S / TMP
            C = C / TMP
            RETURN
         } else if ( ABSGAM.LE.EPS*ABSEST ) {
            S = ZERO
            C = ONE
            SESTPR = ABSGAM
            RETURN
         } else if ( ABSALP.LE.EPS*ABSEST ) {
            S1 = ABSGAM
            S2 = ABSEST
            if ( S1.LE.S2 ) {
               S = ZERO
               C = ONE
               SESTPR = S1
            } else {
               S = ONE
               C = ZERO
               SESTPR = S2
            }
            RETURN
         } else if ( ABSEST.LE.EPS*ABSALP || ABSEST.LE.EPS*ABSGAM ) {
            S1 = ABSGAM
            S2 = ABSALP
            if ( S1.LE.S2 ) {
               TMP = S1 / S2
               C = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST*( TMP / C )
               S = -( GAMMA / S2 ) / C
               C = SIGN( ONE, ALPHA ) / C
            } else {
               TMP = S2 / S1
               S = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST / S
               C = ( ALPHA / S1 ) / S
               S = -SIGN( ONE, GAMMA ) / S
            }
            RETURN
         } else {

            // normal case

            ZETA1 = ALPHA / ABSEST
            ZETA2 = GAMMA / ABSEST

            NORMA = MAX( ONE+ZETA1*ZETA1+ABS( ZETA1*ZETA2 ), ABS( ZETA1*ZETA2 )+ZETA2*ZETA2 )

            // See if root is closer to zero or to ONE

            TEST = ONE + TWO*( ZETA1-ZETA2 )*( ZETA1+ZETA2 )
            if ( TEST.GE.ZERO ) {

               // root is close to zero, compute directly

               B = ( ZETA1*ZETA1+ZETA2*ZETA2+ONE )*HALF
               C = ZETA2*ZETA2
               T = C / ( B+SQRT( ABS( B*B-C ) ) )
               SINE = ZETA1 / ( ONE-T )
               COSINE = -ZETA2 / T
               SESTPR = SQRT( T+FOUR*EPS*EPS*NORMA )*ABSEST
            } else {

               // root is closer to ONE, shift by that amount

               B = ( ZETA2*ZETA2+ZETA1*ZETA1-ONE )*HALF
               C = ZETA1*ZETA1
               if ( B.GE.ZERO ) {
                  T = -C / ( B+SQRT( B*B+C ) )
               } else {
                  T = B - SQRT( B*B+C )
               }
               SINE = -ZETA1 / T
               COSINE = -ZETA2 / ( ONE+T )
               SESTPR = SQRT( ONE+T+FOUR*EPS*EPS*NORMA )*ABSEST
            }
            TMP = SQRT( SINE*SINE+COSINE*COSINE )
            S = SINE / TMP
            C = COSINE / TMP
            RETURN

         }
      }
      RETURN

      // End of SLAIC1

      }
