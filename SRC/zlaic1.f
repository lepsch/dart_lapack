      SUBROUTINE ZLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                J, JOB;
      double             SEST, SESTPR;
      COMPLEX*16         C, GAMMA, S
      // ..
      // .. Array Arguments ..
      COMPLEX*16         W( J ), X( J )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 ;
      double             HALF, FOUR;
      const              HALF = 0.5D0, FOUR = 4.0D0 ;
      // ..
      // .. Local Scalars ..
      double             ABSALP, ABSEST, ABSGAM, B, EPS, NORMA, S1, S2, SCL, T, TEST, TMP, ZETA1, ZETA2;
      COMPLEX*16         ALPHA, COSINE, SINE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, SQRT
      // ..
      // .. External Functions ..
      double             DLAMCH;
      COMPLEX*16         ZDOTC
      // EXTERNAL DLAMCH, ZDOTC
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' )
      ALPHA = ZDOTC( J, X, 1, W, 1 )

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
               TMP = DBLE( SQRT( S*DCONJG( S )+C*DCONJG( C ) ) )
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
               SCL = SQRT( ONE+TMP*TMP )
               SESTPR = S2*SCL
               S = ( ALPHA / S2 ) / SCL
               C = ( GAMMA / S2 ) / SCL
            } else {
               TMP = S2 / S1
               SCL = SQRT( ONE+TMP*TMP )
               SESTPR = S1*SCL
               S = ( ALPHA / S1 ) / SCL
               C = ( GAMMA / S1 ) / SCL
            }
            RETURN
         } else {

            // normal case

            ZETA1 = ABSALP / ABSEST
            ZETA2 = ABSGAM / ABSEST

            B = ( ONE-ZETA1*ZETA1-ZETA2*ZETA2 )*HALF
            C = ZETA1*ZETA1
            if ( B > ZERO ) {
               T = DBLE( C / ( B+SQRT( B*B+C ) ) )
            } else {
               T = DBLE( SQRT( B*B+C ) - B )
            }

            SINE = -( ALPHA / ABSEST ) / T
            COSINE = -( GAMMA / ABSEST ) / ( ONE+T )
            TMP = DBLE( SQRT( SINE * DCONJG( SINE ) + COSINE * DCONJG( COSINE ) ) )

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
               SINE = -DCONJG( GAMMA )
               COSINE = DCONJG( ALPHA )
            }
            S1 = MAX( ABS( SINE ), ABS( COSINE ) )
            S = SINE / S1
            C = COSINE / S1
            TMP = DBLE( SQRT( S*DCONJG( S )+C*DCONJG( C ) ) )
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
               SCL = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST*( TMP / SCL )
               S = -( DCONJG( GAMMA ) / S2 ) / SCL
               C = ( DCONJG( ALPHA ) / S2 ) / SCL
            } else {
               TMP = S2 / S1
               SCL = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST / SCL
               S = -( DCONJG( GAMMA ) / S1 ) / SCL
               C = ( DCONJG( ALPHA ) / S1 ) / SCL
            }
            RETURN
         } else {

            // normal case

            ZETA1 = ABSALP / ABSEST
            ZETA2 = ABSGAM / ABSEST

            NORMA = MAX( ONE+ZETA1*ZETA1+ZETA1*ZETA2, ZETA1*ZETA2+ZETA2*ZETA2 )

            // See if root is closer to zero or to ONE

            TEST = ONE + TWO*( ZETA1-ZETA2 )*( ZETA1+ZETA2 )
            if ( TEST.GE.ZERO ) {

               // root is close to zero, compute directly

               B = ( ZETA1*ZETA1+ZETA2*ZETA2+ONE )*HALF
               C = ZETA2*ZETA2
               T = DBLE( C / ( B+SQRT( ABS( B*B-C ) ) ) )
               SINE = ( ALPHA / ABSEST ) / ( ONE-T )
               COSINE = -( GAMMA / ABSEST ) / T
               SESTPR = SQRT( T+FOUR*EPS*EPS*NORMA )*ABSEST
            } else {

               // root is closer to ONE, shift by that amount

               B = ( ZETA2*ZETA2+ZETA1*ZETA1-ONE )*HALF
               C = ZETA1*ZETA1
               if ( B.GE.ZERO ) {
                  T = DBLE( -C / ( B+SQRT( B*B+C ) ) )
               } else {
                  T = DBLE( B - SQRT( B*B+C ) )
               }
               SINE = -( ALPHA / ABSEST ) / T
               COSINE = -( GAMMA / ABSEST ) / ( ONE+T )
               SESTPR = SQRT( ONE+T+FOUR*EPS*EPS*NORMA )*ABSEST
            }
            TMP = DBLE( SQRT( SINE * DCONJG( SINE ) + COSINE * DCONJG( COSINE ) ) )
            S = SINE / TMP
            C = COSINE / TMP
            RETURN

         }
      }
      RETURN

      // End of ZLAIC1

      }
