      void zlaic1(JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                J, JOB;
      double             SEST, SESTPR;
      Complex         C, GAMMA, S;
      // ..
      // .. Array Arguments ..
      Complex         W( J ), X( J );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double             HALF, FOUR;
      const              HALF = 0.5, FOUR = 4.0 ;
      // ..
      // .. Local Scalars ..
      double             ABSALP, ABSEST, ABSGAM, B, EPS, NORMA, S1, S2, SCL, T, TEST, TMP, ZETA1, ZETA2;
      Complex         ALPHA, COSINE, SINE;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, SQRT
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      //- Complex         ZDOTC;
      // EXTERNAL DLAMCH, ZDOTC
      // ..
      // .. Executable Statements ..

      EPS = dlamch( 'Epsilon' );
      ALPHA = ZDOTC( J, X, 1, W, 1 );

      ABSALP = ( ALPHA ).abs();
      ABSGAM = ( GAMMA ).abs();
      ABSEST = ( SEST ).abs();

      if ( JOB == 1 ) {

         // Estimating largest singular value

         // special cases

         if ( SEST == ZERO ) {
            S1 = max( ABSGAM, ABSALP );
            if ( S1 == ZERO ) {
               S = ZERO;
               C = ONE;
               SESTPR = ZERO;
            } else {
               S = ALPHA / S1;
               C = GAMMA / S1;
               TMP = DBLE( sqrt( S*DCONJG( S )+C*DCONJG( C ) ) );
               S = S / TMP;
               C = C / TMP;
               SESTPR = S1*TMP;
            }
            return;
         } else if ( ABSGAM <= EPS*ABSEST ) {
            S = ONE;
            C = ZERO;
            TMP = max( ABSEST, ABSALP );
            S1 = ABSEST / TMP;
            S2 = ABSALP / TMP;
            SESTPR = TMP*sqrt( S1*S1+S2*S2 );
            return;
         } else if ( ABSALP <= EPS*ABSEST ) {
            S1 = ABSGAM;
            S2 = ABSEST;
            if ( S1 <= S2 ) {
               S = ONE;
               C = ZERO;
               SESTPR = S2;
            } else {
               S = ZERO;
               C = ONE;
               SESTPR = S1;
            }
            return;
         } else if ( ABSEST <= EPS*ABSALP || ABSEST <= EPS*ABSGAM ) {
            S1 = ABSGAM;
            S2 = ABSALP;
            if ( S1 <= S2 ) {
               TMP = S1 / S2;
               SCL = sqrt( ONE+TMP*TMP );
               SESTPR = S2*SCL;
               S = ( ALPHA / S2 ) / SCL;
               C = ( GAMMA / S2 ) / SCL;
            } else {
               TMP = S2 / S1;
               SCL = sqrt( ONE+TMP*TMP );
               SESTPR = S1*SCL;
               S = ( ALPHA / S1 ) / SCL;
               C = ( GAMMA / S1 ) / SCL;
            }
            return;
         } else {

            // normal case

            ZETA1 = ABSALP / ABSEST;
            ZETA2 = ABSGAM / ABSEST;

            B = ( ONE-ZETA1*ZETA1-ZETA2*ZETA2 )*HALF;
            C = ZETA1*ZETA1;
            if ( B > ZERO ) {
               T = DBLE( C / ( B+sqrt( B*B+C ) ) );
            } else {
               T = (sqrt( B*B+C ) - B).toDouble();
            }

            SINE = -( ALPHA / ABSEST ) / T;
            COSINE = -( GAMMA / ABSEST ) / ( ONE+T );
            TMP = DBLE( sqrt( SINE * DCONJG( SINE ) + COSINE * DCONJG( COSINE ) ) );

            S = SINE / TMP;
            C = COSINE / TMP;
            SESTPR = sqrt( T+ONE )*ABSEST;
            return;
         }

      } else if ( JOB == 2 ) {

         // Estimating smallest singular value

         // special cases

         if ( SEST == ZERO ) {
            SESTPR = ZERO;
            if ( max( ABSGAM, ABSALP ) == ZERO ) {
               SINE = ONE;
               COSINE = ZERO;
            } else {
               SINE = -DCONJG( GAMMA );
               COSINE = DCONJG( ALPHA );
            }
            S1 = max( ( SINE ).abs(), ( COSINE ).abs() );
            S = SINE / S1;
            C = COSINE / S1;
            TMP = DBLE( sqrt( S*DCONJG( S )+C*DCONJG( C ) ) );
            S = S / TMP;
            C = C / TMP;
            return;
         } else if ( ABSGAM <= EPS*ABSEST ) {
            S = ZERO;
            C = ONE;
            SESTPR = ABSGAM;
            return;
         } else if ( ABSALP <= EPS*ABSEST ) {
            S1 = ABSGAM;
            S2 = ABSEST;
            if ( S1 <= S2 ) {
               S = ZERO;
               C = ONE;
               SESTPR = S1;
            } else {
               S = ONE;
               C = ZERO;
               SESTPR = S2;
            }
            return;
         } else if ( ABSEST <= EPS*ABSALP || ABSEST <= EPS*ABSGAM ) {
            S1 = ABSGAM;
            S2 = ABSALP;
            if ( S1 <= S2 ) {
               TMP = S1 / S2;
               SCL = sqrt( ONE+TMP*TMP );
               SESTPR = ABSEST*( TMP / SCL );
               S = -( DCONJG( GAMMA ) / S2 ) / SCL;
               C = ( DCONJG( ALPHA ) / S2 ) / SCL;
            } else {
               TMP = S2 / S1;
               SCL = sqrt( ONE+TMP*TMP );
               SESTPR = ABSEST / SCL;
               S = -( DCONJG( GAMMA ) / S1 ) / SCL;
               C = ( DCONJG( ALPHA ) / S1 ) / SCL;
            }
            return;
         } else {

            // normal case

            ZETA1 = ABSALP / ABSEST;
            ZETA2 = ABSGAM / ABSEST;

            NORMA = max( ONE+ZETA1*ZETA1+ZETA1*ZETA2, ZETA1*ZETA2+ZETA2*ZETA2 );

            // See if root is closer to zero or to ONE

            TEST = ONE + TWO*( ZETA1-ZETA2 )*( ZETA1+ZETA2 );
            if ( TEST >= ZERO ) {

               // root is close to zero, compute directly

               B = ( ZETA1*ZETA1+ZETA2*ZETA2+ONE )*HALF;
               C = ZETA2*ZETA2;
               T = DBLE( C / ( B+sqrt( ( B*B-C ).abs() ) ) );
               SINE = ( ALPHA / ABSEST ) / ( ONE-T );
               COSINE = -( GAMMA / ABSEST ) / T;
               SESTPR = sqrt( T+FOUR*EPS*EPS*NORMA )*ABSEST;
            } else {

               // root is closer to ONE, shift by that amount

               B = ( ZETA2*ZETA2+ZETA1*ZETA1-ONE )*HALF;
               C = ZETA1*ZETA1;
               if ( B >= ZERO ) {
                  T = DBLE( -C / ( B+sqrt( B*B+C ) ) );
               } else {
                  T = (B - sqrt( B*B+C )).toDouble();
               }
               SINE = -( ALPHA / ABSEST ) / T;
               COSINE = -( GAMMA / ABSEST ) / ( ONE+T );
               SESTPR = sqrt( ONE+T+FOUR*EPS*EPS*NORMA )*ABSEST;
            }
            TMP = DBLE( sqrt( SINE * DCONJG( SINE ) + COSINE * DCONJG( COSINE ) ) );
            S = SINE / TMP;
            C = COSINE / TMP;
            return;

         }
      }
      return;
      }
