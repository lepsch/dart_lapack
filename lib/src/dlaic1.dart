import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlaic1(JOB, J, X, SEST, W, GAMMA, SESTPR, S, final int C) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                J, JOB;
      double             C, GAMMA, S, SEST, SESTPR;
      double             W( J ), X( J );
      // ..

      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double             HALF, FOUR;
      const              HALF = 0.5, FOUR = 4.0 ;
      double             ABSALP, ABSEST, ABSGAM, ALPHA, B, COSINE, EPS, NORMA, S1, S2, SINE, T, TEST, TMP, ZETA1, ZETA2;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN, SQRT
      // ..
      // .. External Functions ..
      //- double             DDOT, DLAMCH;
      // EXTERNAL DDOT, DLAMCH

      EPS = dlamch( 'Epsilon' );
      ALPHA = ddot( J, X, 1, W, 1 );

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
               TMP = sqrt( S*S+C*C );
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
               S = sqrt( ONE+TMP*TMP );
               SESTPR = S2*S;
               C = ( GAMMA / S2 ) / S;
               S = sign( ONE, ALPHA ) / S;
            } else {
               TMP = S2 / S1;
               C = sqrt( ONE+TMP*TMP );
               SESTPR = S1*C;
               S = ( ALPHA / S1 ) / C;
               C = sign( ONE, GAMMA ) / C;
            }
            return;
         } else {

            // normal case

            ZETA1 = ALPHA / ABSEST;
            ZETA2 = GAMMA / ABSEST;

            B = ( ONE-ZETA1*ZETA1-ZETA2*ZETA2 )*HALF;
            C = ZETA1*ZETA1;
            if ( B > ZERO ) {
               T = C / ( B+sqrt( B*B+C ) );
            } else {
               T = sqrt( B*B+C ) - B;
            }

            SINE = -ZETA1 / T;
            COSINE = -ZETA2 / ( ONE+T );
            TMP = sqrt( SINE*SINE+COSINE*COSINE );
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
               SINE = -GAMMA;
               COSINE = ALPHA;
            }
            S1 = max( ( SINE ).abs(), ( COSINE ).abs() );
            S = SINE / S1;
            C = COSINE / S1;
            TMP = sqrt( S*S+C*C );
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
               C = sqrt( ONE+TMP*TMP );
               SESTPR = ABSEST*( TMP / C );
               S = -( GAMMA / S2 ) / C;
               C = sign( ONE, ALPHA ) / C;
            } else {
               TMP = S2 / S1;
               S = sqrt( ONE+TMP*TMP );
               SESTPR = ABSEST / S;
               C = ( ALPHA / S1 ) / S;
               S = -sign( ONE, GAMMA ) / S;
            }
            return;
         } else {

            // normal case

            ZETA1 = ALPHA / ABSEST;
            ZETA2 = GAMMA / ABSEST;

            NORMA = max( ONE+ZETA1*ZETA1+( ZETA1*ZETA2 ).abs(), ( ZETA1*ZETA2 ).abs()+ZETA2*ZETA2 );

            // See if root is closer to zero or to ONE

            TEST = ONE + TWO*( ZETA1-ZETA2 )*( ZETA1+ZETA2 );
            if ( TEST >= ZERO ) {

               // root is close to zero, compute directly

               B = ( ZETA1*ZETA1+ZETA2*ZETA2+ONE )*HALF;
               C = ZETA2*ZETA2;
               T = C / ( B+sqrt( ( B*B-C ).abs() ) );
               SINE = ZETA1 / ( ONE-T );
               COSINE = -ZETA2 / T;
               SESTPR = sqrt( T+FOUR*EPS*EPS*NORMA )*ABSEST;
            } else {

               // root is closer to ONE, shift by that amount

               B = ( ZETA2*ZETA2+ZETA1*ZETA1-ONE )*HALF;
               C = ZETA1*ZETA1;
               if ( B >= ZERO ) {
                  T = -C / ( B+sqrt( B*B+C ) );
               } else {
                  T = B - sqrt( B*B+C );
               }
               SINE = -ZETA1 / T;
               COSINE = -ZETA2 / ( ONE+T );
               SESTPR = sqrt( ONE+T+FOUR*EPS*EPS*NORMA )*ABSEST;
            }
            TMP = sqrt( SINE*SINE+COSINE*COSINE );
            S = SINE / TMP;
            C = COSINE / TMP;
            return;

         }
      }
      }
