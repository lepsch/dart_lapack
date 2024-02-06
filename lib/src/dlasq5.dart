import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlasq5(I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, IEEE, EPS ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               IEEE;
      int                I0, N0, PP;
      double             DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, SIGMA, EPS;
      double             Z( * );
      // ..

// =====================================================================

      // .. Parameter ..
      double             ZERO, HALF;
      const              ZERO = 0.0, HALF = 0.5 ;
      int                J4, J4P2;
      double             D, EMIN, TEMP, DTHRESH;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      if( ( N0-I0-1 ) <= 0 ) return;

      DTHRESH = EPS*(SIGMA+TAU);
      if (TAU < DTHRESH*HALF) TAU = ZERO;
      if ( TAU != ZERO ) {
      J4 = 4*I0 + PP - 3;
      EMIN = Z( J4+4 );
      D = Z( J4 ) - TAU;
      DMIN = D;
      DMIN1 = -Z( J4 );

      if ( IEEE ) {

         // Code for IEEE arithmetic.

         if ( PP == 0 ) {
            for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 10
               Z[J4-2] = D + Z( J4-1 );
               TEMP = Z( J4+1 ) / Z( J4-2 );
               D = D*TEMP - TAU;
               DMIN = min( DMIN, D );
               Z[J4] = Z( J4-1 )*TEMP;
               EMIN = min( Z( J4 ), EMIN );
            } // 10
         } else {
            for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 20
               Z[J4-3] = D + Z( J4 );
               TEMP = Z( J4+2 ) / Z( J4-3 );
               D = D*TEMP - TAU;
               DMIN = min( DMIN, D );
               Z[J4-1] = Z( J4 )*TEMP;
               EMIN = min( Z( J4-1 ), EMIN );
            } // 20
         }

         // Unroll last two steps.

         DNM2 = D;
         DMIN2 = DMIN;
         J4 = 4*( N0-2 ) - PP;
         J4P2 = J4 + 2*PP - 1;
         Z[J4-2] = DNM2 + Z( J4P2 );
         Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU;
         DMIN = min( DMIN, DNM1 );

         DMIN1 = DMIN;
         J4 = J4 + 4;
         J4P2 = J4 + 2*PP - 1;
         Z[J4-2] = DNM1 + Z( J4P2 );
         Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU;
         DMIN = min( DMIN, DN );

      } else {

         // Code for non IEEE arithmetic.

         if ( PP == 0 ) {
            for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 30
               Z[J4-2] = D + Z( J4-1 );
               if ( D < ZERO ) {
                  return;
               } else {
                  Z[J4] = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) );
                  D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU;
               }
               DMIN = min( DMIN, D );
               EMIN = min( EMIN, Z( J4 ) );
            } // 30
         } else {
            for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 40
               Z[J4-3] = D + Z( J4 );
               if ( D < ZERO ) {
                  return;
               } else {
                  Z[J4-1] = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) );
                  D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU;
               }
               DMIN = min( DMIN, D );
               EMIN = min( EMIN, Z( J4-1 ) );
            } // 40
         }

         // Unroll last two steps.

         DNM2 = D;
         DMIN2 = DMIN;
         J4 = 4*( N0-2 ) - PP;
         J4P2 = J4 + 2*PP - 1;
         Z[J4-2] = DNM2 + Z( J4P2 );
         if ( DNM2 < ZERO ) {
            return;
         } else {
            Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU;
         }
         DMIN = min( DMIN, DNM1 );

         DMIN1 = DMIN;
         J4 = J4 + 4;
         J4P2 = J4 + 2*PP - 1;
         Z[J4-2] = DNM1 + Z( J4P2 );
         if ( DNM1 < ZERO ) {
            return;
         } else {
            Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU;
         }
         DMIN = min( DMIN, DN );

      }
      } else {
      // This is the version that sets d's to zero if they are small enough
         J4 = 4*I0 + PP - 3;
         EMIN = Z( J4+4 );
         D = Z( J4 ) - TAU;
         DMIN = D;
         DMIN1 = -Z( J4 );
         if ( IEEE ) {

      // Code for IEEE arithmetic.

            if ( PP == 0 ) {
               for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 50
                  Z[J4-2] = D + Z( J4-1 );
                  TEMP = Z( J4+1 ) / Z( J4-2 );
                  D = D*TEMP - TAU;
                  if (D < DTHRESH) D = ZERO;
                  DMIN = min( DMIN, D );
                  Z[J4] = Z( J4-1 )*TEMP;
                  EMIN = min( Z( J4 ), EMIN );
               } // 50
            } else {
               for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 60
                  Z[J4-3] = D + Z( J4 );
                  TEMP = Z( J4+2 ) / Z( J4-3 );
                  D = D*TEMP - TAU;
                  if (D < DTHRESH) D = ZERO;
                  DMIN = min( DMIN, D );
                  Z[J4-1] = Z( J4 )*TEMP;
                  EMIN = min( Z( J4-1 ), EMIN );
               } // 60
            }

      // Unroll last two steps.

            DNM2 = D;
            DMIN2 = DMIN;
            J4 = 4*( N0-2 ) - PP;
            J4P2 = J4 + 2*PP - 1;
            Z[J4-2] = DNM2 + Z( J4P2 );
            Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU;
            DMIN = min( DMIN, DNM1 );

            DMIN1 = DMIN;
            J4 = J4 + 4;
            J4P2 = J4 + 2*PP - 1;
            Z[J4-2] = DNM1 + Z( J4P2 );
            Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU;
            DMIN = min( DMIN, DN );

         } else {

      // Code for non IEEE arithmetic.

            if ( PP == 0 ) {
               for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 70
                  Z[J4-2] = D + Z( J4-1 );
                  if ( D < ZERO ) {
                     return;
                  } else {
                     Z[J4] = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) );
                     D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU;
                  }
                  if (D < DTHRESH) D = ZERO;
                  DMIN = min( DMIN, D );
                  EMIN = min( EMIN, Z( J4 ) );
               } // 70
            } else {
               for (J4 = 4*I0; J4 <= 4*( N0-3 ); J4 += 4) { // 80
                  Z[J4-3] = D + Z( J4 );
                  if ( D < ZERO ) {
                     return;
                  } else {
                     Z[J4-1] = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) );
                     D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU;
                  }
                  if (D < DTHRESH) D = ZERO;
                  DMIN = min( DMIN, D );
                  EMIN = min( EMIN, Z( J4-1 ) );
               } // 80
            }

      // Unroll last two steps.

            DNM2 = D;
            DMIN2 = DMIN;
            J4 = 4*( N0-2 ) - PP;
            J4P2 = J4 + 2*PP - 1;
            Z[J4-2] = DNM2 + Z( J4P2 );
            if ( DNM2 < ZERO ) {
               return;
            } else {
               Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
               DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU;
            }
            DMIN = min( DMIN, DNM1 );

            DMIN1 = DMIN;
            J4 = J4 + 4;
            J4P2 = J4 + 2*PP - 1;
            Z[J4-2] = DNM1 + Z( J4P2 );
            if ( DNM1 < ZERO ) {
               return;
            } else {
               Z[J4] = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) );
               DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU;
            }
            DMIN = min( DMIN, DN );

         }
      }

      Z[J4+2] = DN;
      Z[4*N0-PP] = EMIN;
      return;
      }
