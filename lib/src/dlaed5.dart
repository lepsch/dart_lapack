import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlaed5(I, D, Z, DELTA, RHO, DLAM ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I;
      double             DLAM, RHO;
      // ..
      // .. Array Arguments ..
      double             D( 2 ), DELTA( 2 ), Z( 2 );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, FOUR;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, FOUR = 4.0 ;
      // ..
      // .. Local Scalars ..
      double             B, C, DEL, TAU, TEMP, W;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      DEL = D( 2 ) - D( 1 );
      if ( I == 1 ) {
         W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL;
         if ( W > ZERO ) {
            B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) );
            C = RHO*Z( 1 )*Z( 1 )*DEL;

            // B > ZERO, always

            TAU = TWO*C / ( B+sqrt( ( B*B-FOUR*C ).abs() ) );
            DLAM = D( 1 ) + TAU;
            DELTA[1] = -Z( 1 ) / TAU;
            DELTA[2] = Z( 2 ) / ( DEL-TAU );
         } else {
            B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) );
            C = RHO*Z( 2 )*Z( 2 )*DEL;
            if ( B > ZERO ) {
               TAU = -TWO*C / ( B+sqrt( B*B+FOUR*C ) );
            } else {
               TAU = ( B-sqrt( B*B+FOUR*C ) ) / TWO;
            }
            DLAM = D( 2 ) + TAU;
            DELTA[1] = -Z( 1 ) / ( DEL+TAU );
            DELTA[2] = -Z( 2 ) / TAU;
         }
         TEMP = sqrt( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) );
         DELTA[1] = DELTA( 1 ) / TEMP;
         DELTA[2] = DELTA( 2 ) / TEMP;
      } else {

      // Now I=2

         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) );
         C = RHO*Z( 2 )*Z( 2 )*DEL;
         if ( B > ZERO ) {
            TAU = ( B+sqrt( B*B+FOUR*C ) ) / TWO;
         } else {
            TAU = TWO*C / ( -B+sqrt( B*B+FOUR*C ) );
         }
         DLAM = D( 2 ) + TAU;
         DELTA[1] = -Z( 1 ) / ( DEL+TAU );
         DELTA[2] = -Z( 2 ) / TAU;
         TEMP = sqrt( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) );
         DELTA[1] = DELTA( 1 ) / TEMP;
         DELTA[2] = DELTA( 2 ) / TEMP;
      }
      return;
      }
