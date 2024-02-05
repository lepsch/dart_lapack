import 'package:lapack/src/complex.dart';

import 'common.dart';

      bool zslect(Z ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex         Z;
      // ..

// =====================================================================

      // .. Parameters ..
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             RMIN, X;
      // ..
      // .. Scalars in Common ..
      // int                sslct.SELDIM, sslct.SELOPT;
      // // ..
      // // .. Arrays in Common ..
      // bool               sslct.SELVAL( 20 );
      // double             sslct.SELWI( 20 ), sslct.SELWR( 20 );
      // ..
      // .. Common blocks ..
      // COMMON / SSLCT / sslct.SELOPT, sslct.SELDIM, sslct.SELVAL, sslct.SELWR, sslct.SELWI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      if ( sslct.SELOPT == 0 ) {
         ZSLECT = ( Z.toDouble() < ZERO );
      } else {
         RMIN = ABS( Z-DCMPLX( sslct.SELWR( 1 ), sslct.SELWI( 1 ) ) );
         ZSLECT = sslct.SELVAL( 1 );
         for (I = 2; I <= sslct.SELDIM; I++) { // 10
            X = ABS( Z-DCMPLX( sslct.SELWR( I ), sslct.SELWI( I ) ) );
            if ( X <= RMIN ) {
               RMIN = X;
               ZSLECT = sslct.SELVAL( I );
            }
         } // 10
      }
      return;
      }
