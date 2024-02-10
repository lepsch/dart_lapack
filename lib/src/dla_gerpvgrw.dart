import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      double dla_gerpvgrw(final int N, final int NCOLS, final Matrix<double> A, final int LDA, final int AF, final int LDAF) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N, NCOLS, LDA, LDAF;
      double             A( LDA, * ), AF( LDAF, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             AMAX, UMAX, RPVGRW;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

      RPVGRW = 1.0;

      for (J = 1; J <= NCOLS; J++) {
         AMAX = 0.0;
         UMAX = 0.0;
         for (I = 1; I <= N; I++) {
            AMAX = max( ( A( I, J ) ).abs(), AMAX );
         }
         for (I = 1; I <= J; I++) {
            UMAX = max( ( AF( I, J ) ).abs(), UMAX );
         }
         if ( UMAX /= 0.0 ) {
            RPVGRW = min( AMAX / UMAX, RPVGRW );
         }
      }
      DLA_GERPVGRW = RPVGRW;
      }
