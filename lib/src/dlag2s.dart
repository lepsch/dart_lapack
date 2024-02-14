import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlag2s(final int M, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> SA_, final int LDSA, final Box<int> INFO,) {
  final A = A_.dim();
  final SA = SA_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDSA, M, N;
      double               SA( LDSA, * );
      double             A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      RMAX = SLAMCH( 'O' );
      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            if ( ( A( I, J ) < -RMAX ) || ( A( I, J ) > RMAX ) ) {
               INFO = 1;
               GO TO 30;
            }
            SA[I][J] = double( A( I, J ) );
         } // 10
      } // 20
      INFO = 0;
      } // 30
      }
