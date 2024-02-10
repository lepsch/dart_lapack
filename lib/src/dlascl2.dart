import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlascl2(final int M, final int N, final int D, final int X, final int LDX) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                M, N, LDX;
      double             D( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;

      for (J = 1; J <= N; J++) {
         for (I = 1; I <= M; I++) {
            X[I][J] = X( I, J ) * D( I );
         }
      }

      }
