import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dla_wwaddw(N, X, Y, final int W) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N;
      double             X( * ), Y( * ), W( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      double             S;
      int                I;

      for (I = 1; I <= N; I++) { // 10
        S = X(I) + W(I);
        S = (S + S) - S;
        Y[I] = ((X(I) - S) + W(I)) + Y(I);
        X[I] = S;
      } // 10
      }
