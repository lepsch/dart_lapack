import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlascl2(M, N, D, X, LDX ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                M, N, LDX;
      // ..
      // .. Array Arguments ..
      double             D( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. Executable Statements ..

      for (J = 1; J <= N; J++) {
         for (I = 1; I <= M; I++) {
            X[I][J] = X( I, J ) * D( I );
         }
      }

      return;
      }
