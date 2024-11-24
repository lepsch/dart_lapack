// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/nint.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';

double _ESTOLD = 0;
int _ITER = 0, _J = 0, _JUMP = 0;

void dlacon(
  final int N,
  final Array<double> V_,
  final Array<double> X_,
  final Array<int> ISGN_,
  final Box<double> EST,
  final Box<int> KASE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having();
  final X = X_.having();
  final ISGN = ISGN_.having();
  const ITMAX = 5;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;

  // .. Save statement ..
  // SAVE;

  if (KASE.value == 0) {
    for (var I = 1; I <= N; I++) {
      X[I] = ONE / N;
    }
    KASE.value = 1;
    _JUMP = 1;
    return;
  }

  switch (_JUMP) {
    case 1:
      // ................ ENTRY   (_JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      if (N == 1) {
        V[1] = X[1];
        EST.value = V[1].abs();
        KASE.value = 0;
        return;
      }
      EST.value = dasum(N, X, 1);

      for (var I = 1; I <= N; I++) {
        X[I] = sign(ONE, X[I]);
        ISGN[I] = nint(X[I]);
      }
      KASE.value = 2;
      _JUMP = 2;
      return;

    case 2:
      // ................ ENTRY   (_JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      _J = idamax(N, X, 1);
      _ITER = 2;
      continue L50;
    L50:
    case 50:
      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
      for (var I = 1; I <= N; I++) {
        X[I] = ZERO;
      }
      X[_J] = ONE;
      KASE.value = 1;
      _JUMP = 3;
      return;
    case 3:
      // ................ ENTRY   (_JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      dcopy(N, X, 1, V, 1);
      _ESTOLD = EST.value;
      EST.value = dasum(N, V, 1);
      var isConverged = true;
      for (var I = 1; I <= N; I++) {
        if (nint(sign(ONE, X[I])) != ISGN[I]) {
          isConverged = false;
          break;
        }
      }
      // REPEATED sign VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      if (isConverged) continue L120;

      // TEST FOR CYCLING.
      if (EST.value <= _ESTOLD) continue L120;

      for (var I = 1; I <= N; I++) {
        X[I] = sign(ONE, X[I]);
        ISGN[I] = nint(X[I]);
      }
      KASE.value = 2;
      _JUMP = 4;
      return;

    case 4:
      // ................ ENTRY   (_JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      final JLAST = _J;
      _J = idamax(N, X, 1);
      if ((X[JLAST] != X[_J].abs()) && (_ITER < ITMAX)) {
        _ITER++;
        continue L50;
      }

      continue L120;
    L120:
    case 120:
      // ITERATION COMPLETE.  FINAL STAGE.

      var ALTSGN = ONE;
      for (var I = 1; I <= N; I++) {
        X[I] = ALTSGN * (ONE + (I - 1) / (N - 1));
        ALTSGN = -ALTSGN;
      }
      KASE.value = 1;
      _JUMP = 5;
      return;
    case 5:
      // ................ ENTRY   (_JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      double TEMP = TWO * (dasum(N, X, 1) / (3 * N));
      if (TEMP > EST.value) {
        dcopy(N, X, 1, V, 1);
        EST.value = TEMP;
      }
  }
  KASE.value = 0;
}
