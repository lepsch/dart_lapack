// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dasum.dart';
import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/intrinsics/nint.dart';
import 'package:dart_lapack/src/intrinsics/sign.dart';
import 'package:dart_lapack/src/matrix.dart';

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
  final V = V_.having();
  final X = X_.having();
  final ISGN = ISGN_.having();
  const ITMAX = 5;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;

  if (KASE.value == 0) {
    for (var I = 1; I <= N; I++) {
      X[I] = ONE / N;
    }
    KASE.value = 1;
    _JUMP = 1;
    return;
  }

  switch (_JUMP) {
    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
    case 1:
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

    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
    case 2:
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
    // X HAS BEEN OVERWRITTEN BY A*X.
    case 3:
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

    // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
    case 4:
      final JLAST = _J;
      _J = idamax(N, X, 1);
      if ((X[JLAST] != X[_J].abs()) && (_ITER < ITMAX)) {
        _ITER++;
        continue L50;
      }

      continue L120;
    // ITERATION COMPLETE.  FINAL STAGE.
    L120:
    case 120:
      var ALTSGN = ONE;
      for (var I = 1; I <= N; I++) {
        X[I] = ALTSGN * (ONE + (I - 1) / (N - 1));
        ALTSGN = -ALTSGN;
      }
      KASE.value = 1;
      _JUMP = 5;
      return;
    // X HAS BEEN OVERWRITTEN BY A*X.
    case 5:
      double TEMP = TWO * (dasum(N, X, 1) / (3 * N));
      if (TEMP > EST.value) {
        dcopy(N, X, 1, V, 1);
        EST.value = TEMP;
      }
  }
  KASE.value = 0;
}
