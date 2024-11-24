// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void dlacn2(
  final int N,
  final Array<double> V_,
  final Array<double> X_,
  final Array<int> ISGN_,
  final Box<double> EST,
  final Box<int> KASE,
  final Array<int> ISAVE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having();
  final X = X_.having();
  final ISGN = ISGN_.having();
  const ITMAX = 5;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  int I, JLAST;
  double ALTSGN, ESTOLD, TEMP, XS;

  if (KASE.value == 0) {
    for (I = 1; I <= N; I++) {
      X[I] = ONE / N;
    }
    KASE.value = 1;
    ISAVE[1] = 1;
    return;
  }

  void iterationLoop() {
    // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
    for (I = 1; I <= N; I++) {
      X[I] = ZERO;
    }
    X[ISAVE[2]] = ONE;
    KASE.value = 1;
    ISAVE[1] = 3;
  }

  void iterationComplete() {
    // ITERATION COMPLETE.  FINAL STAGE.
    ALTSGN = ONE;
    for (I = 1; I <= N; I++) {
      X[I] = ALTSGN * (ONE + (I - 1) / (N - 1));
      ALTSGN = -ALTSGN;
    }
    KASE.value = 1;
    ISAVE[1] = 5;
  }

  switch (ISAVE[1]) {
    case 1:
      // ................ ENTRY   (ISAVE( 1 ) = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      if (N == 1) {
        V[1] = X[1];
        EST.value = V[1].abs();

        break; // ... QUIT
      }
      EST.value = dasum(N, X, 1);

      for (I = 1; I <= N; I++) {
        if (X[I] >= ZERO) {
          X[I] = ONE;
        } else {
          X[I] = -ONE;
        }
        ISGN[I] = X[I].round();
      }
      KASE.value = 2;
      ISAVE[1] = 2;
      return;

    case 2:
      // ................ ENTRY   (ISAVE( 1 ) = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      ISAVE[2] = idamax(N, X, 1);
      ISAVE[3] = 2;
      iterationLoop();
      return;

    case 3:
      // ................ ENTRY   (ISAVE( 1 ) = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      dcopy(N, X, 1, V, 1);
      ESTOLD = EST.value;
      EST.value = dasum(N, V, 1);
      var hasConverged = true;
      for (I = 1; I <= N; I++) {
        if (X[I] >= ZERO) {
          XS = ONE;
        } else {
          XS = -ONE;
        }
        if (XS.round() != ISGN[I]) {
          hasConverged = false;
          break;
        }
      }
      if ( // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
          hasConverged ||
              // TEST FOR CYCLING.
              EST.value <= ESTOLD) {
        iterationComplete();
        return;
      }

      for (I = 1; I <= N; I++) {
        if (X[I] >= ZERO) {
          X[I] = ONE;
        } else {
          X[I] = -ONE;
        }
        ISGN[I] = X[I].round();
      }
      KASE.value = 2;
      ISAVE[1] = 4;
      return;

    case 4:
      // ................ ENTRY   (ISAVE( 1 ) = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      JLAST = ISAVE[2];
      ISAVE[2] = idamax(N, X, 1);
      if ((X[JLAST] != X[ISAVE[2]].abs()) && (ISAVE[3] < ITMAX)) {
        ISAVE[3]++;
        iterationLoop();
        return;
      }
      iterationComplete();
      return;

    case 5:
      // ................ ENTRY   (ISAVE( 1 ) = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      TEMP = TWO * (dasum(N, X, 1) / (3 * N));
      if (TEMP > EST.value) {
        dcopy(N, X, 1, V, 1);
        EST.value = TEMP;
      }
  }

  KASE.value = 0;
}
