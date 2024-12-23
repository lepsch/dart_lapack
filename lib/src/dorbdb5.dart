// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlassq.dart';
import 'package:dart_lapack/src/dorbdb6.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dorbdb5(
  final int M1,
  final int M2,
  final int N,
  final Array<double> X1_,
  final int INCX1,
  final Array<double> X2_,
  final int INCX2,
  final Matrix<double> Q1_,
  final int LDQ1,
  final Matrix<double> Q2_,
  final int LDQ2,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final X1 = X1_.having();
  final X2 = X2_.having();
  final Q1 = Q1_.having(ld: LDQ1);
  final Q2 = Q2_.having(ld: LDQ2);
  final WORK = WORK_.having();
  const REALZERO = 0.0;
  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  double EPS, NORM;
  final CHILDINFO = Box(0);
  final SCL = Box(0.0), SSQ = Box(0.0);

  // Test input arguments

  INFO.value = 0;
  if (M1 < 0) {
    INFO.value = -1;
  } else if (M2 < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (INCX1 < 1) {
    INFO.value = -5;
  } else if (INCX2 < 1) {
    INFO.value = -7;
  } else if (LDQ1 < max(1, M1)) {
    INFO.value = -9;
  } else if (LDQ2 < max(1, M2)) {
    INFO.value = -11;
  } else if (LWORK < N) {
    INFO.value = -13;
  }

  if (INFO.value != 0) {
    xerbla('DORBDB5', -INFO.value);
    return;
  }

  EPS = dlamch('Precision');

  // Project X onto the orthogonal complement of Q if X is nonzero

  SCL.value = REALZERO;
  SSQ.value = REALZERO;
  dlassq(M1, X1, INCX1, SCL, SSQ);
  dlassq(M2, X2, INCX2, SCL, SSQ);
  NORM = SCL.value * sqrt(SSQ.value);

  if (NORM > N * EPS) {
    // Scale vector to unit norm to avoid problems in the caller code.
    // Computing the reciprocal is undesirable but
    //  * xLASCL cannot be used because of the vector increments and
    //  * the round-off error has a negligible impact on
    //    orthogonalization.
    dscal(M1, ONE / NORM, X1, INCX1);
    dscal(M2, ONE / NORM, X2, INCX2);
    dorbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK,
        CHILDINFO);

    // If the projection is nonzero, then return;

    if (dnrm2(M1, X1, INCX1) != REALZERO || dnrm2(M2, X2, INCX2) != REALZERO) {
      return;
    }
  }

  // Project each standard basis vector e_1,...,e_M1 in turn, stopping
  // when a nonzero projection is found

  for (I = 1; I <= M1; I++) {
    for (J = 1; J <= M1; J++) {
      X1[J] = ZERO;
    }
    X1[I] = ONE;
    for (J = 1; J <= M2; J++) {
      X2[J] = ZERO;
    }
    dorbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK,
        CHILDINFO);
    if (dnrm2(M1, X1, INCX1) != REALZERO || dnrm2(M2, X2, INCX2) != REALZERO) {
      return;
    }
  }

  // Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
  // stopping when a nonzero projection is found

  for (I = 1; I <= M2; I++) {
    for (J = 1; J <= M1; J++) {
      X1[J] = ZERO;
    }
    for (J = 1; J <= M2; J++) {
      X2[J] = ZERO;
    }
    X2[I] = ONE;
    dorbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK,
        CHILDINFO);
    if (dnrm2(M1, X1, INCX1) != REALZERO || dnrm2(M2, X2, INCX2) != REALZERO) {
      return;
    }
  }
}
