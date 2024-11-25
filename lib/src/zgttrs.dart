// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgtts2.dart';

void zgttrs(
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Array<Complex> DU2_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  bool NOTRAN;
  int ITRANS, J, JB, NB;

  INFO.value = 0;
  NOTRAN = (TRANS == 'N' || TRANS == 'n');
  if (!NOTRAN &&
      !(TRANS == 'T' || TRANS == 't') &&
      !(TRANS == 'C' || TRANS == 'c')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(N, 1)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZGTTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  // Decode TRANS

  if (NOTRAN) {
    ITRANS = 0;
  } else if (TRANS == 'T' || TRANS == 't') {
    ITRANS = 1;
  } else {
    ITRANS = 2;
  }

  // Determine the number of right-hand sides to solve at a time.

  if (NRHS == 1) {
    NB = 1;
  } else {
    NB = max(1, ilaenv(1, 'ZGTTRS', TRANS, N, NRHS, -1, -1));
  }

  if (NB >= NRHS) {
    zgtts2(ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
  } else {
    for (J = 1; J <= NRHS; J += NB) {
      JB = min(NRHS - J + 1, NB);
      zgtts2(ITRANS, N, JB, DL, D, DU, DU2, IPIV, B(1, J), LDB);
    }
  }
}
