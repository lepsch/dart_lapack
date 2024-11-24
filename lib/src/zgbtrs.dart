// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zgeru.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/blas/ztbsv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';

void zgbtrs(
  final String TRANS,
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  bool LNOTI, NOTRAN;
  int I, J, KD, L, LM;

  // Test the input parameters.

  INFO.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0) {
    INFO.value = -3;
  } else if (KU < 0) {
    INFO.value = -4;
  } else if (NRHS < 0) {
    INFO.value = -5;
  } else if (LDAB < (2 * KL + KU + 1)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZGBTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  KD = KU + KL + 1;
  LNOTI = KL > 0;

  if (NOTRAN) {
    // Solve  A*X = B.

    // Solve L*X = B, overwriting B with X.

    // L is represented as a product of permutations and unit lower
    // triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
    // where each transformation L(i) is a rank-one modification of
    // the identity matrix.

    if (LNOTI) {
      for (J = 1; J <= N - 1; J++) {
        LM = min(KL, N - J);
        L = IPIV[J];
        if (L != J) zswap(NRHS, B(L, 1).asArray(), LDB, B(J, 1).asArray(), LDB);
        zgeru(LM, NRHS, -Complex.one, AB(KD + 1, J).asArray(), 1,
            B(J, 1).asArray(), LDB, B(J + 1, 1), LDB);
      }
    }

    for (I = 1; I <= NRHS; I++) {
      // Solve U*X = B, overwriting B with X.

      ztbsv('Upper', 'No transpose', 'Non-unit', N, KL + KU, AB, LDAB,
          B(1, I).asArray(), 1);
    }
  } else if (lsame(TRANS, 'T')) {
    // Solve A**T * X = B.

    for (I = 1; I <= NRHS; I++) {
      // Solve U**T * X = B, overwriting B with X.

      ztbsv('Upper', 'Transpose', 'Non-unit', N, KL + KU, AB, LDAB,
          B(1, I).asArray(), 1);
    }

    // Solve L**T * X = B, overwriting B with X.

    if (LNOTI) {
      for (J = N - 1; J >= 1; J--) {
        LM = min(KL, N - J);
        zgemv('Transpose', LM, NRHS, -Complex.one, B(J + 1, 1), LDB,
            AB(KD + 1, J).asArray(), 1, Complex.one, B(J, 1).asArray(), LDB);
        L = IPIV[J];
        if (L != J) zswap(NRHS, B(L, 1).asArray(), LDB, B(J, 1).asArray(), LDB);
      }
    }
  } else {
    // Solve A**H * X = B.

    for (I = 1; I <= NRHS; I++) {
      // Solve U**H * X = B, overwriting B with X.

      ztbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KL + KU, AB, LDAB,
          B(1, I).asArray(), 1);
    }

    // Solve L**H * X = B, overwriting B with X.

    if (LNOTI) {
      for (J = N - 1; J >= 1; J--) {
        LM = min(KL, N - J);
        zlacgv(NRHS, B(J, 1).asArray(), LDB);
        zgemv('Conjugate transpose', LM, NRHS, -Complex.one, B(J + 1, 1), LDB,
            AB(KD + 1, J).asArray(), 1, Complex.one, B(J, 1).asArray(), LDB);
        zlacgv(NRHS, B(J, 1).asArray(), LDB);
        L = IPIV[J];
        if (L != J) zswap(NRHS, B(L, 1).asArray(), LDB, B(J, 1).asArray(), LDB);
      }
    }
  }
}
