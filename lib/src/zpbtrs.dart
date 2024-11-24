// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/ztbsv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zpbtrs(
  final String UPLO,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final B = B_.having(ld: LDB);
  bool UPPER;
  int J;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDAB < KD + 1) {
    INFO.value = -6;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZPBTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B where A = U**H *U.

    for (J = 1; J <= NRHS; J++) {
      // Solve U**H *X = B, overwriting B with X.

      ztbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);

      // Solve U*X = B, overwriting B with X.

      ztbsv('Upper', 'No transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);
    }
  } else {
    // Solve A*X = B where A = L*L**H.

    for (J = 1; J <= NRHS; J++) {
      // Solve L*X = B, overwriting B with X.

      ztbsv('Lower', 'No transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);

      // Solve L**H *X = B, overwriting B with X.

      ztbsv('Lower', 'Conjugate transpose', 'Non-unit', N, KD, AB, LDAB,
          B(1, J).asArray(), 1);
    }
  }
}
