// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zher.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';

void zpbtf2(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int J, KLD, KN;
  double AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (LDAB < KD + 1) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZPBTF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  KLD = max(1, LDAB - 1);

  if (UPPER) {
    // Compute the Cholesky factorization A = U**H * U.

    for (J = 1; J <= N; J++) {
      // Compute U(J,J) and test for non-positive-definiteness.

      AJJ = AB[KD + 1][J].real;
      if (AJJ <= ZERO) {
        AB[KD + 1][J] = AJJ.toComplex();
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[KD + 1][J] = AJJ.toComplex();

      // Compute elements J+1:J+KN of row J and update the
      // trailing submatrix within the band.

      KN = min(KD, N - J);
      if (KN > 0) {
        zdscal(KN, ONE / AJJ, AB(KD, J + 1).asArray(), KLD);
        zlacgv(KN, AB(KD, J + 1).asArray(), KLD);
        zher('Upper', KN, -ONE, AB(KD, J + 1).asArray(), KLD, AB(KD + 1, J + 1),
            KLD);
        zlacgv(KN, AB(KD, J + 1).asArray(), KLD);
      }
    }
  } else {
    // Compute the Cholesky factorization A = L*L**H.

    for (J = 1; J <= N; J++) {
      // Compute L(J,J) and test for non-positive-definiteness.

      AJJ = AB[1][J].real;
      if (AJJ <= ZERO) {
        AB[1][J] = AJJ.toComplex();
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[1][J] = AJJ.toComplex();

      // Compute elements J+1:J+KN of column J and update the
      // trailing submatrix within the band.

      KN = min(KD, N - J);
      if (KN > 0) {
        zdscal(KN, ONE / AJJ, AB(2, J).asArray(), 1);
        zher('Lower', KN, -ONE, AB(2, J).asArray(), 1, AB(1, J + 1), KLD);
      }
    }
  }
}
