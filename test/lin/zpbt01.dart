// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zher.dart';
import 'package:dart_lapack/src/blas/ztrmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlanhb.dart';

void zpbt01(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AFAC_,
  final int LDAFAC,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  // Quick exit if N = 0.

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlanhb('1', UPLO, N, KD, A, LDA, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Check the imaginary parts of the diagonal elements and return with
  // an error code if any are nonzero.

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      if (AFAC[KD + 1][J].imaginary != ZERO) {
        RESID.value = ONE / EPS;
        return;
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      if (AFAC[1][J].imaginary != ZERO) {
        RESID.value = ONE / EPS;
        return;
      }
    }
  }

  // Compute the product U'*U, overwriting U.

  if (lsame(UPLO, 'U')) {
    for (var K = N; K >= 1; K--) {
      final KC = max(1, KD + 2 - K);
      final KLEN = KD + 1 - KC;

      // Compute the (K,K) element of the result.

      final AKK =
          zdotc(KLEN + 1, AFAC(KC, K).asArray(), 1, AFAC(KC, K).asArray(), 1)
              .real;
      AFAC[KD + 1][K] = AKK.toComplex();

      // Compute the rest of column K.

      if (KLEN > 0) {
        ztrmv('Upper', 'Conjugate', 'Non-unit', KLEN, AFAC(KD + 1, K - KLEN),
            LDAFAC - 1, AFAC(KC, K).asArray(), 1);
      }
    }

    // UPLO = 'L':  Compute the product L*L', overwriting L.
  } else {
    for (var K = N; K >= 1; K--) {
      final KLEN = min(KD, N - K);

      // Add a multiple of column K of the factor L to each of
      // columns K+1 through N.

      if (KLEN > 0) {
        zher('Lower', KLEN, ONE, AFAC(2, K).asArray(), 1, AFAC(1, K + 1),
            LDAFAC - 1);
      }

      // Scale column K by the diagonal element.

      final AKK = AFAC[1][K].real;
      zdscal(KLEN + 1, AKK, AFAC(1, K).asArray(), 1);
    }
  }

  // Compute the difference  L*L' - A  or  U'*U - A.

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      final MU = max(1, KD + 2 - J);
      for (var I = MU; I <= KD + 1; I++) {
        AFAC[I][J] -= A[I][J];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      final ML = min(KD + 1, N - J + 1);
      for (var I = 1; I <= ML; I++) {
        AFAC[I][J] -= A[I][J];
      }
    }
  }

  // Compute norm( L*L' - A ) / ( N * norm(A) * EPS )

  RESID.value = zlanhb('1', UPLO, N, KD, AFAC, LDAFAC, RWORK);

  RESID.value = ((RESID.value / N) / ANORM) / EPS;
}
