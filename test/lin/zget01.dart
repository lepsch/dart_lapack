// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zdotu.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/blas/ztrmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlaswp.dart';

void zget01(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AFAC_,
  final int LDAFAC,
  final Array<int> IPIV_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final IPIV = IPIV_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if M = 0 or N = 0.

  if (M <= 0 || N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Determine EPS and the norm of A.

  final EPS = dlamch('Epsilon');
  final ANORM = zlange('1', M, N, A, LDA, RWORK);

  // Compute the product L*U and overwrite AFAC with the result.
  // A column at a time of the product is obtained, starting with
  // column N.

  for (var K = N; K >= 1; K--) {
    if (K > M) {
      ztrmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC,
          AFAC(1, K).asArray(), 1);
    } else {
      // Compute elements (K+1:M,K)

      final T = AFAC[K][K];
      if (K + 1 <= M) {
        zscal(M - K, T, AFAC(K + 1, K).asArray(), 1);
        zgemv('No transpose', M - K, K - 1, Complex.one, AFAC(K + 1, 1), LDAFAC,
            AFAC(1, K).asArray(), 1, Complex.one, AFAC(K + 1, K).asArray(), 1);
      }

      // Compute the (K,K) element

      AFAC[K][K] = T +
          zdotu(K - 1, AFAC(K, 1).asArray(), LDAFAC, AFAC(1, K).asArray(), 1);

      // Compute elements (1:K-1,K)

      ztrmv('Lower', 'No transpose', 'Unit', K - 1, AFAC, LDAFAC,
          AFAC(1, K).asArray(), 1);
    }
  }
  zlaswp(N, AFAC, LDAFAC, 1, min(M, N), IPIV, -1);

  // Compute the difference  L*U - A  and store in AFAC.

  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= M; I++) {
      AFAC[I][J] -= A[I][J];
    }
  }

  // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

  RESID.value = zlange('1', M, N, AFAC, LDAFAC, RWORK);

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }
}
