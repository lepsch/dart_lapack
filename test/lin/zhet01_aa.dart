// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/blas/ztrmm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlacgv.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zhet01_aa(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AFAC_,
  final int LDAFAC,
  final Array<int> IPIV_,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final IPIV = IPIV_.having();
  final C = C_.having(ld: LDC);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Determine EPS and the norm of A.

  final EPS = dlamch('Epsilon');
  final ANORM = zlanhe('1', UPLO, N, A, LDA, RWORK);

  // Initialize C to the tridiagonal matrix T.

  zlaset('Full', N, N, Complex.zero, Complex.zero, C, LDC);
  zlacpy('F', 1, N, AFAC(1, 1), LDAFAC + 1, C(1, 1), LDC + 1);
  if (N > 1) {
    if (lsame(UPLO, 'U')) {
      zlacpy('F', 1, N - 1, AFAC(1, 2), LDAFAC + 1, C(1, 2), LDC + 1);
      zlacpy('F', 1, N - 1, AFAC(1, 2), LDAFAC + 1, C(2, 1), LDC + 1);
      zlacgv(N - 1, C(2, 1).asArray(), LDC + 1);
    } else {
      zlacpy('F', 1, N - 1, AFAC(2, 1), LDAFAC + 1, C(1, 2), LDC + 1);
      zlacpy('F', 1, N - 1, AFAC(2, 1), LDAFAC + 1, C(2, 1), LDC + 1);
      zlacgv(N - 1, C(1, 2).asArray(), LDC + 1);
    }

    // Call ZTRMM to form the product U' * D (or L * D ).

    if (lsame(UPLO, 'U')) {
      ztrmm('Left', UPLO, 'Conjugate transpose', 'Unit', N - 1, N, Complex.one,
          AFAC(1, 2), LDAFAC, C(2, 1), LDC);
    } else {
      ztrmm('Left', UPLO, 'No transpose', 'Unit', N - 1, N, Complex.one,
          AFAC(2, 1), LDAFAC, C(2, 1), LDC);
    }

    // Call ZTRMM again to multiply by U (or L ).

    if (lsame(UPLO, 'U')) {
      ztrmm('Right', UPLO, 'No transpose', 'Unit', N, N - 1, Complex.one,
          AFAC(1, 2), LDAFAC, C(1, 2), LDC);
    } else {
      ztrmm('Right', UPLO, 'Conjugate transpose', 'Unit', N, N - 1, Complex.one,
          AFAC(2, 1), LDAFAC, C(1, 2), LDC);
    }

    // Apply hermitian pivots

    for (var J = N; J >= 1; J--) {
      final I = IPIV[J];
      if (I != J) zswap(N, C(J, 1).asArray(), LDC, C(I, 1).asArray(), LDC);
    }
    for (var J = N; J >= 1; J--) {
      final I = IPIV[J];
      if (I != J) zswap(N, C(1, J).asArray(), 1, C(1, I).asArray(), 1);
    }
  }

  // Compute the difference  C - A .

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J; I++) {
        C[I][J] -= A[I][J];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = J; I <= N; I++) {
        C[I][J] -= A[I][J];
      }
    }
  }

  // Compute norm( C - A ) / ( N * norm(A) * EPS )

  RESID.value = zlanhe('1', UPLO, N, C, LDC, RWORK);

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }
}
