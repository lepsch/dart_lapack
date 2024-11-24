// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/blas/dtrmm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

void dsyt01_aa(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AFAC_,
  final int LDAFAC,
  final Array<int> IPIV_,
  final Matrix<double> C_,
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
  final ANORM = dlansy('1', UPLO, N, A, LDA, RWORK);

  // Initialize C to the tridiagonal matrix T.

  dlaset('Full', N, N, ZERO, ZERO, C, LDC);
  dlacpy('F', 1, N, AFAC(1, 1), LDAFAC + 1, C(1, 1), LDC + 1);
  if (N > 1) {
    if (lsame(UPLO, 'U')) {
      dlacpy('F', 1, N - 1, AFAC(1, 2), LDAFAC + 1, C(1, 2), LDC + 1);
      dlacpy('F', 1, N - 1, AFAC(1, 2), LDAFAC + 1, C(2, 1), LDC + 1);
    } else {
      dlacpy('F', 1, N - 1, AFAC(2, 1), LDAFAC + 1, C(1, 2), LDC + 1);
      dlacpy('F', 1, N - 1, AFAC(2, 1), LDAFAC + 1, C(2, 1), LDC + 1);
    }

    // Call DTRMM to form the product U' * D (or L * D ).

    if (lsame(UPLO, 'U')) {
      dtrmm('Left', UPLO, 'Transpose', 'Unit', N - 1, N, ONE, AFAC(1, 2),
          LDAFAC, C(2, 1), LDC);
    } else {
      dtrmm('Left', UPLO, 'No transpose', 'Unit', N - 1, N, ONE, AFAC(2, 1),
          LDAFAC, C(2, 1), LDC);
    }

    // Call DTRMM again to multiply by U (or L ).

    if (lsame(UPLO, 'U')) {
      dtrmm('Right', UPLO, 'No transpose', 'Unit', N, N - 1, ONE, AFAC(1, 2),
          LDAFAC, C(1, 2), LDC);
    } else {
      dtrmm('Right', UPLO, 'Transpose', 'Unit', N, N - 1, ONE, AFAC(2, 1),
          LDAFAC, C(1, 2), LDC);
    }
  }

  // Apply symmetric pivots

  for (var J = N; J >= 1; J--) {
    var I = IPIV[J];
    if (I != J) dswap(N, C(J, 1).asArray(), LDC, C(I, 1).asArray(), LDC);
  }
  for (var J = N; J >= 1; J--) {
    var I = IPIV[J];
    if (I != J) dswap(N, C(1, J).asArray(), 1, C(1, I).asArray(), 1);
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

  RESID.value = dlansy('1', UPLO, N, C, LDC, RWORK);

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }
}
