// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dspr.dart';
import 'package:lapack/src/blas/dtpmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dppt01(
  final String UPLO,
  final int N,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final AFAC = AFAC_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = dlansp('1', UPLO, N, A, RWORK);
  if (ANORM <= ZERO) {
    RESID.value = ONE / EPS;
    return;
  }

  // Compute the product U'*U, overwriting U.

  if (lsame(UPLO, 'U')) {
    var KC = (N * (N - 1)) ~/ 2 + 1;
    for (var K = N; K >= 1; K--) {
      // Compute the (K,K) element of the result.

      final T = ddot(K, AFAC(KC), 1, AFAC(KC), 1);
      AFAC[KC + K - 1] = T;

      // Compute the rest of column K.

      if (K > 1) {
        dtpmv('Upper', 'Transpose', 'Non-unit', K - 1, AFAC, AFAC(KC), 1);
        KC -= (K - 1);
      }
    }

    // Compute the product L*L', overwriting L.
  } else {
    var KC = (N * (N + 1)) ~/ 2;
    for (var K = N; K >= 1; K--) {
      // Add a multiple of column K of the factor L to each of
      // columns K+1 through N.

      if (K < N) {
        dspr('Lower', N - K, ONE, AFAC(KC + 1), 1, AFAC(KC + N - K + 1));
      }

      // Scale column K by the diagonal element.

      final T = AFAC[KC];
      dscal(N - K + 1, T, AFAC(KC), 1);

      KC -= (N - K + 2);
    }
  }

  // Compute the difference  L*L' - A (or U'*U - A).

  final NPP = N * (N + 1) ~/ 2;
  for (var I = 1; I <= NPP; I++) {
    AFAC[I] -= A[I];
  }

  // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

  RESID.value = dlansp('1', UPLO, N, AFAC, RWORK);

  RESID.value = ((RESID.value / N) / ANORM) / EPS;
}
