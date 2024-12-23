// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/ztpmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlantp.dart';

void ztpt01(
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<Complex> AP_,
  final Array<Complex> AINVP_,
  final Box<double> RCOND,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final AP = AP_.having();
  final AINVP = AINVP_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if N = 0.

  if (N <= 0) {
    RCOND.value = ONE;
    RESID.value = ZERO;
    return;
  }

  // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

  final EPS = dlamch('Epsilon');
  final ANORM = zlantp('1', UPLO, DIAG, N, AP, RWORK);
  final AINVNM = zlantp('1', UPLO, DIAG, N, AINVP, RWORK);
  if (ANORM <= ZERO || AINVNM <= ZERO) {
    RCOND.value = ZERO;
    RESID.value = ONE / EPS;
    return;
  }
  RCOND.value = (ONE / ANORM) / AINVNM;

  // Compute A * AINV, overwriting AINV.

  final UNITD = lsame(DIAG, 'U');
  if (lsame(UPLO, 'U')) {
    var JC = 1;
    for (var J = 1; J <= N; J++) {
      if (UNITD) AINVP[JC + J - 1] = Complex.one;

      // Form the j-th column of A*AINV.

      ztpmv('Upper', 'No transpose', DIAG, J, AP, AINVP(JC), 1);

      // Subtract 1 from the diagonal to form A*AINV - I.

      AINVP[JC + J - 1] -= Complex.one;
      JC += J;
    }
  } else {
    var JC = 1;
    for (var J = 1; J <= N; J++) {
      if (UNITD) AINVP[JC] = Complex.one;

      // Form the j-th column of A*AINV.

      ztpmv('Lower', 'No transpose', DIAG, N - J + 1, AP(JC), AINVP(JC), 1);

      // Subtract 1 from the diagonal to form A*AINV - I.

      AINVP[JC] -= Complex.one;
      JC += N - J + 1;
    }
  }

  // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

  RESID.value = zlantp('1', UPLO, 'Non-unit', N, AINVP, RWORK);

  RESID.value = ((RESID.value * RCOND.value) / N) / EPS;
}
