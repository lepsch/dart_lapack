// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

void dgtsvx(
  final String FACT,
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Array<double> DLF_,
  final Array<double> DF_,
  final Array<double> DUF_,
  final Array<double> DU2_,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Box<double> RCOND,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DLF = DLF_.having();
  final DF = DF_.having();
  final DUF = DUF_.having();
  final DU2 = DU2_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0;

  INFO.value = 0;
  final NOFACT = lsame(FACT, 'N');
  final NOTRAN = lsame(TRANS, 'N');
  if (!NOFACT && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -14;
  } else if (LDX < max(1, N)) {
    INFO.value = -16;
  }
  if (INFO.value != 0) {
    xerbla('DGTSVX', -INFO.value);
    return;
  }

  if (NOFACT) {
    // Compute the LU factorization of A.
    dcopy(N, D, 1, DF, 1);
    if (N > 1) {
      dcopy(N - 1, DL, 1, DLF, 1);
      dcopy(N - 1, DU, 1, DUF, 1);
    }
    dgttrf(N, DLF, DF, DUF, DU2, IPIV, INFO);

    // Return if INFO is non-zero.
    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.
  final NORM = NOTRAN ? '1' : 'I';
  final ANORM = dlangt(NORM, N, DL, D, DU);

  // Compute the reciprocal of the condition number of A.
  dgtcon(NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO);

  // Compute the solution vectors X.
  dlacpy('Full', N, NRHS, B, LDB, X, LDX);
  dgttrs(TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solutions and
  // compute error bounds and backward error estimates for them.
  dgtrfs(TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX,
      FERR, BERR, WORK, IWORK, INFO);

  // Set INFO = N+1 if the matrix is singular to working precision.
  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
