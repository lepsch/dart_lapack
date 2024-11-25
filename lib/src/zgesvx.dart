// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zgetrf.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgecon.dart';
import 'package:dart_lapack/src/zgeequ.dart';
import 'package:dart_lapack/src/zgerfs.dart';
import 'package:dart_lapack/src/zgetrs.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlantr.dart';
import 'package:dart_lapack/src/zlaqge.dart';

void zgesvx(
  final String FACT,
  final String TRANS,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
  final Array<int> IPIV_,
  final Box<String> EQUED,
  final Array<double> R_,
  final Array<double> C_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Box<double> RCOND,
  final Array<double> FERR_,
  final Array<double> BERR_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final R = R_.having();
  final C = C_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
  String NORM;
  int I, J;
  double ANORM, BIGNUM = 0, RCMAX, RCMIN, RPVGRW, SMLNUM = 0;
  final INFEQU = Box(0);
  final AMAX = Box(0.0), COLCND = Box(0.0), ROWCND = Box(0.0);

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  EQUIL = lsame(FACT, 'E');
  NOTRAN = lsame(TRANS, 'N');
  if (NOFACT || EQUIL) {
    EQUED.value = 'N';
    ROWEQU = false;
    COLEQU = false;
  } else {
    ROWEQU = lsame(EQUED.value, 'R') || lsame(EQUED.value, 'B');
    COLEQU = lsame(EQUED.value, 'C') || lsame(EQUED.value, 'B');
    SMLNUM = dlamch('Safe minimum');
    BIGNUM = ONE / SMLNUM;
  }

  // Test the input parameters.

  if (!NOFACT && !EQUIL && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDAF < max(1, N)) {
    INFO.value = -8;
  } else if (lsame(FACT, 'F') &&
      !(ROWEQU || COLEQU || lsame(EQUED.value, 'N'))) {
    INFO.value = -10;
  } else {
    if (ROWEQU) {
      RCMIN = BIGNUM;
      RCMAX = ZERO;
      for (J = 1; J <= N; J++) {
        RCMIN = min(RCMIN, R[J]);
        RCMAX = max(RCMAX, R[J]);
      }
      if (RCMIN <= ZERO) {
        INFO.value = -11;
      } else if (N > 0) {
        ROWCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
      } else {
        ROWCND.value = ONE;
      }
    }
    if (COLEQU && INFO.value == 0) {
      RCMIN = BIGNUM;
      RCMAX = ZERO;
      for (J = 1; J <= N; J++) {
        RCMIN = min(RCMIN, C[J]);
        RCMAX = max(RCMAX, C[J]);
      }
      if (RCMIN <= ZERO) {
        INFO.value = -12;
      } else if (N > 0) {
        COLCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
      } else {
        COLCND.value = ONE;
      }
    }
    if (INFO.value == 0) {
      if (LDB < max(1, N)) {
        INFO.value = -14;
      } else if (LDX < max(1, N)) {
        INFO.value = -16;
      }
    }
  }

  if (INFO.value != 0) {
    xerbla('ZGESVX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    zgeequ(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      zlaqge(N, N, A, LDA, R, C, ROWCND.value, COLCND.value, AMAX.value, EQUED);
      ROWEQU = lsame(EQUED.value, 'R') || lsame(EQUED.value, 'B');
      COLEQU = lsame(EQUED.value, 'C') || lsame(EQUED.value, 'B');
    }
  }

  // Scale the right hand side.

  if (NOTRAN) {
    if (ROWEQU) {
      for (J = 1; J <= NRHS; J++) {
        for (I = 1; I <= N; I++) {
          B[I][J] = R[I].toComplex() * B[I][J];
        }
      }
    }
  } else if (COLEQU) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        B[I][J] = C[I].toComplex() * B[I][J];
      }
    }
  }

  if (NOFACT || EQUIL) {
    // Compute the LU factorization of A.

    zlacpy('Full', N, N, A, LDA, AF, LDAF);
    zgetrf(N, N, AF, LDAF, IPIV, INFO);

    // Return if INFO is non-zero.

    if (INFO.value > 0) {
      // Compute the reciprocal pivot growth factor of the
      // leading rank-deficient INFO columns of A.

      RPVGRW = zlantr('M', 'U', 'N', INFO.value, INFO.value, AF, LDAF, RWORK);
      if (RPVGRW == ZERO) {
        RPVGRW = ONE;
      } else {
        RPVGRW = zlange('M', N, INFO.value, A, LDA, RWORK) / RPVGRW;
      }
      RWORK[1] = RPVGRW;
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A and the
  // reciprocal pivot growth factor RPVGRW.

  if (NOTRAN) {
    NORM = '1';
  } else {
    NORM = 'I';
  }
  ANORM = zlange(NORM, N, N, A, LDA, RWORK);
  RPVGRW = zlantr('M', 'U', 'N', N, N, AF, LDAF, RWORK);
  if (RPVGRW == ZERO) {
    RPVGRW = ONE;
  } else {
    RPVGRW = zlange('M', N, N, A, LDA, RWORK) / RPVGRW;
  }

  // Compute the reciprocal of the condition number of A.

  zgecon(NORM, N, AF, LDAF, ANORM, RCOND, WORK, RWORK, INFO);

  // Compute the solution matrix X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zgetrs(TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  zgerfs(TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR,
      WORK, RWORK, INFO);

  // Transform the solution matrix X to a solution of the original
  // system.

  if (NOTRAN) {
    if (COLEQU) {
      for (J = 1; J <= NRHS; J++) {
        for (I = 1; I <= N; I++) {
          X[I][J] = C[I].toComplex() * X[I][J];
        }
      }
      for (J = 1; J <= NRHS; J++) {
        FERR[J] /= COLCND.value;
      }
    }
  } else if (ROWEQU) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        X[I][J] = R[I].toComplex() * X[I][J];
      }
    }
    for (J = 1; J <= NRHS; J++) {
      FERR[J] /= ROWCND.value;
    }
  }

  // Set INFO = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;

  RWORK[1] = RPVGRW;
}
