import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbcon.dart';
import 'package:lapack/src/dgbequ.dart';
import 'package:lapack/src/dgbrfs.dart';
import 'package:lapack/src/dgbtrf.dart';
import 'package:lapack/src/dgbtrs.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlangb.dart';
import 'package:lapack/src/dlantb.dart';
import 'package:lapack/src/dlaqgb.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgbsvx(
  final String FACT,
  final String TRANS,
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<double> AB_,
  final int LDAB,
  final Matrix<double> AFB_,
  final int LDAFB,
  final Array<int> IPIV_,
  final Box<String> EQUED,
  final Array<double> R_,
  final Array<double> C_,
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
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAFB);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final R = R_.having();
  final C = C_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
  String NORM;
  int I, J, J1, J2;
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
  } else if (KL < 0) {
    INFO.value = -4;
  } else if (KU < 0) {
    INFO.value = -5;
  } else if (NRHS < 0) {
    INFO.value = -6;
  } else if (LDAB < KL + KU + 1) {
    INFO.value = -8;
  } else if (LDAFB < 2 * KL + KU + 1) {
    INFO.value = -10;
  } else if (lsame(FACT, 'F') &&
      !(ROWEQU || COLEQU || lsame(EQUED.value, 'N'))) {
    INFO.value = -12;
  } else {
    if (ROWEQU) {
      RCMIN = BIGNUM;
      RCMAX = ZERO;
      for (J = 1; J <= N; J++) {
        RCMIN = min(RCMIN, R[J]);
        RCMAX = max(RCMAX, R[J]);
      }
      if (RCMIN <= ZERO) {
        INFO.value = -13;
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
        INFO.value = -14;
      } else if (N > 0) {
        COLCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
      } else {
        COLCND.value = ONE;
      }
    }
    if (INFO.value == 0) {
      if (LDB < max(1, N)) {
        INFO.value = -16;
      } else if (LDX < max(1, N)) {
        INFO.value = -18;
      }
    }
  }

  if (INFO.value != 0) {
    xerbla('DGBSVX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    dgbequ(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      dlaqgb(N, N, KL, KU, AB, LDAB, R, C, ROWCND.value, COLCND.value,
          AMAX.value, EQUED);
      ROWEQU = lsame(EQUED.value, 'R') || lsame(EQUED.value, 'B');
      COLEQU = lsame(EQUED.value, 'C') || lsame(EQUED.value, 'B');
    }
  }

  // Scale the right hand side.

  if (NOTRAN) {
    if (ROWEQU) {
      for (J = 1; J <= NRHS; J++) {
        for (I = 1; I <= N; I++) {
          B[I][J] = R[I] * B[I][J];
        }
      }
    }
  } else if (COLEQU) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        B[I][J] = C[I] * B[I][J];
      }
    }
  }

  if (NOFACT || EQUIL) {
    // Compute the LU factorization of the band matrix A.

    for (J = 1; J <= N; J++) {
      J1 = max(J - KU, 1);
      J2 = min(J + KL, N);
      dcopy(J2 - J1 + 1, AB(KU + 1 - J + J1, J).asArray(), 1,
          AFB(KL + KU + 1 - J + J1, J).asArray(), 1);
    }

    dgbtrf(N, N, KL, KU, AFB, LDAFB, IPIV, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      // Compute the reciprocal pivot growth factor of the
      // leading rank-deficient INFO.value columns of A.

      ANORM = ZERO;
      for (J = 1; J <= INFO.value; J++) {
        for (I = max(KU + 2 - J, 1);
            I <= min(N + KU + 1 - J, KL + KU + 1);
            I++) {
          ANORM = max(ANORM, (AB[I][J]).abs());
        }
      }
      RPVGRW = dlantb('M', 'U', 'N', INFO.value, min(INFO.value - 1, KL + KU),
          AFB(max(1, KL + KU + 2 - INFO.value), 1), LDAFB, WORK);
      if (RPVGRW == ZERO) {
        RPVGRW = ONE;
      } else {
        RPVGRW = ANORM / RPVGRW;
      }
      WORK[1] = RPVGRW;
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
  ANORM = dlangb(NORM, N, KL, KU, AB, LDAB, WORK);
  RPVGRW = dlantb('M', 'U', 'N', N, KL + KU, AFB, LDAFB, WORK);
  if (RPVGRW == ZERO) {
    RPVGRW = ONE;
  } else {
    RPVGRW = dlangb('M', N, KL, KU, AB, LDAB, WORK) / RPVGRW;
  }

  // Compute the reciprocal of the condition number of A.

  dgbcon(NORM, N, KL, KU, AFB, LDAFB, IPIV, ANORM, RCOND, WORK, IWORK, INFO);

  // Compute the solution matrix X.

  dlacpy('Full', N, NRHS, B, LDB, X, LDX);
  dgbtrs(TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  dgbrfs(TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX,
      FERR, BERR, WORK, IWORK, INFO);

  // Transform the solution matrix X to a solution of the original
  // system.

  if (NOTRAN) {
    if (COLEQU) {
      for (J = 1; J <= NRHS; J++) {
        for (I = 1; I <= N; I++) {
          X[I][J] = C[I] * X[I][J];
        }
      }
      for (J = 1; J <= NRHS; J++) {
        FERR[J] /= COLCND.value;
      }
    }
  } else if (ROWEQU) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        X[I][J] = R[I] * X[I][J];
      }
    }
    for (J = 1; J <= NRHS; J++) {
      FERR[J] /= ROWCND.value;
    }
  }

  // Set INFO.value = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;

  WORK[1] = RPVGRW;
}
