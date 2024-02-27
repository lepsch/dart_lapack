import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlanhb.dart';
import 'package:lapack/src/zlaqhb.dart';
import 'package:lapack/src/zpbcon.dart';
import 'package:lapack/src/zpbequ.dart';
import 'package:lapack/src/zpbrfs.dart';
import 'package:lapack/src/zpbtrf.dart';
import 'package:lapack/src/zpbtrs.dart';

void zpbsvx(
  final String FACT,
  final String UPLO,
  final int N,
  final int KD,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> AFB_,
  final int LDAFB,
  final Box<String> EQUED,
  final Array<double> S_,
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
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.dim(LDAB);
  final AFB = AFB_.dim(LDAFB);
  final B = B_.dim(LDB);
  final X = X_.dim(LDX);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final S = S_.dim();
  final FERR = FERR_.dim();
  final BERR = BERR_.dim();
  const ZERO = 0.0, ONE = 1.0;
  bool EQUIL, NOFACT, RCEQU, UPPER;
  int I, J, J1, J2;
  double BIGNUM = 0, SMAX, SMIN, SMLNUM = 0;
  final INFEQU = Box(0);
  final AMAX = Box(0.0), ANORM = Box(0.0), SCOND = Box(0.0);

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  EQUIL = lsame(FACT, 'E');
  UPPER = lsame(UPLO, 'U');
  if (NOFACT || EQUIL) {
    EQUED.value = 'N';
    RCEQU = false;
  } else {
    RCEQU = lsame(EQUED.value, 'Y');
    SMLNUM = dlamch('Safe minimum');
    BIGNUM = ONE / SMLNUM;
  }

  // Test the input parameters.

  if (!NOFACT && !EQUIL && !lsame(FACT, 'F')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (KD < 0) {
    INFO.value = -4;
  } else if (NRHS < 0) {
    INFO.value = -5;
  } else if (LDAB < KD + 1) {
    INFO.value = -7;
  } else if (LDAFB < KD + 1) {
    INFO.value = -9;
  } else if (lsame(FACT, 'F') && !(RCEQU || lsame(EQUED.value, 'N'))) {
    INFO.value = -10;
  } else {
    if (RCEQU) {
      SMIN = BIGNUM;
      SMAX = ZERO;
      for (J = 1; J <= N; J++) {
        // 10
        SMIN = min(SMIN, S[J]);
        SMAX = max(SMAX, S[J]);
      } // 10
      if (SMIN <= ZERO) {
        INFO.value = -11;
      } else if (N > 0) {
        SCOND.value = max(SMIN, SMLNUM) / min(SMAX, BIGNUM);
      } else {
        SCOND.value = ONE;
      }
    }
    if (INFO.value == 0) {
      if (LDB < max(1, N)) {
        INFO.value = -13;
      } else if (LDX < max(1, N)) {
        INFO.value = -15;
      }
    }
  }

  if (INFO.value != 0) {
    xerbla('ZPBSVX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    zpbequ(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      zlaqhb(UPLO, N, KD, AB, LDAB, S, SCOND.value, AMAX.value, EQUED);
      RCEQU = lsame(EQUED.value, 'Y');
    }
  }

  // Scale the right-hand side.

  if (RCEQU) {
    for (J = 1; J <= NRHS; J++) {
      // 30
      for (I = 1; I <= N; I++) {
        // 20
        B[I][J] = S[I].toComplex() * B[I][J];
      } // 20
    } // 30
  }

  if (NOFACT || EQUIL) {
    // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

    if (UPPER) {
      for (J = 1; J <= N; J++) {
        // 40
        J1 = max(J - KD, 1);
        zcopy(J - J1 + 1, AB(KD + 1 - J + J1, J).asArray(), 1,
            AFB(KD + 1 - J + J1, J).asArray(), 1);
      } // 40
    } else {
      for (J = 1; J <= N; J++) {
        // 50
        J2 = min(J + KD, N);
        zcopy(J2 - J + 1, AB(1, J).asArray(), 1, AFB(1, J).asArray(), 1);
      } // 50
    }

    zpbtrf(UPLO, N, KD, AFB, LDAFB, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM.value = zlanhb('1', UPLO, N, KD, AB, LDAB, RWORK);

  // Compute the reciprocal of the condition number of A.

  zpbcon(UPLO, N, KD, AFB, LDAFB, ANORM, RCOND, WORK, RWORK, INFO);

  // Compute the solution matrix X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zpbtrs(UPLO, N, KD, NRHS, AFB, LDAFB, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  zpbrfs(UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR,
      WORK, RWORK, INFO);

  // Transform the solution matrix X to a solution of the original
  // system.

  if (RCEQU) {
    for (J = 1; J <= NRHS; J++) {
      // 70
      for (I = 1; I <= N; I++) {
        // 60
        X[I][J] = S[I].toComplex() * X[I][J];
      } // 60
    } // 70
    for (J = 1; J <= NRHS; J++) {
      // 80
      FERR[J] = FERR[J] / SCOND.value;
    } // 80
  }

  // Set INFO.value = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
