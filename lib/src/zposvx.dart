import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/cholesky/top/zpotrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlaqhe.dart';
import 'package:lapack/src/zpocon.dart';
import 'package:lapack/src/zpoequ.dart';
import 'package:lapack/src/zporfs.dart';
import 'package:lapack/src/zpotrs.dart';

void zposvx(
  final String FACT,
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
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
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final S = S_.having();
  final FERR = FERR_.having();
  final BERR = BERR_.having();
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool EQUIL, NOFACT, RCEQU;
  int I, J;
  double ANORM, BIGNUM = 0, SMAX, SMIN, SMLNUM = 0;
  final INFEQU = Box(0);
  final AMAX = Box(0.0), SCOND = Box(0.0);

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  EQUIL = lsame(FACT, 'E');
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
  } else if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDAF < max(1, N)) {
    INFO.value = -8;
  } else if (lsame(FACT, 'F') && !(RCEQU || lsame(EQUED.value, 'N'))) {
    INFO.value = -9;
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
        INFO.value = -10;
      } else if (N > 0) {
        SCOND.value = max(SMIN, SMLNUM) / min(SMAX, BIGNUM);
      } else {
        SCOND.value = ONE;
      }
    }
    if (INFO.value == 0) {
      if (LDB < max(1, N)) {
        INFO.value = -12;
      } else if (LDX < max(1, N)) {
        INFO.value = -14;
      }
    }
  }

  if (INFO.value != 0) {
    xerbla('ZPOSVX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    zpoequ(N, A, LDA, S, SCOND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      zlaqhe(UPLO, N, A, LDA, S, SCOND.value, AMAX.value, EQUED);
      RCEQU = lsame(EQUED.value, 'Y');
    }
  }

  // Scale the right hand side.

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

    zlacpy(UPLO, N, N, A, LDA, AF, LDAF);
    zpotrf(UPLO, N, AF, LDAF, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      RCOND.value = ZERO;
      return;
    }
  }

  // Compute the norm of the matrix A.

  ANORM = zlanhe('1', UPLO, N, A, LDA, RWORK);

  // Compute the reciprocal of the condition number of A.

  zpocon(UPLO, N, AF, LDAF, ANORM, RCOND, WORK, RWORK, INFO);

  // Compute the solution matrix X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zpotrs(UPLO, N, NRHS, AF, LDAF, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  zporfs(UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, FERR, BERR, WORK,
      RWORK, INFO);

  // Transform the solution matrix X to a solution of the original
  // system.

  if (RCEQU) {
    for (J = 1; J <= NRHS; J++) {
      // 50
      for (I = 1; I <= N; I++) {
        // 40
        X[I][J] = S[I].toComplex() * X[I][J];
      } // 40
    } // 50
    for (J = 1; J <= NRHS; J++) {
      // 60
      FERR[J] = FERR[J] / SCOND.value;
    } // 60
  }

  // Set INFO.value = N+1 if the matrix is singular to working precision.

  if (RCOND.value < dlamch('Epsilon')) INFO.value = N + 1;
}
