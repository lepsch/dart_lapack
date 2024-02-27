import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/cholesky/top/zpotrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zla_porpvgrw.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaqhe.dart';
import 'package:lapack/src/zlascl2.dart';
import 'package:lapack/src/zpoequb.dart';
import 'package:lapack/src/zporfsx.dart';
import 'package:lapack/src/zpotrs.dart';

void zposvxx(
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
  final Box<double> RPVGRW,
  final Array<double> BERR_,
  final int N_ERR_BNDS,
  final Matrix<double> ERR_BNDS_NORM_,
  final Matrix<double> ERR_BNDS_COMP_,
  final int NPARAMS,
  final Array<double> PARAMS_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDAF);
  final B = B_.dim(LDB);
  final X = X_.dim(LDX);
  final ERR_BNDS_NORM = ERR_BNDS_NORM_.dim(NRHS);
  final ERR_BNDS_COMP = ERR_BNDS_COMP_.dim(NRHS);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final S = S_.dim();
  final BERR = BERR_.dim();
  final PARAMS = PARAMS_.dim();

  // .. Parameters ..
  const ZERO = 0.0, ONE = 1.0;
  bool EQUIL, NOFACT, RCEQU;
  int J;
  double BIGNUM, SMIN, SMAX, SMLNUM;
  final INFEQU = Box(0);
  final AMAX = Box(0.0), SCOND = Box(0.0);

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  EQUIL = lsame(FACT, 'E');
  SMLNUM = dlamch('Safe minimum');
  BIGNUM = ONE / SMLNUM;
  if (NOFACT || EQUIL) {
    EQUED.value = 'N';
    RCEQU = false;
  } else {
    RCEQU = lsame(EQUED.value, 'Y');
  }

  // Default is failure.  If an input parameter is wrong or
  // factorization fails, make everything look horrible.  Only the
  // pivot growth is set here, the rest is initialized in ZPORFSX.

  RPVGRW.value = ZERO;

  // Test the input parameters.  PARAMS is not tested until ZPORFSX.

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
    xerbla('ZPOSVXX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    zpoequb(N, A, LDA, S, SCOND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      zlaqhe(UPLO, N, A, LDA, S, SCOND.value, AMAX.value, EQUED);
      RCEQU = lsame(EQUED.value, 'Y');
    }
  }

  // Scale the right-hand side.

  if (RCEQU) zlascl2(N, NRHS, S, B, LDB);

  if (NOFACT || EQUIL) {
    // Compute the Cholesky factorization of A.

    zlacpy(UPLO, N, N, A, LDA, AF, LDAF);
    zpotrf(UPLO, N, AF, LDAF, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      // Pivot in column INFO is exactly 0
      // Compute the reciprocal pivot growth factor of the
      // leading rank-deficient INFO columns of A.

      RPVGRW.value = zla_porpvgrw(UPLO, N, A, LDA, AF, LDAF, RWORK);
      return;
    }
  }

  // Compute the reciprocal pivot growth factor RPVGRW.

  RPVGRW.value = zla_porpvgrw(UPLO, N, A, LDA, AF, LDAF, RWORK);

  // Compute the solution matrix X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zpotrs(UPLO, N, NRHS, AF, LDAF, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  zporfsx(
      UPLO,
      EQUED.value,
      N,
      NRHS,
      A,
      LDA,
      AF,
      LDAF,
      S,
      B,
      LDB,
      X,
      LDX,
      RCOND,
      BERR,
      N_ERR_BNDS,
      ERR_BNDS_NORM,
      ERR_BNDS_COMP,
      NPARAMS,
      PARAMS,
      WORK,
      RWORK,
      INFO);

  // Scale solutions.

  if (RCEQU) {
    zlascl2(N, NRHS, S, X, LDX);
  }
}
