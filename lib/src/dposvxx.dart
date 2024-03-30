import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dla_porpvgrw.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaqsy.dart';
import 'package:lapack/src/dlascl2.dart';
import 'package:lapack/src/dpoequb.dart';
import 'package:lapack/src/dporfsx.dart';
import 'package:lapack/src/dpotrs.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/dpotrf.dart';
import 'package:lapack/src/xerbla.dart';

void dposvxx(
  final String FACT,
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AF_,
  final int LDAF,
  final Box<String> EQUED,
  final Array<double> S_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> X_,
  final int LDX,
  final Box<double> RCOND,
  final Box<double> RPVGRW,
  final Array<double> BERR_,
  final int N_ERR_BNDS,
  final Matrix<double> ERR_BNDS_NORM_,
  final Matrix<double> ERR_BNDS_COMP_,
  final int NPARAMS,
  final Array<double> PARAMS_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final ERR_BNDS_NORM = ERR_BNDS_NORM_.having(ld: NRHS);
  final ERR_BNDS_COMP = ERR_BNDS_COMP_.having(ld: NRHS);
  final S = S_.having();
  final BERR = BERR_.having();
  final PARAMS = PARAMS_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  // const FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3;
  // const RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6;
  // const CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9;
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
  // pivot growth is set here, the rest is initialized in DPORFSX.

  RPVGRW.value = ZERO;

  // Test the input parameters.  PARAMS is not tested until DPORFSX.

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
        SMIN = min(SMIN, S[J]);
        SMAX = max(SMAX, S[J]);
      }
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
    xerbla('DPOSVXX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    dpoequb(N, A, LDA, S, SCOND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      dlaqsy(UPLO, N, A, LDA, S, SCOND.value, AMAX.value, EQUED);
      RCEQU = lsame(EQUED.value, 'Y');
    }
  }

  // Scale the right-hand side.

  if (RCEQU) dlascl2(N, NRHS, S, B, LDB);

  if (NOFACT || EQUIL) {
    // Compute the Cholesky factorization of A.

    dlacpy(UPLO, N, N, A, LDA, AF, LDAF);
    dpotrf(UPLO, N, AF, LDAF, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value != 0) {
      // Pivot in column INFO.value is exactly 0
      // Compute the reciprocal pivot growth factor of the
      // leading rank-deficient INFO.value columns of A.

      RPVGRW.value = dla_porpvgrw(UPLO, INFO.value, A, LDA, AF, LDAF, WORK);
      return;
    }
  }

  // Compute the reciprocal growth factor RPVGRW.value.

  RPVGRW.value = dla_porpvgrw(UPLO, N, A, LDA, AF, LDAF, WORK);

  // Compute the solution matrix X.

  dlacpy('Full', N, NRHS, B, LDB, X, LDX);
  dpotrs(UPLO, N, NRHS, AF, LDAF, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  dporfsx(
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
      IWORK,
      INFO);

  // Scale solutions.

  if (RCEQU) {
    dlascl2(N, NRHS, S, X, LDX);
  }
}
