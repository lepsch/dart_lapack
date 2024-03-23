import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgetrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeequb.dart';
import 'package:lapack/src/zgerfsx.dart';
import 'package:lapack/src/zgetrs.dart';
import 'package:lapack/src/zla_gerpvgrw.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaqge.dart';
import 'package:lapack/src/zlascl2.dart';

void zgesvxx(
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
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final R = R_.having();
  final C = C_.having();
  final BERR = BERR_.having();
  final ERR_BNDS_NORM = ERR_BNDS_NORM_.having(ld: NRHS);
  final ERR_BNDS_COMP = ERR_BNDS_COMP_.having(ld: NRHS);
  final PARAMS = PARAMS_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const ZERO = 0.0, ONE = 1.0;
  bool COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
  int J;
  double BIGNUM, RCMAX, RCMIN, SMLNUM;
  final INFEQU = Box(0);
  final AMAX = Box(0.0), COLCND = Box(0.0), ROWCND = Box(0.0);

  INFO.value = 0;
  NOFACT = lsame(FACT, 'N');
  EQUIL = lsame(FACT, 'E');
  NOTRAN = lsame(TRANS, 'N');
  SMLNUM = dlamch('Safe minimum');
  BIGNUM = ONE / SMLNUM;
  if (NOFACT || EQUIL) {
    EQUED.value = 'N';
    ROWEQU = false;
    COLEQU = false;
  } else {
    ROWEQU = lsame(EQUED.value, 'R') || lsame(EQUED.value, 'B');
    COLEQU = lsame(EQUED.value, 'C') || lsame(EQUED.value, 'B');
  }

  // Default is failure.  If an input parameter is wrong or
  // factorization fails, make everything look horrible.  Only the
  // pivot growth is set here, the rest is initialized in ZGERFSX.

  RPVGRW.value = ZERO;

  // Test the input parameters.  PARAMS is not tested until ZGERFSX.

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
    xerbla('ZGESVXX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    zgeequb(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      zlaqge(N, N, A, LDA, R, C, ROWCND.value, COLCND.value, AMAX.value, EQUED);
      ROWEQU = lsame(EQUED.value, 'R') || lsame(EQUED.value, 'B');
      COLEQU = lsame(EQUED.value, 'C') || lsame(EQUED.value, 'B');
    }

    // If the scaling factors are not applied, set them to 1.0.

    if (!ROWEQU) {
      for (J = 1; J <= N; J++) {
        R[J] = 1.0;
      }
    }
    if (!COLEQU) {
      for (J = 1; J <= N; J++) {
        C[J] = 1.0;
      }
    }
  }

  // Scale the right-hand side.

  if (NOTRAN) {
    if (ROWEQU) zlascl2(N, NRHS, R, B, LDB);
  } else {
    if (COLEQU) zlascl2(N, NRHS, C, B, LDB);
  }

  if (NOFACT || EQUIL) {
    // Compute the LU factorization of A.

    zlacpy('Full', N, N, A, LDA, AF, LDAF);
    zgetrf(N, N, AF, LDAF, IPIV, INFO);

    // Return if INFO is non-zero.

    if (INFO.value > 0) {
      // Pivot in column INFO is exactly 0
      // Compute the reciprocal pivot growth factor of the
      // leading rank-deficient INFO columns of A.

      RPVGRW.value = zla_gerpvgrw(N, INFO.value, A, LDA, AF, LDAF);
      return;
    }
  }

  // Compute the reciprocal pivot growth factor RPVGRW.

  RPVGRW.value = zla_gerpvgrw(N, N, A, LDA, AF, LDAF);

  // Compute the solution matrix X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zgetrs(TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  zgerfsx(
      TRANS,
      EQUED.value,
      N,
      NRHS,
      A,
      LDA,
      AF,
      LDAF,
      IPIV,
      R,
      C,
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

  if (COLEQU && NOTRAN) {
    zlascl2(N, NRHS, C, X, LDX);
  } else if (ROWEQU && !NOTRAN) {
    zlascl2(N, NRHS, R, X, LDX);
  }
}
