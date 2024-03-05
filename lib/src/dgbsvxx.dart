import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbequb.dart';
import 'package:lapack/src/dgbrfsx.dart';
import 'package:lapack/src/dgbtrf.dart';
import 'package:lapack/src/dgbtrs.dart';
import 'package:lapack/src/dla_gbrpvgrw.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaqgb.dart';
import 'package:lapack/src/dlascl2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgbsvxx(
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
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAFB);
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
  final IWORK = IWORK_.having();

  const ZERO = 0.0, ONE = 1.0;
  // const FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3;
  // const RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6;
  // const CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9;
  bool COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
  int I, J;
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
  // pivot growth is set here, the rest is initialized in DGBRFSX.

  RPVGRW.value = ZERO;

  // Test the input parameters.  PARAMS is not tested until DGBRFSX.

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
        // 10
        RCMIN = min(RCMIN, R[J]);
        RCMAX = max(RCMAX, R[J]);
      } // 10
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
        // 20
        RCMIN = min(RCMIN, C[J]);
        RCMAX = max(RCMAX, C[J]);
      } // 20
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
        INFO.value = -15;
      } else if (LDX < max(1, N)) {
        INFO.value = -16;
      }
    }
  }

  if (INFO.value != 0) {
    xerbla('DGBSVXX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    dgbequb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      dlaqgb(N, N, KL, KU, AB, LDAB, R, C, ROWCND.value, COLCND.value,
          AMAX.value, EQUED);
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

  // Scale the right hand side.

  if (NOTRAN) {
    if (ROWEQU) dlascl2(N, NRHS, R, B, LDB);
  } else {
    if (COLEQU) dlascl2(N, NRHS, C, B, LDB);
  }

  if (NOFACT || EQUIL) {
    // Compute the LU factorization of A.

    for (J = 1; J <= N; J++) {
      for (I = KL + 1; I <= 2 * KL + KU + 1; I++) {
        AFB[I][J] = AB[I - KL][J];
      } // 30
    } // 40
    dgbtrf(N, N, KL, KU, AFB, LDAFB, IPIV, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      // Pivot in column INFO.value is exactly 0
      // Compute the reciprocal pivot growth factor of the
      // leading rank-deficient INFO.value columns of A.

      RPVGRW.value = dla_gbrpvgrw(N, KL, KU, INFO.value, AB, LDAB, AFB, LDAFB);
      return;
    }
  }

  // Compute the reciprocal pivot growth factor RPVGRW.value.

  RPVGRW.value = dla_gbrpvgrw(N, KL, KU, N, AB, LDAB, AFB, LDAFB);

  // Compute the solution matrix X.

  dlacpy('Full', N, NRHS, B, LDB, X, LDX);
  dgbtrs(TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  dgbrfsx(
      TRANS,
      EQUED.value,
      N,
      KL,
      KU,
      NRHS,
      AB,
      LDAB,
      AFB,
      LDAFB,
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
      IWORK,
      INFO);

  // Scale solutions.

  if (COLEQU && NOTRAN) {
    dlascl2(N, NRHS, C, X, LDX);
  } else if (ROWEQU && !NOTRAN) {
    dlascl2(N, NRHS, R, X, LDX);
  }
}
