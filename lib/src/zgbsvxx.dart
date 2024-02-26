import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgbequb.dart';
import 'package:lapack/src/zgbrfsx.dart';
import 'package:lapack/src/zgbtrf.dart';
import 'package:lapack/src/zgbtrs.dart';
import 'package:lapack/src/zla_gbrpvgrw.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaqgb.dart';
import 'package:lapack/src/zlascl2.dart';

void zgbsvxx(
  final String FACT,
  final String TRANS,
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> AFB_,
  final int LDAFB,
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
  final AB = AB_.dim(LDAB);
  final AFB = AFB_.dim(LDAFB);
  final IPIV = IPIV_.dim();
  final B = B_.dim(LDB);
  final X = X_.dim(LDX);
  final R = R_.dim();
  final C = C_.dim();
  final ERR_BNDS_NORM = ERR_BNDS_NORM_.dim(NRHS);
  final ERR_BNDS_COMP = ERR_BNDS_COMP_.dim(NRHS);
  final BERR = BERR_.dim();
  final PARAMS = PARAMS_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

  const ZERO = 0.0, ONE = 1.0;
  // const              FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 ;
  // const              RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 ;
  // const              CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 ;
  bool COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
  int I, J, KL = 0, KU = 0;
  double BIGNUM, RCMAX, RCMIN, SMLNUM;
  final ROWCND = Box(0.0), COLCND = Box(0.0), AMAX = Box(0.0);
  final INFEQU = Box(0);

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
  // pivot growth is set here, the rest is initialized in ZGBRFSX.

  RPVGRW.value = ZERO;

  // Test the input parameters.  PARAMS is not tested until ZGERFSX.

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
    xerbla('ZGBSVXX', -INFO.value);
    return;
  }

  if (EQUIL) {
    // Compute row and column scalings to equilibrate the matrix A.

    zgbequb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU);
    if (INFEQU.value == 0) {
      // Equilibrate the matrix.

      zlaqgb(N, N, KL, KU, AB, LDAB, R, C, ROWCND.value, COLCND.value,
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

  // Scale the right-hand side.

  if (NOTRAN) {
    if (ROWEQU) zlascl2(N, NRHS, R, B, LDB);
  } else {
    if (COLEQU) zlascl2(N, NRHS, C, B, LDB);
  }

  if (NOFACT || EQUIL) {
    // Compute the LU factorization of A.

    for (J = 1; J <= N; J++) {
      for (I = KL + 1; J <= 2 * KL + KU + 1; J++) {
        AFB[I][J] = AB[I - KL][J];
      } // 30
    } // 40
    zgbtrf(N, N, KL, KU, AFB, LDAFB, IPIV, INFO);

    // Return if INFO.value is non-zero.

    if (INFO.value > 0) {
      // Pivot in column INFO.value is exactly 0
      // Compute the reciprocal pivot growth factor of the
      // leading rank-deficient INFO.value columns of A.

      RPVGRW.value = zla_gbrpvgrw(N, KL, KU, INFO.value, AB, LDAB, AFB, LDAFB);
      return;
    }
  }

  // Compute the reciprocal pivot growth factor RPVGRW.value.

  RPVGRW.value = zla_gbrpvgrw(N, KL, KU, N, AB, LDAB, AFB, LDAFB);

  // Compute the solution matrix X.

  zlacpy('Full', N, NRHS, B, LDB, X, LDX);
  zgbtrs(TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO);

  // Use iterative refinement to improve the computed solution and
  // compute error bounds and backward error estimates for it.

  zgbrfsx(
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
      RWORK,
      INFO);

  // Scale solutions.

  if (COLEQU && NOTRAN) {
    zlascl2(N, NRHS, C, X, LDX);
  } else if (ROWEQU && !NOTRAN) {
    zlascl2(N, NRHS, R, X, LDX);
  }
}
