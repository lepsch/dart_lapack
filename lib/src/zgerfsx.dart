import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaprec.dart';
import 'package:lapack/src/ilatrans.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgecon.dart';
import 'package:lapack/src/zla_gercond_c.dart';
import 'package:lapack/src/zla_gercond_x.dart';
import 'package:lapack/src/zla_gerfsx_extended.dart';
import 'package:lapack/src/zlange.dart';

void zgerfsx(
  final String TRANS,
  final String EQUED,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
  final Array<int> IPIV_,
  final Array<double> R_,
  final Array<double> C_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> X_,
  final int LDX,
  final Box<double> RCOND,
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
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final R = R_.having();
  final C = C_.having();
  final BERR = BERR_.having();
  final PARAMS = PARAMS_.having();
  final ERR_BNDS_NORM = ERR_BNDS_NORM_.having(ld: NRHS);
  final ERR_BNDS_COMP = ERR_BNDS_COMP_.having(ld: NRHS);
  const ITREF_DEFAULT = 1.0;
  const ITHRESH_DEFAULT = 10.0;
  const COMPONENTWISE_DEFAULT = 1.0;
  const RTHRESH_DEFAULT = 0.5;
  const DZTHRESH_DEFAULT = 0.25;
  const LA_LINRX_ITREF_I = 1, LA_LINRX_ITHRESH_I = 2;
  const LA_LINRX_CWISE_I = 3;
  const LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2;
  const LA_LINRX_RCOND_I = 3;
  String NORM;
  bool ROWEQU, COLEQU, NOTRAN;
  int J, TRANS_TYPE, PREC_TYPE, REF_TYPE;
  int N_NORMS;
  double ANORM, RCOND_TMP;
  double ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG;
  bool IGNORE_CWISE;
  int ITHRESH;
  double RTHRESH, UNSTABLE_THRESH;

  // Check the input parameters.

  INFO.value = 0;
  TRANS_TYPE = ilatrans(TRANS);
  REF_TYPE = ITREF_DEFAULT.toInt();
  if (NPARAMS >= LA_LINRX_ITREF_I) {
    if (PARAMS[LA_LINRX_ITREF_I] < 0.0) {
      PARAMS[LA_LINRX_ITREF_I] = ITREF_DEFAULT;
    } else {
      REF_TYPE = PARAMS[LA_LINRX_ITREF_I].toInt();
    }
  }

  // Set default parameters.

  ILLRCOND_THRESH = N * dlamch('Epsilon');
  ITHRESH = ITHRESH_DEFAULT.toInt();
  RTHRESH = RTHRESH_DEFAULT;
  UNSTABLE_THRESH = DZTHRESH_DEFAULT;
  IGNORE_CWISE = COMPONENTWISE_DEFAULT == 0.0;

  if (NPARAMS >= LA_LINRX_ITHRESH_I) {
    if (PARAMS[LA_LINRX_ITHRESH_I] < 0.0) {
      PARAMS[LA_LINRX_ITHRESH_I] = ITHRESH.toDouble();
    } else {
      ITHRESH = PARAMS[LA_LINRX_ITHRESH_I].toInt();
    }
  }
  if (NPARAMS >= LA_LINRX_CWISE_I) {
    if (PARAMS[LA_LINRX_CWISE_I] < 0.0) {
      if (IGNORE_CWISE) {
        PARAMS[LA_LINRX_CWISE_I] = 0.0;
      } else {
        PARAMS[LA_LINRX_CWISE_I] = 1.0;
      }
    } else {
      IGNORE_CWISE = PARAMS[LA_LINRX_CWISE_I] == 0.0;
    }
  }
  if (REF_TYPE == 0 || N_ERR_BNDS == 0) {
    N_NORMS = 0;
  } else if (IGNORE_CWISE) {
    N_NORMS = 1;
  } else {
    N_NORMS = 2;
  }

  NOTRAN = lsame(TRANS, 'N');
  ROWEQU = lsame(EQUED, 'R') || lsame(EQUED, 'B');
  COLEQU = lsame(EQUED, 'C') || lsame(EQUED, 'B');

  // Test input parameters.

  if (TRANS_TYPE == -1) {
    INFO.value = -1;
  } else if (!ROWEQU && !COLEQU && !lsame(EQUED, 'N')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDAF < max(1, N)) {
    INFO.value = -8;
  } else if (LDB < max(1, N)) {
    INFO.value = -13;
  } else if (LDX < max(1, N)) {
    INFO.value = -15;
  }
  if (INFO.value != 0) {
    xerbla('ZGERFSX', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0 || NRHS == 0) {
    RCOND.value = 1.0;
    for (J = 1; J <= NRHS; J++) {
      BERR[J] = 0.0;
      if (N_ERR_BNDS >= 1) {
        ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 1.0;
        ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 1.0;
      }
      if (N_ERR_BNDS >= 2) {
        ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = 0.0;
        ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = 0.0;
      }
      if (N_ERR_BNDS >= 3) {
        ERR_BNDS_NORM[J][LA_LINRX_RCOND_I] = 1.0;
        ERR_BNDS_COMP[J][LA_LINRX_RCOND_I] = 1.0;
      }
    }
    return;
  }

  // Default to failure.

  RCOND.value = 0.0;
  for (J = 1; J <= NRHS; J++) {
    BERR[J] = 1.0;
    if (N_ERR_BNDS >= 1) {
      ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 1.0;
      ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 1.0;
    }
    if (N_ERR_BNDS >= 2) {
      ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = 1.0;
      ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = 1.0;
    }
    if (N_ERR_BNDS >= 3) {
      ERR_BNDS_NORM[J][LA_LINRX_RCOND_I] = 0.0;
      ERR_BNDS_COMP[J][LA_LINRX_RCOND_I] = 0.0;
    }
  }

  // Compute the norm of A and the reciprocal of the condition
  // number of A.

  if (NOTRAN) {
    NORM = 'I';
  } else {
    NORM = '1';
  }
  ANORM = zlange(NORM, N, N, A, LDA, RWORK);
  zgecon(NORM, N, AF, LDAF, ANORM, RCOND, WORK, RWORK, INFO);

  // Perform refinement on each right-hand side

  if (REF_TYPE != 0) {
    PREC_TYPE = ilaprec('E');

    if (NOTRAN) {
      zla_gerfsx_extended(
          PREC_TYPE,
          TRANS_TYPE,
          N,
          NRHS,
          A,
          LDA,
          AF,
          LDAF,
          IPIV,
          COLEQU,
          C,
          B,
          LDB,
          X,
          LDX,
          BERR,
          N_NORMS,
          ERR_BNDS_NORM,
          ERR_BNDS_COMP,
          WORK,
          RWORK,
          WORK(N + 1),
          RWORK.cast<Complex>(),
          RCOND.value,
          ITHRESH,
          RTHRESH,
          UNSTABLE_THRESH,
          IGNORE_CWISE,
          INFO);
    } else {
      zla_gerfsx_extended(
          PREC_TYPE,
          TRANS_TYPE,
          N,
          NRHS,
          A,
          LDA,
          AF,
          LDAF,
          IPIV,
          ROWEQU,
          R,
          B,
          LDB,
          X,
          LDX,
          BERR,
          N_NORMS,
          ERR_BNDS_NORM,
          ERR_BNDS_COMP,
          WORK,
          RWORK,
          WORK(N + 1),
          RWORK.cast<Complex>(),
          RCOND.value,
          ITHRESH,
          RTHRESH,
          UNSTABLE_THRESH,
          IGNORE_CWISE,
          INFO);
    }
  }

  ERR_LBND = max(10.0, sqrt(N)) * dlamch('Epsilon');
  if (N_ERR_BNDS >= 1 && N_NORMS >= 1) {
    // Compute scaled normwise condition number cond(A*C).

    if (COLEQU && NOTRAN) {
      RCOND_TMP = zla_gercond_c(
          TRANS, N, A, LDA, AF, LDAF, IPIV, C, true, INFO, WORK, RWORK);
    } else if (ROWEQU && !NOTRAN) {
      RCOND_TMP = zla_gercond_c(
          TRANS, N, A, LDA, AF, LDAF, IPIV, R, true, INFO, WORK, RWORK);
    } else {
      RCOND_TMP = zla_gercond_c(
          TRANS, N, A, LDA, AF, LDAF, IPIV, C, false, INFO, WORK, RWORK);
    }
    for (J = 1; J <= NRHS; J++) {
      // Cap the error at 1.0.

      if (N_ERR_BNDS >= LA_LINRX_ERR_I &&
          ERR_BNDS_NORM[J][LA_LINRX_ERR_I] > 1.0) {
        ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = 1.0;
      }

      // Threshold the error (see LAWN).

      if (RCOND_TMP < ILLRCOND_THRESH) {
        ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = 1.0;
        ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 0.0;
        if (INFO.value <= N) INFO.value = N + J;
      } else if (ERR_BNDS_NORM[J][LA_LINRX_ERR_I] < ERR_LBND) {
        ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = ERR_LBND;
        ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 1.0;
      }

      // Save the condition number.

      if (N_ERR_BNDS >= LA_LINRX_RCOND_I) {
        ERR_BNDS_NORM[J][LA_LINRX_RCOND_I] = RCOND_TMP;
      }
    }
  }

  if (N_ERR_BNDS >= 1 && N_NORMS >= 2) {
    // Compute componentwise condition number cond(A*diag(Y(:,J))) for
    // each right-hand side using the current solution as an estimate of
    // the true solution.  If the componentwise error estimate is too
    // large, then the solution is a lousy estimate of truth and the
    // estimated RCOND.value may be too optimistic.  To avoid misleading users,
    // the inverse condition number is set to 0.0 when the estimated
    // cwise error is at least CWISE_WRONG.

    CWISE_WRONG = sqrt(dlamch('Epsilon'));
    for (J = 1; J <= NRHS; J++) {
      if (ERR_BNDS_COMP[J][LA_LINRX_ERR_I] < CWISE_WRONG) {
        RCOND_TMP = zla_gercond_x(TRANS, N, A, LDA, AF, LDAF, IPIV,
            X(1, J).asArray(), INFO, WORK, RWORK);
      } else {
        RCOND_TMP = 0.0;
      }

      // Cap the error at 1.0.

      if (N_ERR_BNDS >= LA_LINRX_ERR_I &&
          ERR_BNDS_COMP[J][LA_LINRX_ERR_I] > 1.0) {
        ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = 1.0;
      }

      // Threshold the error (see LAWN).

      if (RCOND_TMP < ILLRCOND_THRESH) {
        ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = 1.0;
        ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 0.0;
        if (PARAMS[LA_LINRX_CWISE_I] == 1.0 && INFO.value < N + J) {
          INFO.value = N + J;
        }
      } else if (ERR_BNDS_COMP[J][LA_LINRX_ERR_I] < ERR_LBND) {
        ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = ERR_LBND;
        ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 1.0;
      }

      // Save the condition number.

      if (N_ERR_BNDS >= LA_LINRX_RCOND_I) {
        ERR_BNDS_COMP[J][LA_LINRX_RCOND_I] = RCOND_TMP;
      }
    }
  }
}
