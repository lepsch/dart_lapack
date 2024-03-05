import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dsymv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dla_lin_berr.dart';
import 'package:lapack/src/dla_syamv.dart';
import 'package:lapack/src/dla_wwaddw.dart';
import 'package:lapack/src/dpotrs.dart';
import 'package:lapack/src/ilauplo.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xblas/blas_dsymv2_x.dart';
import 'package:lapack/src/xblas/blas_dsymv_x.dart';

void dla_porfsx_extended(
  final int PREC_TYPE,
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AF_,
  final int LDAF,
  final bool COLEQU,
  final Array<double> C_,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Y_,
  final int LDY,
  final Array<double> BERR_OUT_,
  final int N_NORMS,
  final Matrix<double> ERR_BNDS_NORM_,
  final Matrix<double> ERR_BNDS_COMP_,
  final Array<double> RES_,
  final Array<double> AYB_,
  final Array<double> DY_,
  final Array<double> Y_TAIL_,
  final double RCOND,
  final int ITHRESH,
  final double RTHRESH,
  final double DZ_UB,
  final bool IGNORE_CWISE,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final B = B_.having(ld: LDB);
  final Y = Y_.having(ld: LDY);
  final ERR_BNDS_NORM = ERR_BNDS_NORM_.having(ld: NRHS);
  final ERR_BNDS_COMP = ERR_BNDS_COMP_.having(ld: NRHS);
  final C = C_.having();
  final BERR_OUT = BERR_OUT_.having();
  final RES = RES_.having();
  final AYB = AYB_.having();
  final DY = DY_.having();
  final Y_TAIL = Y_TAIL_.having();
  int UPLO2, CNT, I, J, X_STATE, Z_STATE;
  double YK,
      DYK,
      YMIN,
      NORMY,
      NORMX,
      NORMDX,
      DXRAT,
      DZRAT,
      PREVNORMDX,
      PREV_DZ_Z,
      DXRATMAX,
      DZRATMAX,
      DX_X,
      DZ_Z,
      FINAL_DX_X,
      FINAL_DZ_Z,
      EPS,
      HUGEVAL,
      INCR_THRESH;
  bool INCR_PREC;
  int Y_PREC_STATE;
  const UNSTABLE_STATE = 0, WORKING_STATE = 1, CONV_STATE = 2, NOPROG_STATE = 3;
  const BASE_RESIDUAL = 0, EXTRA_RESIDUAL = 1, EXTRA_Y = 2;
  // const FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3;
  // const RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6;
  // const CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9;
  // const LA_LINRX_ITREF_I = 1, LA_LINRX_ITHRESH_I = 2;
  // const LA_LINRX_CWISE_I = 3;
  // const LA_LINRX_TRUST_I = 1;
  const LA_LINRX_ERR_I = 2;
  // const LA_LINRX_RCOND_I = 3;

  if (INFO.value != 0) return;
  EPS = dlamch('Epsilon');
  HUGEVAL = dlamch('Overflow');
  // Force HUGEVAL to Inf
  HUGEVAL = HUGEVAL * HUGEVAL;
  // Using HUGEVAL may lead to spurious underflows.
  INCR_THRESH = N.toDouble() * EPS;

  if (lsame(UPLO, 'L')) {
    UPLO2 = ilauplo('L');
  } else {
    UPLO2 = ilauplo('U');
  }

  for (J = 1; J <= NRHS; J++) {
    Y_PREC_STATE = EXTRA_RESIDUAL;
    if (Y_PREC_STATE == EXTRA_Y) {
      for (I = 1; I <= N; I++) {
        Y_TAIL[I] = 0.0;
      }
    }

    DXRAT = 0.0;
    DXRATMAX = 0.0;
    DZRAT = 0.0;
    DZRATMAX = 0.0;
    FINAL_DX_X = HUGEVAL;
    FINAL_DZ_Z = HUGEVAL;
    PREVNORMDX = HUGEVAL;
    PREV_DZ_Z = HUGEVAL;
    DZ_Z = HUGEVAL;
    DX_X = HUGEVAL;

    X_STATE = WORKING_STATE;
    Z_STATE = UNSTABLE_STATE;
    INCR_PREC = false;

    for (CNT = 1; CNT <= ITHRESH; CNT++) {
      // Compute residual RES = B_s - op(A_s) * Y,
      //     op(A) = A, A**T, or A**H depending on TRANS (and type).

      dcopy(N, B(1, J).asArray(), 1, RES, 1);
      if (Y_PREC_STATE == BASE_RESIDUAL) {
        dsymv(UPLO, N, -1.0, A, LDA, Y(1, J).asArray(), 1, 1.0, RES, 1);
      } else if (Y_PREC_STATE == EXTRA_RESIDUAL) {
        blas_dsymv_x(UPLO2, N, -1.0, A, LDA, Y(1, J).asArray(), 1, 1.0, RES, 1,
            PREC_TYPE);
      } else {
        blas_dsymv2_x(UPLO2, N, -1.0, A, LDA, Y(1, J).asArray(), Y_TAIL, 1, 1.0,
            RES, 1, PREC_TYPE);
      }

      // XXX: RES is no longer needed.
      dcopy(N, RES, 1, DY, 1);
      dpotrs(UPLO, N, 1, AF, LDAF, DY.asMatrix(N), N, INFO);

      // Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.

      NORMX = 0.0;
      NORMY = 0.0;
      NORMDX = 0.0;
      DZ_Z = 0.0;
      YMIN = HUGEVAL;

      for (I = 1; I <= N; I++) {
        YK = Y[I][J].abs();
        DYK = DY[I].abs();

        if (YK != 0.0) {
          DZ_Z = max(DZ_Z, DYK / YK);
        } else if (DYK != 0.0) {
          DZ_Z = HUGEVAL;
        }

        YMIN = min(YMIN, YK);

        NORMY = max(NORMY, YK);

        if (COLEQU) {
          NORMX = max(NORMX, YK * C[I]);
          NORMDX = max(NORMDX, DYK * C[I]);
        } else {
          NORMX = NORMY;
          NORMDX = max(NORMDX, DYK);
        }
      }

      if (NORMX != 0.0) {
        DX_X = NORMDX / NORMX;
      } else if (NORMDX == 0.0) {
        DX_X = 0.0;
      } else {
        DX_X = HUGEVAL;
      }

      DXRAT = NORMDX / PREVNORMDX;
      DZRAT = DZ_Z / PREV_DZ_Z;

      // Check termination criteria.

      if (YMIN * RCOND < INCR_THRESH * NORMY && Y_PREC_STATE < EXTRA_Y) {
        INCR_PREC = true;
      }
      if (X_STATE == NOPROG_STATE && DXRAT <= RTHRESH) X_STATE = WORKING_STATE;
      if (X_STATE == WORKING_STATE) {
        if (DX_X <= EPS) {
          X_STATE = CONV_STATE;
        } else if (DXRAT > RTHRESH) {
          if (Y_PREC_STATE != EXTRA_Y) {
            INCR_PREC = true;
          } else {
            X_STATE = NOPROG_STATE;
          }
        } else {
          if (DXRAT > DXRATMAX) DXRATMAX = DXRAT;
        }
        if (X_STATE > WORKING_STATE) FINAL_DX_X = DX_X;
      }
      if (Z_STATE == UNSTABLE_STATE && DZ_Z <= DZ_UB) Z_STATE = WORKING_STATE;
      if (Z_STATE == NOPROG_STATE && DZRAT <= RTHRESH) Z_STATE = WORKING_STATE;
      if (Z_STATE == WORKING_STATE) {
        if (DZ_Z <= EPS) {
          Z_STATE = CONV_STATE;
        } else if (DZ_Z > DZ_UB) {
          Z_STATE = UNSTABLE_STATE;
          DZRATMAX = 0.0;
          FINAL_DZ_Z = HUGEVAL;
        } else if (DZRAT > RTHRESH) {
          if (Y_PREC_STATE != EXTRA_Y) {
            INCR_PREC = true;
          } else {
            Z_STATE = NOPROG_STATE;
          }
        } else {
          if (DZRAT > DZRATMAX) DZRATMAX = DZRAT;
        }
        if (Z_STATE > WORKING_STATE) FINAL_DZ_Z = DZ_Z;
      }
      if (X_STATE != WORKING_STATE &&
          (IGNORE_CWISE || Z_STATE != WORKING_STATE)) break;

      if (INCR_PREC) {
        INCR_PREC = false;
        Y_PREC_STATE = Y_PREC_STATE + 1;
        for (I = 1; I <= N; I++) {
          Y_TAIL[I] = 0.0;
        }
      }

      PREVNORMDX = NORMDX;
      PREV_DZ_Z = DZ_Z;

      // Update solution.

      if (Y_PREC_STATE < EXTRA_Y) {
        daxpy(N, 1.0, DY, 1, Y(1, J).asArray(), 1);
      } else {
        dla_wwaddw(N, Y(1, J).asArray(), Y_TAIL, DY);
      }
    }
    // Target of "IF (Z_STOP && X_STOP)".  Sun's f77 won't EXIT.

    // Set final_* when cnt hits ithresh.

    if (X_STATE == WORKING_STATE) FINAL_DX_X = DX_X;
    if (Z_STATE == WORKING_STATE) FINAL_DZ_Z = DZ_Z;

    // Compute error bounds.

    if (N_NORMS >= 1) {
      ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = FINAL_DX_X / (1 - DXRATMAX);
    }
    if (N_NORMS >= 2) {
      ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = FINAL_DZ_Z / (1 - DZRATMAX);
    }

    // Compute componentwise relative backward error from formula
    //     max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
    // where abs(Z) is the componentwise absolute value of the matrix
    // or vector Z.

    // Compute residual RES = B_s - op(A_s) * Y,
    //     op(A) = A, A**T, or A**H depending on TRANS (and type).

    dcopy(N, B(1, J).asArray(), 1, RES, 1);
    dsymv(UPLO, N, -1.0, A, LDA, Y(1, J).asArray(), 1, 1.0, RES, 1);

    for (I = 1; I <= N; I++) {
      AYB[I] = B[I][J].abs();
    }

    // Compute abs(op(A_s))*abs(Y) + abs(B_s).

    dla_syamv(UPLO2, N, 1.0, A, LDA, Y(1, J).asArray(), 1, 1.0, AYB, 1);

    dla_lin_berr(N, N, 1, RES.asMatrix(N), AYB.asMatrix(N), BERR_OUT(J));

    // End of loop for each RHS.
  }
}
