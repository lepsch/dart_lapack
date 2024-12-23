// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zgbmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/chla_transtype.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xblas/blas_zgbmv2_x.dart';
import 'package:dart_lapack/src/xblas/blas_zgbmv_x.dart';
import 'package:dart_lapack/src/zgbtrs.dart';
import 'package:dart_lapack/src/zla_gbamv.dart';
import 'package:dart_lapack/src/zla_lin_berr.dart';
import 'package:dart_lapack/src/zla_wwaddw.dart';

void zla_gbrfsx_extended(
  final int PREC_TYPE,
  final int TRANS_TYPE,
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> AFB_,
  final int LDAFB,
  final Array<int> IPIV_,
  final bool COLEQU,
  final Array<double> C_,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> Y_,
  final int LDY,
  final Array<double> BERR_OUT_,
  final int N_NORMS,
  final Matrix<double> ERR_BNDS_NORM_,
  final Matrix<double> ERR_BNDS_COMP_,
  final Array<Complex> RES_,
  final Array<double> AYB_,
  final Array<Complex> DY_,
  final Array<Complex> Y_TAIL_,
  final double RCOND,
  final int ITHRESH,
  final double RTHRESH,
  final double DZ_UB,
  final bool IGNORE_CWISE,
  final Box<int> INFO,
) {
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAFB);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final Y = Y_.having(ld: LDY);
  final RES = RES_.having();
  final DY = DY_.having();
  final Y_TAIL = Y_TAIL_.having();
  final C = C_.having();
  final BERR_OUT = BERR_OUT_.having();
  final AYB = AYB_.having();
  final ERR_BNDS_NORM = ERR_BNDS_NORM_.having(ld: NRHS);
  final ERR_BNDS_COMP = ERR_BNDS_COMP_.having(ld: NRHS);
  String TRANS;
  int CNT, I, J, M, X_STATE = 0, Z_STATE = 0, Y_PREC_STATE;
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
      DXRATMAX = 0,
      DZRATMAX = 0,
      DX_X = 0,
      DZ_Z = 0,
      FINAL_DX_X = 0,
      FINAL_DZ_Z = 0,
      EPS,
      HUGEVAL,
      INCR_THRESH;
  bool INCR_PREC;

  const UNSTABLE_STATE = 0, WORKING_STATE = 1, CONV_STATE = 2, NOPROG_STATE = 3;
  const BASE_RESIDUAL = 0, EXTRA_RESIDUAL = 1, EXTRA_Y = 2;
  const LA_LINRX_ERR_I = 2;

  if (INFO.value != 0) return;
  TRANS = chla_transtype(TRANS_TYPE);
  EPS = dlamch('Epsilon');
  HUGEVAL = dlamch('Overflow');
  // Force HUGEVAL to Inf
  HUGEVAL *= HUGEVAL;
  // Using HUGEVAL may lead to spurious underflows.
  INCR_THRESH = N * EPS;
  M = KL + KU + 1;

  for (J = 1; J <= NRHS; J++) {
    Y_PREC_STATE = EXTRA_RESIDUAL;
    if (Y_PREC_STATE == EXTRA_Y) {
      for (I = 1; I <= N; I++) {
        Y_TAIL[I] = Complex.zero;
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

      zcopy(N, B(1, J).asArray(), 1, RES, 1);
      if (Y_PREC_STATE == BASE_RESIDUAL) {
        zgbmv(TRANS, M, N, KL, KU, -Complex.one, AB, LDAB, Y(1, J).asArray(), 1,
            Complex.one, RES, 1);
      } else if (Y_PREC_STATE == EXTRA_RESIDUAL) {
        blas_zgbmv_x(TRANS_TYPE, N, N, KL, KU, -Complex.one, AB, LDAB,
            Y(1, J).asArray(), 1, Complex.one, RES, 1, PREC_TYPE);
      } else {
        blas_zgbmv2_x(TRANS_TYPE, N, N, KL, KU, -Complex.one, AB, LDAB,
            Y(1, J).asArray(), Y_TAIL, 1, Complex.one, RES, 1, PREC_TYPE);
      }

      // XXX: RES is no longer needed.
      zcopy(N, RES, 1, DY, 1);
      zgbtrs(TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, DY.asMatrix(N), N, INFO);

      // Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.

      NORMX = 0.0;
      NORMY = 0.0;
      NORMDX = 0.0;
      DZ_Z = 0.0;
      YMIN = HUGEVAL;

      for (I = 1; I <= N; I++) {
        YK = Y[I][J].cabs1();
        DYK = DY[I].cabs1();

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

      if (!IGNORE_CWISE &&
          YMIN * RCOND < INCR_THRESH * NORMY &&
          Y_PREC_STATE < EXTRA_Y) INCR_PREC = true;
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

      // Exit if both normwise and componentwise stopped working,
      // but if componentwise is unstable, let it go at least two
      // iterations.

      if (X_STATE != WORKING_STATE) {
        if (IGNORE_CWISE) break;
        if (Z_STATE == NOPROG_STATE || Z_STATE == CONV_STATE) break;
        if (Z_STATE == UNSTABLE_STATE && CNT > 1) break;
      }

      if (INCR_PREC) {
        INCR_PREC = false;
        Y_PREC_STATE++;
        for (I = 1; I <= N; I++) {
          Y_TAIL[I] = Complex.zero;
        }
      }

      PREVNORMDX = NORMDX;
      PREV_DZ_Z = DZ_Z;

      // Update solution.

      if (Y_PREC_STATE < EXTRA_Y) {
        zaxpy(N, Complex.one, DY, 1, Y(1, J).asArray(), 1);
      } else {
        zla_wwaddw(N, Y(1, J).asArray(), Y_TAIL, DY);
      }
    }

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

    zcopy(N, B(1, J).asArray(), 1, RES, 1);
    zgbmv(TRANS, N, N, KL, KU, -Complex.one, AB, LDAB, Y(1, J).asArray(), 1,
        Complex.one, RES, 1);

    for (I = 1; I <= N; I++) {
      AYB[I] = B[I][J].cabs1();
    }

    // Compute abs(op(A_s))*abs(Y) + abs(B_s).

    zla_gbamv(TRANS_TYPE, N, N, KL, KU, 1.0, AB, LDAB, Y(1, J).asArray(), 1,
        1.0, AYB, 1);

    zla_lin_berr(N, N, 1, RES.asMatrix(N), AYB.asMatrix(N), BERR_OUT(J));

    // End of loop for each RHS.
  }
}
