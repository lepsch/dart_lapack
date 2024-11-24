// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dladiv.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlatrs.dart';
import 'package:lapack/src/matrix.dart';

void dlaein(
  final bool RIGHTV,
  final bool NOINIT,
  final int N,
  final Matrix<double> H_,
  final int LDH,
  final double WR,
  final double WI,
  final Array<double> VR_,
  final Array<double> VI_,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final double EPS3,
  final double SMLNUM,
  final double BIGNUM,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final VR = VR_.having();
  final VI = VI_.having();
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0, TENTH = 1.0e-1;
  String NORMIN, TRANS;
  int I, I1, I2, I3, ITS, J;
  double ABSBII,
      ABSBJJ,
      EI,
      EJ,
      GROWTO,
      NORM,
      NRMSML,
      REC,
      ROOTN,
      TEMP,
      VCRIT,
      VMAX,
      VNORM,
      W,
      W1,
      X,
      XI,
      XR,
      Y;
  final IERR = Box(0);
  final SCALE = Box(0.0);

  INFO.value = 0;

  // GROWTO is the threshold used in the acceptance test for an
  // eigenvector.

  ROOTN = sqrt(N);
  GROWTO = TENTH / ROOTN;
  NRMSML = max(ONE, EPS3 * ROOTN) * SMLNUM;

  // Form B = H - (WR,WI)*I (except that the subdiagonal elements and
  // the imaginary parts of the diagonal elements are not stored).

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= J - 1; I++) {
      B[I][J] = H[I][J];
    }
    B[J][J] = H[J][J] - WR;
  }

  if (WI == ZERO) {
    // Real eigenvalue.

    if (NOINIT) {
      // Set initial vector.

      for (I = 1; I <= N; I++) {
        VR[I] = EPS3;
      }
    } else {
      // Scale supplied initial vector.

      VNORM = dnrm2(N, VR, 1);
      dscal(N, (EPS3 * ROOTN) / max(VNORM, NRMSML), VR, 1);
    }

    if (RIGHTV) {
      // LU decomposition with partial pivoting of B, replacing zero
      // pivots by EPS3.

      for (I = 1; I <= N - 1; I++) {
        EI = H[I + 1][I];
        if (B[I][I].abs() < EI.abs()) {
          // Interchange rows and eliminate.

          X = B[I][I] / EI;
          B[I][I] = EI;
          for (J = I + 1; J <= N; J++) {
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - X * TEMP;
            B[I][J] = TEMP;
          }
        } else {
          // Eliminate without interchange.

          if (B[I][I] == ZERO) B[I][I] = EPS3;
          X = EI / B[I][I];
          if (X != ZERO) {
            for (J = I + 1; J <= N; J++) {
              B[I + 1][J] -= X * B[I][J];
            }
          }
        }
      }
      if (B[N][N] == ZERO) B[N][N] = EPS3;

      TRANS = 'N';
    } else {
      // UL decomposition with partial pivoting of B, replacing zero
      // pivots by EPS3.

      for (J = N; J >= 2; J--) {
        EJ = H[J][J - 1];
        if (B[J][J].abs() < EJ.abs()) {
          // Interchange columns and eliminate.

          X = B[J][J] / EJ;
          B[J][J] = EJ;
          for (I = 1; I <= J - 1; I++) {
            TEMP = B[I][J - 1];
            B[I][J - 1] = B[I][J] - X * TEMP;
            B[I][J] = TEMP;
          }
        } else {
          // Eliminate without interchange.

          if (B[J][J] == ZERO) B[J][J] = EPS3;
          X = EJ / B[J][J];
          if (X != ZERO) {
            for (I = 1; I <= J - 1; I++) {
              B[I][J - 1] -= X * B[I][J];
            }
          }
        }
      }
      if (B[1][1] == ZERO) B[1][1] = EPS3;

      TRANS = 'T';
    }

    NORMIN = 'N';
    var eigenvectorFound = false;
    for (ITS = 1; ITS <= N; ITS++) {
      // Solve U*x = scale*v for a right eigenvector
      //    or U**T*x = scale*v for a left eigenvector,
      // overwriting x on v.

      dlatrs(
          'Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, VR, SCALE, WORK, IERR);
      NORMIN = 'Y';

      // Test for sufficient growth in the norm of v.

      VNORM = dasum(N, VR, 1);
      if (VNORM >= GROWTO * SCALE.value) {
        eigenvectorFound = true;
        break;
      }

      // Choose new orthogonal starting vector and try again.

      TEMP = EPS3 / (ROOTN + ONE);
      VR[1] = EPS3;
      for (I = 2; I <= N; I++) {
        VR[I] = TEMP;
      }
      VR[N - ITS + 1] -= EPS3 * ROOTN;
    }

    // Failure to find eigenvector in N iterations.
    if (!eigenvectorFound) {
      INFO.value = 1;
    }

    // Normalize eigenvector.

    I = idamax(N, VR, 1);
    dscal(N, ONE / VR[I].abs(), VR, 1);
  } else {
    // Complex eigenvalue.

    if (NOINIT) {
      // Set initial vector.

      for (I = 1; I <= N; I++) {
        VR[I] = EPS3;
        VI[I] = ZERO;
      }
    } else {
      // Scale supplied initial vector.

      NORM = dlapy2(dnrm2(N, VR, 1), dnrm2(N, VI, 1));
      REC = (EPS3 * ROOTN) / max(NORM, NRMSML);
      dscal(N, REC, VR, 1);
      dscal(N, REC, VI, 1);
    }

    if (RIGHTV) {
      // LU decomposition with partial pivoting of B, replacing zero
      // pivots by EPS3.

      // The imaginary part of the (i,j)-th element of U is stored in
      // B[j+1][i].

      B[2][1] = -WI;
      for (I = 2; I <= N; I++) {
        B[I + 1][1] = ZERO;
      }

      for (I = 1; I <= N - 1; I++) {
        ABSBII = dlapy2(B[I][I], B[I + 1][I]);
        EI = H[I + 1][I];
        if (ABSBII < EI.abs()) {
          // Interchange rows and eliminate.

          XR = B[I][I] / EI;
          XI = B[I + 1][I] / EI;
          B[I][I] = EI;
          B[I + 1][I] = ZERO;
          for (J = I + 1; J <= N; J++) {
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - XR * TEMP;
            B[J + 1][I + 1] = B[J + 1][I] - XI * TEMP;
            B[I][J] = TEMP;
            B[J + 1][I] = ZERO;
          }
          B[I + 2][I] = -WI;
          B[I + 1][I + 1] -= XI * WI;
          B[I + 2][I + 1] += XR * WI;
        } else {
          // Eliminate without interchanging rows.

          if (ABSBII == ZERO) {
            B[I][I] = EPS3;
            B[I + 1][I] = ZERO;
            ABSBII = EPS3;
          }
          EI = (EI / ABSBII) / ABSBII;
          XR = B[I][I] * EI;
          XI = -B[I + 1][I] * EI;
          for (J = I + 1; J <= N; J++) {
            B[I + 1][J] -= XR * B[I][J] - XI * B[J + 1][I];
            B[J + 1][I + 1] = -XR * B[J + 1][I] - XI * B[I][J];
          }
          B[I + 2][I + 1] -= WI;
        }

        // Compute 1-norm of offdiagonal elements of i-th row.

        WORK[I] = dasum(N - I, B(I, I + 1).asArray(), LDB) +
            dasum(N - I, B(I + 2, I).asArray(), 1);
      }
      if (B[N][N] == ZERO && B[N + 1][N] == ZERO) B[N][N] = EPS3;
      WORK[N] = ZERO;

      I1 = N;
      I2 = 1;
      I3 = -1;
    } else {
      // UL decomposition with partial pivoting of conjg(B),
      // replacing zero pivots by EPS3.

      // The imaginary part of the (i,j)-th element of U is stored in
      // B[j+1][i].

      B[N + 1][N] = WI;
      for (J = 1; J <= N - 1; J++) {
        B[N + 1][J] = ZERO;
      }

      for (J = N; J >= 2; J--) {
        EJ = H[J][J - 1];
        ABSBJJ = dlapy2(B[J][J], B[J + 1][J]);
        if (ABSBJJ < EJ.abs()) {
          // Interchange columns and eliminate

          XR = B[J][J] / EJ;
          XI = B[J + 1][J] / EJ;
          B[J][J] = EJ;
          B[J + 1][J] = ZERO;
          for (I = 1; I <= J - 1; I++) {
            TEMP = B[I][J - 1];
            B[I][J - 1] = B[I][J] - XR * TEMP;
            B[J][I] = B[J + 1][I] - XI * TEMP;
            B[I][J] = TEMP;
            B[J + 1][I] = ZERO;
          }
          B[J + 1][J - 1] = WI;
          B[J - 1][J - 1] += XI * WI;
          B[J][J - 1] -= XR * WI;
        } else {
          // Eliminate without interchange.

          if (ABSBJJ == ZERO) {
            B[J][J] = EPS3;
            B[J + 1][J] = ZERO;
            ABSBJJ = EPS3;
          }
          EJ = (EJ / ABSBJJ) / ABSBJJ;
          XR = B[J][J] * EJ;
          XI = -B[J + 1][J] * EJ;
          for (I = 1; I <= J - 1; I++) {
            B[I][J - 1] -= XR * B[I][J] - XI * B[J + 1][I];
            B[J][I] = -XR * B[J + 1][I] - XI * B[I][J];
          }
          B[J][J - 1] += WI;
        }

        // Compute 1-norm of offdiagonal elements of j-th column.

        WORK[J] = dasum(J - 1, B(1, J).asArray(), 1) +
            dasum(J - 1, B(J + 1, 1).asArray(), LDB);
      }
      if (B[1][1] == ZERO && B[2][1] == ZERO) B[1][1] = EPS3;
      WORK[1] = ZERO;

      I1 = 1;
      I2 = N;
      I3 = 1;
    }

    var eigenvectorFound = false;
    for (ITS = 1; ITS <= N; ITS++) {
      SCALE.value = ONE;
      VMAX = ONE;
      VCRIT = BIGNUM;

      // Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
      // or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector,
      // overwriting (xr,xi) on (vr,vi).

      for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
        if (WORK[I] > VCRIT) {
          REC = ONE / VMAX;
          dscal(N, REC, VR, 1);
          dscal(N, REC, VI, 1);
          SCALE.value *= REC;
          VMAX = ONE;
          VCRIT = BIGNUM;
        }

        XR = VR[I];
        XI = VI[I];
        if (RIGHTV) {
          for (J = I + 1; J <= N; J++) {
            XR -= B[I][J] * VR[J] - B[J + 1][I] * VI[J];
            XI -= B[I][J] * VI[J] + B[J + 1][I] * VR[J];
          }
        } else {
          for (J = 1; J <= I - 1; J++) {
            XR -= B[J][I] * VR[J] - B[I + 1][J] * VI[J];
            XI -= B[J][I] * VI[J] + B[I + 1][J] * VR[J];
          }
        }

        W = B[I][I].abs() + B[I + 1][I].abs();
        if (W > SMLNUM) {
          if (W < ONE) {
            W1 = XR.abs() + XI.abs();
            if (W1 > W * BIGNUM) {
              REC = ONE / W1;
              dscal(N, REC, VR, 1);
              dscal(N, REC, VI, 1);
              XR = VR[I];
              XI = VI[I];
              SCALE.value *= REC;
              VMAX *= REC;
            }
          }

          // Divide by diagonal element of B.

          dladiv(XR, XI, B[I][I], B[I + 1][I], VR(I), VI(I));
          VMAX = max(VR[I].abs() + VI[I].abs(), VMAX);
          VCRIT = BIGNUM / VMAX;
        } else {
          for (J = 1; J <= N; J++) {
            VR[J] = ZERO;
            VI[J] = ZERO;
          }
          VR[I] = ONE;
          VI[I] = ONE;
          SCALE.value = ZERO;
          VMAX = ONE;
          VCRIT = BIGNUM;
        }
      }

      // Test for sufficient growth in the norm of (VR,VI).

      VNORM = dasum(N, VR, 1) + dasum(N, VI, 1);
      if (VNORM >= GROWTO * SCALE.value) {
        eigenvectorFound = true;
        break;
      }

      // Choose a new orthogonal starting vector and try again.

      Y = EPS3 / (ROOTN + ONE);
      VR[1] = EPS3;
      VI[1] = ZERO;

      for (I = 2; I <= N; I++) {
        VR[I] = Y;
        VI[I] = ZERO;
      }
      VR[N - ITS + 1] -= EPS3 * ROOTN;
    }

    // Failure to find eigenvector in N iterations
    if (!eigenvectorFound) {
      INFO.value = 1;
    }

    // Normalize eigenvector.

    VNORM = ZERO;
    for (I = 1; I <= N; I++) {
      VNORM = max(VNORM, VR[I].abs() + VI[I].abs());
    }
    dscal(N, ONE / VNORM, VR, 1);
    dscal(N, ONE / VNORM, VI, 1);
  }
}
