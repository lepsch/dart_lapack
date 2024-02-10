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
  final Matrix<double> H,
  final int LDH,
  final double WR,
  final double WI,
  final Array<double> VR,
  final Array<double> VI,
  final Matrix<double> B,
  final int LDB,
  final Array<double> WORK,
  final double EPS3,
  final double SMLNUM,
  final double BIGNUM,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
      SCALE = 0,
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

  INFO.value = 0;

  // GROWTO is the threshold used in the acceptance test for an
  // eigenvector.

  ROOTN = sqrt(N.toDouble());
  GROWTO = TENTH / ROOTN;
  NRMSML = max(ONE, EPS3 * ROOTN) * SMLNUM;

  // Form B = H - (WR,WI)*I (except that the subdiagonal elements and
  // the imaginary parts of the diagonal elements are not stored).

  for (J = 1; J <= N; J++) {
    // 20
    for (I = 1; I <= J - 1; I++) {
      // 10
      B[I][J] = H[I][J];
    } // 10
    B[J][J] = H[J][J] - WR;
  } // 20

  if (WI == ZERO) {
    // Real eigenvalue.

    if (NOINIT) {
      // Set initial vector.

      for (I = 1; I <= N; I++) {
        // 30
        VR[I] = EPS3;
      } // 30
    } else {
      // Scale supplied initial vector.

      VNORM = dnrm2(N, VR, 1);
      dscal(N, (EPS3 * ROOTN) / max(VNORM, NRMSML), VR, 1);
    }

    if (RIGHTV) {
      // LU decomposition with partial pivoting of B, replacing zero
      // pivots by EPS3.

      for (I = 1; I <= N - 1; I++) {
        // 60
        EI = H[I + 1][I];
        if ((B[I][I]).abs() < (EI).abs()) {
          // Interchange rows and eliminate.

          X = B[I][I] / EI;
          B[I][I] = EI;
          for (J = I + 1; J <= N; J++) {
            // 40
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - X * TEMP;
            B[I][J] = TEMP;
          } // 40
        } else {
          // Eliminate without interchange.

          if (B[I][I] == ZERO) B[I][I] = EPS3;
          X = EI / B[I][I];
          if (X != ZERO) {
            for (J = I + 1; J <= N; J++) {
              // 50
              B[I + 1][J] = B[I + 1][J] - X * B[I][J];
            } // 50
          }
        }
      } // 60
      if (B[N][N] == ZERO) B[N][N] = EPS3;

      TRANS = 'N';
    } else {
      // UL decomposition with partial pivoting of B, replacing zero
      // pivots by EPS3.

      for (J = N; J >= 2; J--) {
        // 90
        EJ = H[J][J - 1];
        if ((B[J][J]).abs() < (EJ).abs()) {
          // Interchange columns and eliminate.

          X = B[J][J] / EJ;
          B[J][J] = EJ;
          for (I = 1; I <= J - 1; I++) {
            // 70
            TEMP = B[I][J - 1];
            B[I][J - 1] = B[I][J] - X * TEMP;
            B[I][J] = TEMP;
          } // 70
        } else {
          // Eliminate without interchange.

          if (B[J][J] == ZERO) B[J][J] = EPS3;
          X = EJ / B[J][J];
          if (X != ZERO) {
            for (I = 1; I <= J - 1; I++) {
              // 80
              B[I][J - 1] = B[I][J - 1] - X * B[I][J];
            } // 80
          }
        }
      } // 90
      if (B[1][1] == ZERO) B[1][1] = EPS3;

      TRANS = 'T';
    }

    NORMIN = 'N';
    var eigenvectorFound = false;
    for (ITS = 1; ITS <= N; ITS++) {
      // 110

      // Solve U*x = scale*v for a right eigenvector
      //    or U**T*x = scale*v for a left eigenvector,
      // overwriting x on v.

      dlatrs(
          'Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, VR, SCALE, WORK, IERR);
      NORMIN = 'Y';

      // Test for sufficient growth in the norm of v.

      VNORM = dasum(N, VR, 1);
      if (VNORM < GROWTO * SCALE) {
        eigenvectorFound = true;
        break; // 120
      }

      // Choose new orthogonal starting vector and try again.

      TEMP = EPS3 / (ROOTN + ONE);
      VR[1] = EPS3;
      for (I = 2; I <= N; I++) {
        // 100
        VR[I] = TEMP;
      } // 100
      VR[N - ITS + 1] = VR[N - ITS + 1] - EPS3 * ROOTN;
    } // 110

    // Failure to find eigenvector in N iterations.
    if (!eigenvectorFound) {
      INFO.value = 1;
    } // 120

    // Normalize eigenvector.

    I = idamax(N, VR, 1);
    dscal(N, ONE / (VR[I]).abs(), VR, 1);
  } else {
    // Complex eigenvalue.

    if (NOINIT) {
      // Set initial vector.

      for (I = 1; I <= N; I++) {
        // 130
        VR[I] = EPS3;
        VI[I] = ZERO;
      } // 130
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
        // 140
        B[I + 1][1] = ZERO;
      } // 140

      for (I = 1; I <= N - 1; I++) {
        // 170
        ABSBII = dlapy2(B[I][I], B[I + 1][I]);
        EI = H[I + 1][I];
        if (ABSBII < (EI).abs()) {
          // Interchange rows and eliminate.

          XR = B[I][I] / EI;
          XI = B[I + 1][I] / EI;
          B[I][I] = EI;
          B[I + 1][I] = ZERO;
          for (J = I + 1; J <= N; J++) {
            // 150
            TEMP = B[I + 1][J];
            B[I + 1][J] = B[I][J] - XR * TEMP;
            B[J + 1][I + 1] = B[J + 1][I] - XI * TEMP;
            B[I][J] = TEMP;
            B[J + 1][I] = ZERO;
          } // 150
          B[I + 2][I] = -WI;
          B[I + 1][I + 1] = B[I + 1][I + 1] - XI * WI;
          B[I + 2][I + 1] = B[I + 2][I + 1] + XR * WI;
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
            // 160
            B[I + 1][J] = B[I + 1][J] - XR * B[I][J] + XI * B[J + 1][I];
            B[J + 1][I + 1] = -XR * B[J + 1][I] - XI * B[I][J];
          } // 160
          B[I + 2][I + 1] = B[I + 2][I + 1] - WI;
        }

        // Compute 1-norm of offdiagonal elements of i-th row.

        WORK[I] = dasum(N - I, B(I, I + 1).asArray(), LDB) +
            dasum(N - I, B(I + 2, I).asArray(), 1);
      } // 170
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
        // 180
        B[N + 1][J] = ZERO;
      } // 180

      for (J = N; J >= 2; J--) {
        // 210
        EJ = H[J][J - 1];
        ABSBJJ = dlapy2(B[J][J], B[J + 1][J]);
        if (ABSBJJ < (EJ).abs()) {
          // Interchange columns and eliminate

          XR = B[J][J] / EJ;
          XI = B[J + 1][J] / EJ;
          B[J][J] = EJ;
          B[J + 1][J] = ZERO;
          for (I = 1; I <= J - 1; I++) {
            // 190
            TEMP = B[I][J - 1];
            B[I][J - 1] = B[I][J] - XR * TEMP;
            B[J][I] = B[J + 1][I] - XI * TEMP;
            B[I][J] = TEMP;
            B[J + 1][I] = ZERO;
          } // 190
          B[J + 1][J - 1] = WI;
          B[J - 1][J - 1] = B[J - 1][J - 1] + XI * WI;
          B[J][J - 1] = B[J][J - 1] - XR * WI;
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
            // 200
            B[I][J - 1] = B[I][J - 1] - XR * B[I][J] + XI * B[J + 1][I];
            B[J][I] = -XR * B[J + 1][I] - XI * B[I][J];
          } // 200
          B[J][J - 1] = B[J][J - 1] + WI;
        }

        // Compute 1-norm of offdiagonal elements of j-th column.

        WORK[J] = dasum(J - 1, B(1, J).asArray(), 1) +
            dasum(J - 1, B(J + 1, 1).asArray(), LDB);
      } // 210
      if (B[1][1] == ZERO && B[2][1] == ZERO) B[1][1] = EPS3;
      WORK[1] = ZERO;

      I1 = 1;
      I2 = N;
      I3 = 1;
    }

    var eigenvectorFound = false;
    for (ITS = 1; ITS <= N; ITS++) {
      // 270
      SCALE = ONE;
      VMAX = ONE;
      VCRIT = BIGNUM;

      // Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
      // or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector,
      // overwriting (xr,xi) on (vr,vi).

      for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
        // 250

        if (WORK[I] > VCRIT) {
          REC = ONE / VMAX;
          dscal(N, REC, VR, 1);
          dscal(N, REC, VI, 1);
          SCALE = SCALE * REC;
          VMAX = ONE;
          VCRIT = BIGNUM;
        }

        XR = VR[I];
        XI = VI[I];
        if (RIGHTV) {
          for (J = I + 1; J <= N; J++) {
            // 220
            XR = XR - B[I][J] * VR[J] + B[J + 1][I] * VI[J];
            XI = XI - B[I][J] * VI[J] - B[J + 1][I] * VR[J];
          } // 220
        } else {
          for (J = 1; J <= I - 1; J++) {
            // 230
            XR = XR - B[J][I] * VR[J] + B[I + 1][J] * VI[J];
            XI = XI - B[J][I] * VI[J] - B[I + 1][J] * VR[J];
          } // 230
        }

        W = (B[I][I]).abs() + (B[I + 1][I]).abs();
        if (W > SMLNUM) {
          if (W < ONE) {
            W1 = (XR).abs() + (XI).abs();
            if (W1 > W * BIGNUM) {
              REC = ONE / W1;
              dscal(N, REC, VR, 1);
              dscal(N, REC, VI, 1);
              XR = VR[I];
              XI = VI[I];
              SCALE = SCALE * REC;
              VMAX = VMAX * REC;
            }
          }

          // Divide by diagonal element of B.

          dladiv(XR, XI, B[I][I], B[I + 1][I], VR[I], VI[I]);
          VMAX = max((VR[I]).abs() + (VI[I]).abs(), VMAX);
          VCRIT = BIGNUM / VMAX;
        } else {
          for (J = 1; J <= N; J++) {
            // 240
            VR[J] = ZERO;
            VI[J] = ZERO;
          } // 240
          VR[I] = ONE;
          VI[I] = ONE;
          SCALE = ZERO;
          VMAX = ONE;
          VCRIT = BIGNUM;
        }
      } // 250

      // Test for sufficient growth in the norm of (VR,VI).

      VNORM = dasum(N, VR, 1) + dasum(N, VI, 1);
      if (VNORM >= GROWTO * SCALE) {
        eigenvectorFound = true;
        break;
      }

      // Choose a new orthogonal starting vector and try again.

      Y = EPS3 / (ROOTN + ONE);
      VR[1] = EPS3;
      VI[1] = ZERO;

      for (I = 2; I <= N; I++) {
        // 260
        VR[I] = Y;
        VI[I] = ZERO;
      } // 260
      VR[N - ITS + 1] = VR[N - ITS + 1] - EPS3 * ROOTN;
    } // 270

    // Failure to find eigenvector in N iterations
    if (!eigenvectorFound) {
      INFO.value = 1;
    } // 280

    // Normalize eigenvector.

    VNORM = ZERO;
    for (I = 1; I <= N; I++) {
      // 290
      VNORM = max(VNORM, (VR[I]).abs() + (VI[I]).abs());
    } // 290
    dscal(N, ONE / VNORM, VR, 1);
    dscal(N, ONE / VNORM, VI, 1);
  }
}
