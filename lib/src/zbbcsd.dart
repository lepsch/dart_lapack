import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlartgp.dart';
import 'package:lapack/src/dlartgs.dart';
import 'package:lapack/src/dlas2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlasr.dart';

void zbbcsd(
  final String JOBU1,
  final String JOBU2,
  final String JOBV1T,
  final String JOBV2T,
  final String TRANS,
  final int M,
  final int P,
  final int Q,
  final Array<double> THETA_,
  final Array<double> PHI_,
  final Matrix<Complex> U1_,
  final int LDU1,
  final Matrix<Complex> U2_,
  final int LDU2,
  final Matrix<Complex> V1T_,
  final int LDV1T,
  final Matrix<Complex> V2T_,
  final int LDV2T,
  final Array<double> B11D_,
  final Array<double> B11E_,
  final Array<double> B12D_,
  final Array<double> B12E_,
  final Array<double> B21D_,
  final Array<double> B21E_,
  final Array<double> B22D_,
  final Array<double> B22E_,
  final Array<double> RWORK_,
  final int LRWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final U2 = U2_.having(ld: LDU2);
  final U1 = U1_.having(ld: LDU1);
  final V1T = V1T_.having(ld: LDV1T);
  final V2T = V2T_.having(ld: LDV2T);
  final THETA = THETA_.having();
  final PHI = PHI_.having();
  final B11D = B11D_.having();
  final B11E = B11E_.having();
  final B12D = B12D_.having();
  final B12E = B12E_.having();
  final B21D = B21D_.having();
  final B21E = B21E_.having();
  final B22D = B22D_.having();
  final B22E = B22E_.having();
  final RWORK = RWORK_.having();

  const MAXITR = 6;
  const HUNDRED = 100.0, MEIGHTH = -0.125, ONE = 1.0, TEN = 10.0, ZERO = 0.0;
  const NEGONECOMPLEX = Complex(-1.0, 0.0);
  const PIOVER2 = 1.57079632679489661923132169163975144210;
  bool COLMAJOR,
      LQUERY,
      RESTART11,
      RESTART12,
      RESTART21,
      RESTART22,
      WANTU1,
      WANTU2,
      WANTV1T,
      WANTV2T;
  int I,
      IMIN,
      IMAX,
      ITER,
      IU1CS = 0,
      IU1SN = 0,
      IU2CS = 0,
      IU2SN = 0,
      IV1TCS = 0,
      IV1TSN = 0,
      IV2TCS = 0,
      IV2TSN = 0,
      J,
      LRWORKMIN,
      LRWORKOPT,
      MAXIT,
      MINI;
  double B11BULGE,
      B12BULGE,
      B21BULGE,
      B22BULGE,
      EPS,
      MU,
      NU,
      TEMP,
      THETAMAX,
      THETAMIN,
      THRESH,
      TOL,
      TOLMUL,
      UNFL,
      X1,
      X2,
      Y1,
      Y2;
  final SIGMA11 = Box(0.0), SIGMA21 = Box(0.0), DUMMY = Box(0.0), R = Box(0.0);

  // Test input arguments

  INFO.value = 0;
  LQUERY = LRWORK == -1;
  WANTU1 = lsame(JOBU1, 'Y');
  WANTU2 = lsame(JOBU2, 'Y');
  WANTV1T = lsame(JOBV1T, 'Y');
  WANTV2T = lsame(JOBV2T, 'Y');
  COLMAJOR = !lsame(TRANS, 'T');

  if (M < 0) {
    INFO.value = -6;
  } else if (P < 0 || P > M) {
    INFO.value = -7;
  } else if (Q < 0 || Q > M) {
    INFO.value = -8;
  } else if (Q > P || Q > M - P || Q > M - Q) {
    INFO.value = -8;
  } else if (WANTU1 && LDU1 < P) {
    INFO.value = -12;
  } else if (WANTU2 && LDU2 < M - P) {
    INFO.value = -14;
  } else if (WANTV1T && LDV1T < Q) {
    INFO.value = -16;
  } else if (WANTV2T && LDV2T < M - Q) {
    INFO.value = -18;
  }

  // Quick return if Q = 0

  if (INFO.value == 0 && Q == 0) {
    LRWORKMIN = 1;
    RWORK[1] = LRWORKMIN.toDouble();
    return;
  }

  // Compute workspace

  if (INFO.value == 0) {
    IU1CS = 1;
    IU1SN = IU1CS + Q;
    IU2CS = IU1SN + Q;
    IU2SN = IU2CS + Q;
    IV1TCS = IU2SN + Q;
    IV1TSN = IV1TCS + Q;
    IV2TCS = IV1TSN + Q;
    IV2TSN = IV2TCS + Q;
    LRWORKOPT = IV2TSN + Q - 1;
    LRWORKMIN = LRWORKOPT;
    RWORK[1] = LRWORKOPT.toDouble();
    if (LRWORK < LRWORKMIN && !LQUERY) {
      INFO.value = -28;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZBBCSD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Get machine constants

  EPS = dlamch('Epsilon');
  UNFL = dlamch('Safe minimum');
  TOLMUL = max(TEN, min(HUNDRED, pow(EPS, MEIGHTH))).toDouble();
  TOL = TOLMUL * EPS;
  THRESH = max(TOL, MAXITR * Q * Q * UNFL);

  // Test for negligible sines or cosines

  for (I = 1; I <= Q; I++) {
    if (THETA[I] < THRESH) {
      THETA[I] = ZERO;
    } else if (THETA[I] > PIOVER2 - THRESH) {
      THETA[I] = PIOVER2;
    }
  }
  for (I = 1; I <= Q - 1; I++) {
    if (PHI[I] < THRESH) {
      PHI[I] = ZERO;
    } else if (PHI[I] > PIOVER2 - THRESH) {
      PHI[I] = PIOVER2;
    }
  }

  // Initial deflation

  IMAX = Q;
  while (IMAX > 1) {
    if (PHI[IMAX - 1] != ZERO) {
      break;
    }
    IMAX--;
  }
  IMIN = IMAX - 1;
  if (IMIN > 1) {
    while (PHI[IMIN - 1] != ZERO) {
      IMIN--;
      if (IMIN <= 1) break;
    }
  }

  // Initialize iteration counter

  MAXIT = MAXITR * Q * Q;
  ITER = 0;

  // Begin main iteration loop

  while (IMAX > 1) {
    // Compute the matrix entries

    B11D[IMIN] = cos(THETA[IMIN]);
    B21D[IMIN] = -sin(THETA[IMIN]);
    for (I = IMIN; I <= IMAX - 1; I++) {
      B11E[I] = -sin(THETA[I]) * sin(PHI[I]);
      B11D[I + 1] = cos(THETA[I + 1]) * cos(PHI[I]);
      B12D[I] = sin(THETA[I]) * cos(PHI[I]);
      B12E[I] = cos(THETA[I + 1]) * sin(PHI[I]);
      B21E[I] = -cos(THETA[I]) * sin(PHI[I]);
      B21D[I + 1] = -sin(THETA[I + 1]) * cos(PHI[I]);
      B22D[I] = cos(THETA[I]) * cos(PHI[I]);
      B22E[I] = -sin(THETA[I + 1]) * sin(PHI[I]);
    }
    B12D[IMAX] = sin(THETA[IMAX]);
    B22D[IMAX] = cos(THETA[IMAX]);

    // Abort if not converging; otherwise, increment ITER

    if (ITER > MAXIT) {
      INFO.value = 0;
      for (I = 1; I <= Q; I++) {
        if (PHI[I] != ZERO) INFO.value++;
      }
      return;
    }

    ITER += IMAX - IMIN;

    // Compute shifts

    THETAMAX = THETA[IMIN];
    THETAMIN = THETA[IMIN];
    for (I = IMIN + 1; I <= IMAX; I++) {
      if (THETA[I] > THETAMAX) THETAMAX = THETA[I];
      if (THETA[I] < THETAMIN) THETAMIN = THETA[I];
    }

    if (THETAMAX > PIOVER2 - THRESH) {
      // Zero on diagonals of B11 and B22; induce deflation with a
      // zero shift

      MU = ZERO;
      NU = ONE;
    } else if (THETAMIN < THRESH) {
      // Zero on diagonals of B12 and B22; induce deflation with a
      // zero shift

      MU = ONE;
      NU = ZERO;
    } else {
      // Compute shifts for B11 and B21 and use the lesser

      dlas2(B11D[IMAX - 1], B11E[IMAX - 1], B11D[IMAX], SIGMA11, DUMMY);
      dlas2(B21D[IMAX - 1], B21E[IMAX - 1], B21D[IMAX], SIGMA21, DUMMY);

      if (SIGMA11.value <= SIGMA21.value) {
        MU = SIGMA11.value;
        NU = sqrt(ONE - pow(MU, 2));
        if (MU < THRESH) {
          MU = ZERO;
          NU = ONE;
        }
      } else {
        NU = SIGMA21.value;
        MU = sqrt(1.0 - pow(NU, 2));
        if (NU < THRESH) {
          MU = ONE;
          NU = ZERO;
        }
      }
    }

    // Rotate to produce bulges in B11 and B21

    if (MU <= NU) {
      dlartgs(B11D[IMIN], B11E[IMIN], MU, RWORK.box(IV1TCS + IMIN - 1),
          RWORK.box(IV1TSN + IMIN - 1));
    } else {
      dlartgs(B21D[IMIN], B21E[IMIN], NU, RWORK.box(IV1TCS + IMIN - 1),
          RWORK.box(IV1TSN + IMIN - 1));
    }

    TEMP = RWORK[IV1TCS + IMIN - 1] * B11D[IMIN] +
        RWORK[IV1TSN + IMIN - 1] * B11E[IMIN];
    B11E[IMIN] = RWORK[IV1TCS + IMIN - 1] * B11E[IMIN] -
        RWORK[IV1TSN + IMIN - 1] * B11D[IMIN];
    B11D[IMIN] = TEMP;
    B11BULGE = RWORK[IV1TSN + IMIN - 1] * B11D[IMIN + 1];
    B11D[IMIN + 1] = RWORK[IV1TCS + IMIN - 1] * B11D[IMIN + 1];
    TEMP = RWORK[IV1TCS + IMIN - 1] * B21D[IMIN] +
        RWORK[IV1TSN + IMIN - 1] * B21E[IMIN];
    B21E[IMIN] = RWORK[IV1TCS + IMIN - 1] * B21E[IMIN] -
        RWORK[IV1TSN + IMIN - 1] * B21D[IMIN];
    B21D[IMIN] = TEMP;
    B21BULGE = RWORK[IV1TSN + IMIN - 1] * B21D[IMIN + 1];
    B21D[IMIN + 1] = RWORK[IV1TCS + IMIN - 1] * B21D[IMIN + 1];

    // Compute THETA[IMIN]

    THETA[IMIN] = atan2(sqrt(pow(B21D[IMIN], 2) + pow(B21BULGE, 2)),
        sqrt(pow(B11D[IMIN], 2) + pow(B11BULGE, 2)));

    // Chase the bulges in B11[IMIN+1,IMIN] and B21[IMIN+1,IMIN]

    if (pow(B11D[IMIN], 2) + pow(B11BULGE, 2) > pow(THRESH, 2)) {
      dlartgp(B11BULGE, B11D[IMIN], RWORK.box(IU1SN + IMIN - 1),
          RWORK.box(IU1CS + IMIN - 1), R);
    } else if (MU <= NU) {
      dlartgs(B11E[IMIN], B11D[IMIN + 1], MU, RWORK.box(IU1CS + IMIN - 1),
          RWORK.box(IU1SN + IMIN - 1));
    } else {
      dlartgs(B12D[IMIN], B12E[IMIN], NU, RWORK.box(IU1CS + IMIN - 1),
          RWORK.box(IU1SN + IMIN - 1));
    }
    if (pow(B21D[IMIN], 2) + pow(B21BULGE, 2) > pow(THRESH, 2)) {
      dlartgp(B21BULGE, B21D[IMIN], RWORK.box(IU2SN + IMIN - 1),
          RWORK.box(IU2CS + IMIN - 1), R);
    } else if (NU < MU) {
      dlartgs(B21E[IMIN], B21D[IMIN + 1], NU, RWORK.box(IU2CS + IMIN - 1),
          RWORK.box(IU2SN + IMIN - 1));
    } else {
      dlartgs(B22D[IMIN], B22E[IMIN], MU, RWORK.box(IU2CS + IMIN - 1),
          RWORK.box(IU2SN + IMIN - 1));
    }
    RWORK[IU2CS + IMIN - 1] = -RWORK[IU2CS + IMIN - 1];
    RWORK[IU2SN + IMIN - 1] = -RWORK[IU2SN + IMIN - 1];

    TEMP = RWORK[IU1CS + IMIN - 1] * B11E[IMIN] +
        RWORK[IU1SN + IMIN - 1] * B11D[IMIN + 1];
    B11D[IMIN + 1] = RWORK[IU1CS + IMIN - 1] * B11D[IMIN + 1] -
        RWORK[IU1SN + IMIN - 1] * B11E[IMIN];
    B11E[IMIN] = TEMP;
    if (IMAX > IMIN + 1) {
      B11BULGE = RWORK[IU1SN + IMIN - 1] * B11E[IMIN + 1];
      B11E[IMIN + 1] = RWORK[IU1CS + IMIN - 1] * B11E[IMIN + 1];
    }
    TEMP = RWORK[IU1CS + IMIN - 1] * B12D[IMIN] +
        RWORK[IU1SN + IMIN - 1] * B12E[IMIN];
    B12E[IMIN] = RWORK[IU1CS + IMIN - 1] * B12E[IMIN] -
        RWORK[IU1SN + IMIN - 1] * B12D[IMIN];
    B12D[IMIN] = TEMP;
    B12BULGE = RWORK[IU1SN + IMIN - 1] * B12D[IMIN + 1];
    B12D[IMIN + 1] = RWORK[IU1CS + IMIN - 1] * B12D[IMIN + 1];
    TEMP = RWORK[IU2CS + IMIN - 1] * B21E[IMIN] +
        RWORK[IU2SN + IMIN - 1] * B21D[IMIN + 1];
    B21D[IMIN + 1] = RWORK[IU2CS + IMIN - 1] * B21D[IMIN + 1] -
        RWORK[IU2SN + IMIN - 1] * B21E[IMIN];
    B21E[IMIN] = TEMP;
    if (IMAX > IMIN + 1) {
      B21BULGE = RWORK[IU2SN + IMIN - 1] * B21E[IMIN + 1];
      B21E[IMIN + 1] = RWORK[IU2CS + IMIN - 1] * B21E[IMIN + 1];
    }
    TEMP = RWORK[IU2CS + IMIN - 1] * B22D[IMIN] +
        RWORK[IU2SN + IMIN - 1] * B22E[IMIN];
    B22E[IMIN] = RWORK[IU2CS + IMIN - 1] * B22E[IMIN] -
        RWORK[IU2SN + IMIN - 1] * B22D[IMIN];
    B22D[IMIN] = TEMP;
    B22BULGE = RWORK[IU2SN + IMIN - 1] * B22D[IMIN + 1];
    B22D[IMIN + 1] = RWORK[IU2CS + IMIN - 1] * B22D[IMIN + 1];

    // Inner loop: chase bulges from B11[IMIN,IMIN+2],
    // B12[IMIN,IMIN+1], B21[IMIN,IMIN+2], and B22[IMIN,IMIN+1] to
    // bottom-right

    for (I = IMIN + 1; I <= IMAX - 1; I++) {
      // Compute PHI[I-1]

      X1 = sin(THETA[I - 1]) * B11E[I - 1] + cos(THETA[I - 1]) * B21E[I - 1];
      X2 = sin(THETA[I - 1]) * B11BULGE + cos(THETA[I - 1]) * B21BULGE;
      Y1 = sin(THETA[I - 1]) * B12D[I - 1] + cos(THETA[I - 1]) * B22D[I - 1];
      Y2 = sin(THETA[I - 1]) * B12BULGE + cos(THETA[I - 1]) * B22BULGE;

      PHI[I - 1] =
          atan2(sqrt(pow(X1, 2) + pow(X2, 2)), sqrt(pow(Y1, 2) + pow(Y2, 2)));

      // Determine if there are bulges to chase or if a new direct
      // summand has been reached

      RESTART11 = pow(B11E[I - 1], 2) + pow(B11BULGE, 2) <= pow(THRESH, 2);
      RESTART21 = pow(B21E[I - 1], 2) + pow(B21BULGE, 2) <= pow(THRESH, 2);
      RESTART12 = pow(B12D[I - 1], 2) + pow(B12BULGE, 2) <= pow(THRESH, 2);
      RESTART22 = pow(B22D[I - 1], 2) + pow(B22BULGE, 2) <= pow(THRESH, 2);

      // If possible, chase bulges from B11[I-1,I+1], B12[I-1,I],
      // B21[I-1,I+1], and B22[I-1,I]. If necessary, restart bulge-
      // chasing by applying the original shift again.

      if (!RESTART11 && !RESTART21) {
        dlartgp(
            X2, X1, RWORK.box(IV1TSN + I - 1), RWORK.box(IV1TCS + I - 1), R);
      } else if (!RESTART11 && RESTART21) {
        dlartgp(B11BULGE, B11E[I - 1], RWORK.box(IV1TSN + I - 1),
            RWORK.box(IV1TCS + I - 1), R);
      } else if (RESTART11 && !RESTART21) {
        dlartgp(B21BULGE, B21E[I - 1], RWORK.box(IV1TSN + I - 1),
            RWORK.box(IV1TCS + I - 1), R);
      } else if (MU <= NU) {
        dlartgs(B11D[I], B11E[I], MU, RWORK.box(IV1TCS + I - 1),
            RWORK.box(IV1TSN + I - 1));
      } else {
        dlartgs(B21D[I], B21E[I], NU, RWORK.box(IV1TCS + I - 1),
            RWORK.box(IV1TSN + I - 1));
      }
      RWORK[IV1TCS + I - 1] = -RWORK[IV1TCS + I - 1];
      RWORK[IV1TSN + I - 1] = -RWORK[IV1TSN + I - 1];
      if (!RESTART12 && !RESTART22) {
        dlartgp(Y2, Y1, RWORK.box(IV2TSN + I - 1 - 1),
            RWORK.box(IV2TCS + I - 1 - 1), R);
      } else if (!RESTART12 && RESTART22) {
        dlartgp(B12BULGE, B12D[I - 1], RWORK.box(IV2TSN + I - 1 - 1),
            RWORK.box(IV2TCS + I - 1 - 1), R);
      } else if (RESTART12 && !RESTART22) {
        dlartgp(B22BULGE, B22D[I - 1], RWORK.box(IV2TSN + I - 1 - 1),
            RWORK.box(IV2TCS + I - 1 - 1), R);
      } else if (NU < MU) {
        dlartgs(B12E[I - 1], B12D[I], NU, RWORK.box(IV2TCS + I - 1 - 1),
            RWORK.box(IV2TSN + I - 1 - 1));
      } else {
        dlartgs(B22E[I - 1], B22D[I], MU, RWORK.box(IV2TCS + I - 1 - 1),
            RWORK.box(IV2TSN + I - 1 - 1));
      }

      TEMP = RWORK[IV1TCS + I - 1] * B11D[I] + RWORK[IV1TSN + I - 1] * B11E[I];
      B11E[I] =
          RWORK[IV1TCS + I - 1] * B11E[I] - RWORK[IV1TSN + I - 1] * B11D[I];
      B11D[I] = TEMP;
      B11BULGE = RWORK[IV1TSN + I - 1] * B11D[I + 1];
      B11D[I + 1] = RWORK[IV1TCS + I - 1] * B11D[I + 1];
      TEMP = RWORK[IV1TCS + I - 1] * B21D[I] + RWORK[IV1TSN + I - 1] * B21E[I];
      B21E[I] =
          RWORK[IV1TCS + I - 1] * B21E[I] - RWORK[IV1TSN + I - 1] * B21D[I];
      B21D[I] = TEMP;
      B21BULGE = RWORK[IV1TSN + I - 1] * B21D[I + 1];
      B21D[I + 1] = RWORK[IV1TCS + I - 1] * B21D[I + 1];
      TEMP = RWORK[IV2TCS + I - 1 - 1] * B12E[I - 1] +
          RWORK[IV2TSN + I - 1 - 1] * B12D[I];
      B12D[I] = RWORK[IV2TCS + I - 1 - 1] * B12D[I] -
          RWORK[IV2TSN + I - 1 - 1] * B12E[I - 1];
      B12E[I - 1] = TEMP;
      B12BULGE = RWORK[IV2TSN + I - 1 - 1] * B12E[I];
      B12E[I] = RWORK[IV2TCS + I - 1 - 1] * B12E[I];
      TEMP = RWORK[IV2TCS + I - 1 - 1] * B22E[I - 1] +
          RWORK[IV2TSN + I - 1 - 1] * B22D[I];
      B22D[I] = RWORK[IV2TCS + I - 1 - 1] * B22D[I] -
          RWORK[IV2TSN + I - 1 - 1] * B22E[I - 1];
      B22E[I - 1] = TEMP;
      B22BULGE = RWORK[IV2TSN + I - 1 - 1] * B22E[I];
      B22E[I] = RWORK[IV2TCS + I - 1 - 1] * B22E[I];

      // Compute THETA[I]

      X1 = cos(PHI[I - 1]) * B11D[I] + sin(PHI[I - 1]) * B12E[I - 1];
      X2 = cos(PHI[I - 1]) * B11BULGE + sin(PHI[I - 1]) * B12BULGE;
      Y1 = cos(PHI[I - 1]) * B21D[I] + sin(PHI[I - 1]) * B22E[I - 1];
      Y2 = cos(PHI[I - 1]) * B21BULGE + sin(PHI[I - 1]) * B22BULGE;

      THETA[I] =
          atan2(sqrt(pow(Y1, 2) + pow(Y2, 2)), sqrt(pow(X1, 2) + pow(X2, 2)));

      // Determine if there are bulges to chase or if a new direct
      // summand has been reached

      RESTART11 = pow(B11D[I], 2) + pow(B11BULGE, 2) <= pow(THRESH, 2);
      RESTART12 = pow(B12E[I - 1], 2) + pow(B12BULGE, 2) <= pow(THRESH, 2);
      RESTART21 = pow(B21D[I], 2) + pow(B21BULGE, 2) <= pow(THRESH, 2);
      RESTART22 = pow(B22E[I - 1], 2) + pow(B22BULGE, 2) <= pow(THRESH, 2);

      // If possible, chase bulges from B11[I+1,I], B12[I+1,I-1],
      // B21[I+1,I], and B22[I+1,I-1]. If necessary, restart bulge-
      // chasing by applying the original shift again.

      if (!RESTART11 && !RESTART12) {
        dlartgp(X2, X1, RWORK.box(IU1SN + I - 1), RWORK.box(IU1CS + I - 1), R);
      } else if (!RESTART11 && RESTART12) {
        dlartgp(B11BULGE, B11D[I], RWORK.box(IU1SN + I - 1),
            RWORK.box(IU1CS + I - 1), R);
      } else if (RESTART11 && !RESTART12) {
        dlartgp(B12BULGE, B12E[I - 1], RWORK.box(IU1SN + I - 1),
            RWORK.box(IU1CS + I - 1), R);
      } else if (MU <= NU) {
        dlartgs(B11E[I], B11D[I + 1], MU, RWORK.box(IU1CS + I - 1),
            RWORK.box(IU1SN + I - 1));
      } else {
        dlartgs(B12D[I], B12E[I], NU, RWORK.box(IU1CS + I - 1),
            RWORK.box(IU1SN + I - 1));
      }
      if (!RESTART21 && !RESTART22) {
        dlartgp(Y2, Y1, RWORK.box(IU2SN + I - 1), RWORK.box(IU2CS + I - 1), R);
      } else if (!RESTART21 && RESTART22) {
        dlartgp(B21BULGE, B21D[I], RWORK.box(IU2SN + I - 1),
            RWORK.box(IU2CS + I - 1), R);
      } else if (RESTART21 && !RESTART22) {
        dlartgp(B22BULGE, B22E[I - 1], RWORK.box(IU2SN + I - 1),
            RWORK.box(IU2CS + I - 1), R);
      } else if (NU < MU) {
        dlartgs(B21E[I], B21E[I + 1], NU, RWORK.box(IU2CS + I - 1),
            RWORK.box(IU2SN + I - 1));
      } else {
        dlartgs(B22D[I], B22E[I], MU, RWORK.box(IU2CS + I - 1),
            RWORK.box(IU2SN + I - 1));
      }
      RWORK[IU2CS + I - 1] = -RWORK[IU2CS + I - 1];
      RWORK[IU2SN + I - 1] = -RWORK[IU2SN + I - 1];

      TEMP =
          RWORK[IU1CS + I - 1] * B11E[I] + RWORK[IU1SN + I - 1] * B11D[I + 1];
      B11D[I + 1] =
          RWORK[IU1CS + I - 1] * B11D[I + 1] - RWORK[IU1SN + I - 1] * B11E[I];
      B11E[I] = TEMP;
      if (I < IMAX - 1) {
        B11BULGE = RWORK[IU1SN + I - 1] * B11E[I + 1];
        B11E[I + 1] = RWORK[IU1CS + I - 1] * B11E[I + 1];
      }
      TEMP =
          RWORK[IU2CS + I - 1] * B21E[I] + RWORK[IU2SN + I - 1] * B21D[I + 1];
      B21D[I + 1] =
          RWORK[IU2CS + I - 1] * B21D[I + 1] - RWORK[IU2SN + I - 1] * B21E[I];
      B21E[I] = TEMP;
      if (I < IMAX - 1) {
        B21BULGE = RWORK[IU2SN + I - 1] * B21E[I + 1];
        B21E[I + 1] = RWORK[IU2CS + I - 1] * B21E[I + 1];
      }
      TEMP = RWORK[IU1CS + I - 1] * B12D[I] + RWORK[IU1SN + I - 1] * B12E[I];
      B12E[I] = RWORK[IU1CS + I - 1] * B12E[I] - RWORK[IU1SN + I - 1] * B12D[I];
      B12D[I] = TEMP;
      B12BULGE = RWORK[IU1SN + I - 1] * B12D[I + 1];
      B12D[I + 1] = RWORK[IU1CS + I - 1] * B12D[I + 1];
      TEMP = RWORK[IU2CS + I - 1] * B22D[I] + RWORK[IU2SN + I - 1] * B22E[I];
      B22E[I] = RWORK[IU2CS + I - 1] * B22E[I] - RWORK[IU2SN + I - 1] * B22D[I];
      B22D[I] = TEMP;
      B22BULGE = RWORK[IU2SN + I - 1] * B22D[I + 1];
      B22D[I + 1] = RWORK[IU2CS + I - 1] * B22D[I + 1];
    }

    // Compute PHI[IMAX-1]

    X1 = sin(THETA[IMAX - 1]) * B11E[IMAX - 1] +
        cos(THETA[IMAX - 1]) * B21E[IMAX - 1];
    Y1 = sin(THETA[IMAX - 1]) * B12D[IMAX - 1] +
        cos(THETA[IMAX - 1]) * B22D[IMAX - 1];
    Y2 = sin(THETA[IMAX - 1]) * B12BULGE + cos(THETA[IMAX - 1]) * B22BULGE;

    PHI[IMAX - 1] = atan2(X1.abs(), sqrt(pow(Y1, 2) + pow(Y2, 2)));

    // Chase bulges from B12[IMAX-1,IMAX] and B22[IMAX-1,IMAX]

    RESTART12 = pow(B12D[IMAX - 1], 2) + pow(B12BULGE, 2) <= pow(THRESH, 2);
    RESTART22 = pow(B22D[IMAX - 1], 2) + pow(B22BULGE, 2) <= pow(THRESH, 2);

    if (!RESTART12 && !RESTART22) {
      dlartgp(Y2, Y1, RWORK.box(IV2TSN + IMAX - 1 - 1),
          RWORK.box(IV2TCS + IMAX - 1 - 1), R);
    } else if (!RESTART12 && RESTART22) {
      dlartgp(B12BULGE, B12D[IMAX - 1], RWORK.box(IV2TSN + IMAX - 1 - 1),
          RWORK.box(IV2TCS + IMAX - 1 - 1), R);
    } else if (RESTART12 && !RESTART22) {
      dlartgp(B22BULGE, B22D[IMAX - 1], RWORK.box(IV2TSN + IMAX - 1 - 1),
          RWORK.box(IV2TCS + IMAX - 1 - 1), R);
    } else if (NU < MU) {
      dlartgs(B12E[IMAX - 1], B12D[IMAX], NU, RWORK.box(IV2TCS + IMAX - 1 - 1),
          RWORK.box(IV2TSN + IMAX - 1 - 1));
    } else {
      dlartgs(B22E[IMAX - 1], B22D[IMAX], MU, RWORK.box(IV2TCS + IMAX - 1 - 1),
          RWORK.box(IV2TSN + IMAX - 1 - 1));
    }

    TEMP = RWORK[IV2TCS + IMAX - 1 - 1] * B12E[IMAX - 1] +
        RWORK[IV2TSN + IMAX - 1 - 1] * B12D[IMAX];
    B12D[IMAX] = RWORK[IV2TCS + IMAX - 1 - 1] * B12D[IMAX] -
        RWORK[IV2TSN + IMAX - 1 - 1] * B12E[IMAX - 1];
    B12E[IMAX - 1] = TEMP;
    TEMP = RWORK[IV2TCS + IMAX - 1 - 1] * B22E[IMAX - 1] +
        RWORK[IV2TSN + IMAX - 1 - 1] * B22D[IMAX];
    B22D[IMAX] = RWORK[IV2TCS + IMAX - 1 - 1] * B22D[IMAX] -
        RWORK[IV2TSN + IMAX - 1 - 1] * B22E[IMAX - 1];
    B22E[IMAX - 1] = TEMP;

    // Update singular vectors

    if (WANTU1) {
      if (COLMAJOR) {
        zlasr('R', 'V', 'F', P, IMAX - IMIN + 1, RWORK(IU1CS + IMIN - 1),
            RWORK(IU1SN + IMIN - 1), U1(1, IMIN), LDU1);
      } else {
        zlasr('L', 'V', 'F', IMAX - IMIN + 1, P, RWORK(IU1CS + IMIN - 1),
            RWORK(IU1SN + IMIN - 1), U1(IMIN, 1), LDU1);
      }
    }
    if (WANTU2) {
      if (COLMAJOR) {
        zlasr('R', 'V', 'F', M - P, IMAX - IMIN + 1, RWORK(IU2CS + IMIN - 1),
            RWORK(IU2SN + IMIN - 1), U2(1, IMIN), LDU2);
      } else {
        zlasr('L', 'V', 'F', IMAX - IMIN + 1, M - P, RWORK(IU2CS + IMIN - 1),
            RWORK(IU2SN + IMIN - 1), U2(IMIN, 1), LDU2);
      }
    }
    if (WANTV1T) {
      if (COLMAJOR) {
        zlasr('L', 'V', 'F', IMAX - IMIN + 1, Q, RWORK(IV1TCS + IMIN - 1),
            RWORK(IV1TSN + IMIN - 1), V1T(IMIN, 1), LDV1T);
      } else {
        zlasr('R', 'V', 'F', Q, IMAX - IMIN + 1, RWORK(IV1TCS + IMIN - 1),
            RWORK(IV1TSN + IMIN - 1), V1T(1, IMIN), LDV1T);
      }
    }
    if (WANTV2T) {
      if (COLMAJOR) {
        zlasr('L', 'V', 'F', IMAX - IMIN + 1, M - Q, RWORK(IV2TCS + IMIN - 1),
            RWORK(IV2TSN + IMIN - 1), V2T(IMIN, 1), LDV2T);
      } else {
        zlasr('R', 'V', 'F', M - Q, IMAX - IMIN + 1, RWORK(IV2TCS + IMIN - 1),
            RWORK(IV2TSN + IMIN - 1), V2T(1, IMIN), LDV2T);
      }
    }

    // Fix signs on B11[IMAX-1,IMAX] and B21[IMAX-1,IMAX]

    if (B11E[IMAX - 1] + B21E[IMAX - 1] > 0) {
      B11D[IMAX] = -B11D[IMAX];
      B21D[IMAX] = -B21D[IMAX];
      if (WANTV1T) {
        if (COLMAJOR) {
          zscal(Q, NEGONECOMPLEX, V1T(IMAX, 1).asArray(), LDV1T);
        } else {
          zscal(Q, NEGONECOMPLEX, V1T(1, IMAX).asArray(), 1);
        }
      }
    }

    // Compute THETA[IMAX]

    X1 = cos(PHI[IMAX - 1]) * B11D[IMAX] + sin(PHI[IMAX - 1]) * B12E[IMAX - 1];
    Y1 = cos(PHI[IMAX - 1]) * B21D[IMAX] + sin(PHI[IMAX - 1]) * B22E[IMAX - 1];

    THETA[IMAX] = atan2(Y1.abs(), X1.abs());

    // Fix signs on B11[IMAX,IMAX], B12[IMAX,IMAX-1], B21[IMAX,IMAX],
    // and B22[IMAX,IMAX-1]

    if (B11D[IMAX] + B12E[IMAX - 1] < 0) {
      B12D[IMAX] = -B12D[IMAX];
      if (WANTU1) {
        if (COLMAJOR) {
          zscal(P, NEGONECOMPLEX, U1(1, IMAX).asArray(), 1);
        } else {
          zscal(P, NEGONECOMPLEX, U1(IMAX, 1).asArray(), LDU1);
        }
      }
    }
    if (B21D[IMAX] + B22E[IMAX - 1] > 0) {
      B22D[IMAX] = -B22D[IMAX];
      if (WANTU2) {
        if (COLMAJOR) {
          zscal(M - P, NEGONECOMPLEX, U2(1, IMAX).asArray(), 1);
        } else {
          zscal(M - P, NEGONECOMPLEX, U2(IMAX, 1).asArray(), LDU2);
        }
      }
    }

    // Fix signs on B12[IMAX,IMAX] and B22[IMAX,IMAX]

    if (B12D[IMAX] + B22D[IMAX] < 0) {
      if (WANTV2T) {
        if (COLMAJOR) {
          zscal(M - Q, NEGONECOMPLEX, V2T(IMAX, 1).asArray(), LDV2T);
        } else {
          zscal(M - Q, NEGONECOMPLEX, V2T(1, IMAX).asArray(), 1);
        }
      }
    }

    // Test for negligible sines or cosines

    for (I = IMIN; I <= IMAX; I++) {
      if (THETA[I] < THRESH) {
        THETA[I] = ZERO;
      } else if (THETA[I] > PIOVER2 - THRESH) {
        THETA[I] = PIOVER2;
      }
    }
    for (I = IMIN; I <= IMAX - 1; I++) {
      if (PHI[I] < THRESH) {
        PHI[I] = ZERO;
      } else if (PHI[I] > PIOVER2 - THRESH) {
        PHI[I] = PIOVER2;
      }
    }

    // Deflate

    if (IMAX > 1) {
      while (PHI[IMAX - 1] == ZERO) {
        IMAX--;
        if (IMAX <= 1) break;
      }
    }
    if (IMIN > IMAX - 1) IMIN = IMAX - 1;
    if (IMIN > 1) {
      while (PHI[IMIN - 1] != ZERO) {
        IMIN--;
        if (IMIN <= 1) break;
      }
    }

    // Repeat main iteration loop
  }

  // Postprocessing: order THETA from least to greatest

  for (I = 1; I <= Q; I++) {
    MINI = I;
    THETAMIN = THETA[I];
    for (J = I + 1; J <= Q; J++) {
      if (THETA[J] < THETAMIN) {
        MINI = J;
        THETAMIN = THETA[J];
      }
    }

    if (MINI != I) {
      THETA[MINI] = THETA[I];
      THETA[I] = THETAMIN;
      if (COLMAJOR) {
        if (WANTU1) zswap(P, U1(1, I).asArray(), 1, U1(1, MINI).asArray(), 1);
        if (WANTU2) {
          zswap(M - P, U2(1, I).asArray(), 1, U2(1, MINI).asArray(), 1);
        }
        if (WANTV1T) {
          zswap(Q, V1T(I, 1).asArray(), LDV1T, V1T(MINI, 1).asArray(), LDV1T);
        }
        if (WANTV2T) {
          zswap(
              M - Q, V2T(I, 1).asArray(), LDV2T, V2T(MINI, 1).asArray(), LDV2T);
        }
      } else {
        if (WANTU1) {
          zswap(P, U1(I, 1).asArray(), LDU1, U1(MINI, 1).asArray(), LDU1);
        }
        if (WANTU2) {
          zswap(M - P, U2(I, 1).asArray(), LDU2, U2(MINI, 1).asArray(), LDU2);
        }
        if (WANTV1T) {
          zswap(Q, V1T(1, I).asArray(), 1, V1T(1, MINI).asArray(), 1);
        }
        if (WANTV2T) {
          zswap(M - Q, V2T(1, I).asArray(), 1, V2T(1, MINI).asArray(), 1);
        }
      }
    }
  }
}
