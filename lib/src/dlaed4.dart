import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaed5.dart';
import 'package:lapack/src/dlaed6.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaed4(
  final int N,
  final int I,
  final Array<double> D_,
  final Array<double> Z_,
  final Array<double> DELTA_,
  final double RHO,
  final Box<double> DLAM,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Z = Z_.having();
  final DELTA = DELTA_.having();
  const MAXIT = 30;
  const ZERO = 0.0,
      ONE = 1.0,
      TWO = 2.0,
      THREE = 3.0,
      FOUR = 4.0,
      EIGHT = 8.0,
      TEN = 10.0;
  bool ORGATI, SWTCH, SWTCH3;
  int II, IIM1, IIP1, IP1, ITER, J, NITER;
  double A,
      B,
      C,
      DEL,
      DLTLB,
      DLTUB,
      DPHI,
      DPSI,
      DW,
      EPS,
      ERRETM,
      MIDPT,
      PHI,
      PREW,
      PSI,
      RHOINV,
      TAU,
      TEMP,
      TEMP1,
      W;
  final ZZ = Array<double>(3);
  final ETA = Box(0.0);

  // Quick return for N=1 and 2.

  INFO.value = 0;
  if (N == 1) {
    // Presumably, I=1 upon entry

    DLAM.value = D[1] + RHO * Z[1] * Z[1];
    DELTA[1] = ONE;
    return;
  }
  if (N == 2) {
    dlaed5(I, D, Z, DELTA, RHO, DLAM);
    return;
  }

  // Compute machine epsilon

  EPS = dlamch('Epsilon');
  RHOINV = ONE / RHO;

  // The case I = N

  if (I == N) {
    // Initialize some basic variables

    II = N - 1;
    NITER = 1;

    // Calculate initial guess

    MIDPT = RHO / TWO;

    // If ||Z||_2 is not one, then TEMP should be set to
    // RHO * ||Z||_2^2 / TWO

    for (J = 1; J <= N; J++) {
      DELTA[J] = (D[J] - D[I]) - MIDPT;
    }

    PSI = ZERO;
    for (J = 1; J <= N - 2; J++) {
      PSI = PSI + Z[J] * Z[J] / DELTA[J];
    }

    C = RHOINV + PSI;
    W = C + Z[II] * Z[II] / DELTA[II] + Z[N] * Z[N] / DELTA[N];

    if (W <= ZERO) {
      TEMP = Z[N - 1] * Z[N - 1] / (D[N] - D[N - 1] + RHO) + Z[N] * Z[N] / RHO;
      if (C <= TEMP) {
        TAU = RHO;
      } else {
        DEL = D[N] - D[N - 1];
        A = -C * DEL + Z[N - 1] * Z[N - 1] + Z[N] * Z[N];
        B = Z[N] * Z[N] * DEL;
        if (A < ZERO) {
          TAU = TWO * B / (sqrt(A * A + FOUR * B * C) - A);
        } else {
          TAU = (A + sqrt(A * A + FOUR * B * C)) / (TWO * C);
        }
      }

      // It can be proved that
      // D[N]+RHO/2 <= LAMBDA(N) < D[N]+TAU <= D[N]+RHO

      DLTLB = MIDPT;
      DLTUB = RHO;
    } else {
      DEL = D[N] - D[N - 1];
      A = -C * DEL + Z[N - 1] * Z[N - 1] + Z[N] * Z[N];
      B = Z[N] * Z[N] * DEL;
      if (A < ZERO) {
        TAU = TWO * B / (sqrt(A * A + FOUR * B * C) - A);
      } else {
        TAU = (A + sqrt(A * A + FOUR * B * C)) / (TWO * C);
      }

      // It can be proved that
      // D[N] < D[N]+TAU < LAMBDA(N) < D[N]+RHO/2

      DLTLB = ZERO;
      DLTUB = MIDPT;
    }

    for (J = 1; J <= N; J++) {
      DELTA[J] = (D[J] - D[I]) - TAU;
    }

    // Evaluate PSI and the derivative DPSI

    DPSI = ZERO;
    PSI = ZERO;
    ERRETM = ZERO;
    for (J = 1; J <= II; J++) {
      TEMP = Z[J] / DELTA[J];
      PSI = PSI + Z[J] * TEMP;
      DPSI = DPSI + TEMP * TEMP;
      ERRETM = ERRETM + PSI;
    }
    ERRETM = (ERRETM).abs();

    // Evaluate PHI and the derivative DPHI

    TEMP = Z[N] / DELTA[N];
    PHI = Z[N] * TEMP;
    DPHI = TEMP * TEMP;
    ERRETM = EIGHT * (-PHI - PSI) +
        ERRETM -
        PHI +
        RHOINV +
        (TAU).abs() * (DPSI + DPHI);

    W = RHOINV + PHI + PSI;

    // Test for convergence

    if ((W).abs() <= EPS * ERRETM) {
      DLAM.value = D[I] + TAU;
      return;
    }

    if (W <= ZERO) {
      DLTLB = max(DLTLB, TAU);
    } else {
      DLTUB = min(DLTUB, TAU);
    }

    // Calculate the new step

    NITER = NITER + 1;
    C = W - DELTA[N - 1] * DPSI - DELTA[N] * DPHI;
    A = (DELTA[N - 1] + DELTA[N]) * W - DELTA[N - 1] * DELTA[N] * (DPSI + DPHI);
    B = DELTA[N - 1] * DELTA[N] * W;
    if (C < ZERO) C = (C).abs();
    if (C == ZERO) {
      // ETA.value = B/A
      // ETA.value = RHO - TAU
      // ETA.value = DLTUB - TAU

      // Update proposed by Li, Ren-Cang:
      ETA.value = -W / (DPSI + DPHI);
    } else if (A >= ZERO) {
      ETA.value = (A + sqrt((A * A - FOUR * B * C).abs())) / (TWO * C);
    } else {
      ETA.value = TWO * B / (A - sqrt((A * A - FOUR * B * C).abs()));
    }

    // Note, eta should be positive if w is negative, and
    // eta should be negative otherwise. However,
    // if for some reason caused by roundoff, eta*w > 0,
    // we simply use one Newton step instead. This way
    // will guarantee eta*w < 0.

    if (W * ETA.value > ZERO) ETA.value = -W / (DPSI + DPHI);
    TEMP = TAU + ETA.value;
    if (TEMP > DLTUB || TEMP < DLTLB) {
      if (W < ZERO) {
        ETA.value = (DLTUB - TAU) / TWO;
      } else {
        ETA.value = (DLTLB - TAU) / TWO;
      }
    }
    for (J = 1; J <= N; J++) {
      DELTA[J] = DELTA[J] - ETA.value;
    }

    TAU = TAU + ETA.value;

    // Evaluate PSI and the derivative DPSI

    DPSI = ZERO;
    PSI = ZERO;
    ERRETM = ZERO;
    for (J = 1; J <= II; J++) {
      TEMP = Z[J] / DELTA[J];
      PSI = PSI + Z[J] * TEMP;
      DPSI = DPSI + TEMP * TEMP;
      ERRETM = ERRETM + PSI;
    }
    ERRETM = (ERRETM).abs();

    // Evaluate PHI and the derivative DPHI

    TEMP = Z[N] / DELTA[N];
    PHI = Z[N] * TEMP;
    DPHI = TEMP * TEMP;
    ERRETM = EIGHT * (-PHI - PSI) +
        ERRETM -
        PHI +
        RHOINV +
        (TAU).abs() * (DPSI + DPHI);

    W = RHOINV + PHI + PSI;

    // Main loop to update the values of the array   DELTA

    ITER = NITER + 1;

    for (NITER = ITER; NITER <= MAXIT; NITER++) {
      // Test for convergence

      if ((W).abs() <= EPS * ERRETM) {
        DLAM.value = D[I] + TAU;
        return;
      }

      if (W <= ZERO) {
        DLTLB = max(DLTLB, TAU);
      } else {
        DLTUB = min(DLTUB, TAU);
      }

      // Calculate the new step

      C = W - DELTA[N - 1] * DPSI - DELTA[N] * DPHI;
      A = (DELTA[N - 1] + DELTA[N]) * W -
          DELTA[N - 1] * DELTA[N] * (DPSI + DPHI);
      B = DELTA[N - 1] * DELTA[N] * W;
      if (A >= ZERO) {
        ETA.value = (A + sqrt((A * A - FOUR * B * C).abs())) / (TWO * C);
      } else {
        ETA.value = TWO * B / (A - sqrt((A * A - FOUR * B * C).abs()));
      }

      // Note, eta should be positive if w is negative, and
      // eta should be negative otherwise. However,
      // if for some reason caused by roundoff, eta*w > 0,
      // we simply use one Newton step instead. This way
      // will guarantee eta*w < 0.

      if (W * ETA.value > ZERO) ETA.value = -W / (DPSI + DPHI);
      TEMP = TAU + ETA.value;
      if (TEMP > DLTUB || TEMP < DLTLB) {
        if (W < ZERO) {
          ETA.value = (DLTUB - TAU) / TWO;
        } else {
          ETA.value = (DLTLB - TAU) / TWO;
        }
      }
      for (J = 1; J <= N; J++) {
        DELTA[J] = DELTA[J] - ETA.value;
      }

      TAU = TAU + ETA.value;

      // Evaluate PSI and the derivative DPSI

      DPSI = ZERO;
      PSI = ZERO;
      ERRETM = ZERO;
      for (J = 1; J <= II; J++) {
        TEMP = Z[J] / DELTA[J];
        PSI = PSI + Z[J] * TEMP;
        DPSI = DPSI + TEMP * TEMP;
        ERRETM = ERRETM + PSI;
      }
      ERRETM = (ERRETM).abs();

      // Evaluate PHI and the derivative DPHI

      TEMP = Z[N] / DELTA[N];
      PHI = Z[N] * TEMP;
      DPHI = TEMP * TEMP;
      ERRETM = EIGHT * (-PHI - PSI) +
          ERRETM -
          PHI +
          RHOINV +
          (TAU).abs() * (DPSI + DPHI);

      W = RHOINV + PHI + PSI;
    }

    // Return with INFO.value = 1, NITER = MAXIT and not converged

    INFO.value = 1;
    DLAM.value = D[I] + TAU;
    return;

    // End for the case I = N
  } else {
    // The case for I < N

    NITER = 1;
    IP1 = I + 1;

    // Calculate initial guess

    DEL = D[IP1] - D[I];
    MIDPT = DEL / TWO;
    for (J = 1; J <= N; J++) {
      DELTA[J] = (D[J] - D[I]) - MIDPT;
    }

    PSI = ZERO;
    for (J = 1; J <= I - 1; J++) {
      PSI = PSI + Z[J] * Z[J] / DELTA[J];
    }

    PHI = ZERO;
    for (J = N; J >= I + 2; J--) {
      PHI = PHI + Z[J] * Z[J] / DELTA[J];
    }
    C = RHOINV + PSI + PHI;
    W = C + Z[I] * Z[I] / DELTA[I] + Z[IP1] * Z[IP1] / DELTA[IP1];

    if (W > ZERO) {
      // d(i)< the ith eigenvalue < (d(i)+d(i+1))/2

      // We choose d(i) as origin.

      ORGATI = true;
      A = C * DEL + Z[I] * Z[I] + Z[IP1] * Z[IP1];
      B = Z[I] * Z[I] * DEL;
      if (A > ZERO) {
        TAU = TWO * B / (A + sqrt((A * A - FOUR * B * C).abs()));
      } else {
        TAU = (A - sqrt((A * A - FOUR * B * C).abs())) / (TWO * C);
      }
      DLTLB = ZERO;
      DLTUB = MIDPT;
    } else {
      // (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)

      // We choose d(i+1) as origin.

      ORGATI = false;
      A = C * DEL - Z[I] * Z[I] - Z[IP1] * Z[IP1];
      B = Z[IP1] * Z[IP1] * DEL;
      if (A < ZERO) {
        TAU = TWO * B / (A - sqrt((A * A + FOUR * B * C).abs()));
      } else {
        TAU = -(A + sqrt((A * A + FOUR * B * C).abs())) / (TWO * C);
      }
      DLTLB = -MIDPT;
      DLTUB = ZERO;
    }

    if (ORGATI) {
      for (J = 1; J <= N; J++) {
        DELTA[J] = (D[J] - D[I]) - TAU;
      }
    } else {
      for (J = 1; J <= N; J++) {
        DELTA[J] = (D[J] - D[IP1]) - TAU;
      }
    }
    if (ORGATI) {
      II = I;
    } else {
      II = I + 1;
    }
    IIM1 = II - 1;
    IIP1 = II + 1;

    // Evaluate PSI and the derivative DPSI

    DPSI = ZERO;
    PSI = ZERO;
    ERRETM = ZERO;
    for (J = 1; J <= IIM1; J++) {
      TEMP = Z[J] / DELTA[J];
      PSI = PSI + Z[J] * TEMP;
      DPSI = DPSI + TEMP * TEMP;
      ERRETM = ERRETM + PSI;
    }
    ERRETM = (ERRETM).abs();

    // Evaluate PHI and the derivative DPHI

    DPHI = ZERO;
    PHI = ZERO;
    for (J = N; J >= IIP1; J--) {
      TEMP = Z[J] / DELTA[J];
      PHI = PHI + Z[J] * TEMP;
      DPHI = DPHI + TEMP * TEMP;
      ERRETM = ERRETM + PHI;
    }

    W = RHOINV + PHI + PSI;

    // W is the value of the secular function with
    // its ii-th element removed.

    SWTCH3 = false;
    if (ORGATI) {
      if (W < ZERO) SWTCH3 = true;
    } else {
      if (W > ZERO) SWTCH3 = true;
    }
    if (II == 1 || II == N) SWTCH3 = false;

    TEMP = Z[II] / DELTA[II];
    DW = DPSI + DPHI + TEMP * TEMP;
    TEMP = Z[II] * TEMP;
    W = W + TEMP;
    ERRETM = EIGHT * (PHI - PSI) +
        ERRETM +
        TWO * RHOINV +
        THREE * (TEMP).abs() +
        (TAU).abs() * DW;

    // Test for convergence

    if ((W).abs() <= EPS * ERRETM) {
      if (ORGATI) {
        DLAM.value = D[I] + TAU;
      } else {
        DLAM.value = D[IP1] + TAU;
      }
      return;
    }

    if (W <= ZERO) {
      DLTLB = max(DLTLB, TAU);
    } else {
      DLTUB = min(DLTUB, TAU);
    }

    // Calculate the new step

    NITER = NITER + 1;
    if (!SWTCH3) {
      if (ORGATI) {
        C = W - DELTA[IP1] * DW - (D[I] - D[IP1]) * pow((Z[I] / DELTA[I]), 2);
      } else {
        C = W - DELTA[I] * DW - (D[IP1] - D[I]) * pow((Z[IP1] / DELTA[IP1]), 2);
      }
      A = (DELTA[I] + DELTA[IP1]) * W - DELTA[I] * DELTA[IP1] * DW;
      B = DELTA[I] * DELTA[IP1] * W;
      if (C == ZERO) {
        if (A == ZERO) {
          if (ORGATI) {
            A = Z[I] * Z[I] + DELTA[IP1] * DELTA[IP1] * (DPSI + DPHI);
          } else {
            A = Z[IP1] * Z[IP1] + DELTA[I] * DELTA[I] * (DPSI + DPHI);
          }
        }
        ETA.value = B / A;
      } else if (A <= ZERO) {
        ETA.value = (A - sqrt((A * A - FOUR * B * C).abs())) / (TWO * C);
      } else {
        ETA.value = TWO * B / (A + sqrt((A * A - FOUR * B * C).abs()));
      }
    } else {
      // Interpolation using THREE most relevant poles

      TEMP = RHOINV + PSI + PHI;
      if (ORGATI) {
        TEMP1 = Z[IIM1] / DELTA[IIM1];
        TEMP1 = TEMP1 * TEMP1;
        C = TEMP - DELTA[IIP1] * (DPSI + DPHI) - (D[IIM1] - D[IIP1]) * TEMP1;
        ZZ[1] = Z[IIM1] * Z[IIM1];
        ZZ[3] = DELTA[IIP1] * DELTA[IIP1] * ((DPSI - TEMP1) + DPHI);
      } else {
        TEMP1 = Z[IIP1] / DELTA[IIP1];
        TEMP1 = TEMP1 * TEMP1;
        C = TEMP - DELTA[IIM1] * (DPSI + DPHI) - (D[IIP1] - D[IIM1]) * TEMP1;
        ZZ[1] = DELTA[IIM1] * DELTA[IIM1] * (DPSI + (DPHI - TEMP1));
        ZZ[3] = Z[IIP1] * Z[IIP1];
      }
      ZZ[2] = Z[II] * Z[II];
      dlaed6(NITER, ORGATI, C, DELTA(IIM1), ZZ, W, ETA, INFO);
      if (INFO.value != 0) return;
    }

    // Note, eta should be positive if w is negative, and
    // eta should be negative otherwise. However,
    // if for some reason caused by roundoff, eta*w > 0,
    // we simply use one Newton step instead. This way
    // will guarantee eta*w < 0.

    if (W * ETA.value >= ZERO) ETA.value = -W / DW;
    TEMP = TAU + ETA.value;
    if (TEMP > DLTUB || TEMP < DLTLB) {
      if (W < ZERO) {
        ETA.value = (DLTUB - TAU) / TWO;
      } else {
        ETA.value = (DLTLB - TAU) / TWO;
      }
    }

    PREW = W;

    for (J = 1; J <= N; J++) {
      DELTA[J] = DELTA[J] - ETA.value;
    }

    // Evaluate PSI and the derivative DPSI

    DPSI = ZERO;
    PSI = ZERO;
    ERRETM = ZERO;
    for (J = 1; J <= IIM1; J++) {
      TEMP = Z[J] / DELTA[J];
      PSI = PSI + Z[J] * TEMP;
      DPSI = DPSI + TEMP * TEMP;
      ERRETM = ERRETM + PSI;
    }
    ERRETM = (ERRETM).abs();

    // Evaluate PHI and the derivative DPHI

    DPHI = ZERO;
    PHI = ZERO;
    for (J = N; J >= IIP1; J--) {
      TEMP = Z[J] / DELTA[J];
      PHI = PHI + Z[J] * TEMP;
      DPHI = DPHI + TEMP * TEMP;
      ERRETM = ERRETM + PHI;
    }

    TEMP = Z[II] / DELTA[II];
    DW = DPSI + DPHI + TEMP * TEMP;
    TEMP = Z[II] * TEMP;
    W = RHOINV + PHI + PSI + TEMP;
    ERRETM = EIGHT * (PHI - PSI) +
        ERRETM +
        TWO * RHOINV +
        THREE * (TEMP).abs() +
        (TAU + ETA.value).abs() * DW;

    SWTCH = false;
    if (ORGATI) {
      if (-W > (PREW).abs() / TEN) SWTCH = true;
    } else {
      if (W > (PREW).abs() / TEN) SWTCH = true;
    }

    TAU = TAU + ETA.value;

    // Main loop to update the values of the array   DELTA

    ITER = NITER + 1;

    for (NITER = ITER; NITER <= MAXIT; NITER++) {
      // Test for convergence

      if ((W).abs() <= EPS * ERRETM) {
        if (ORGATI) {
          DLAM.value = D[I] + TAU;
        } else {
          DLAM.value = D[IP1] + TAU;
        }
        return;
      }

      if (W <= ZERO) {
        DLTLB = max(DLTLB, TAU);
      } else {
        DLTUB = min(DLTUB, TAU);
      }

      // Calculate the new step

      if (!SWTCH3) {
        if (!SWTCH) {
          if (ORGATI) {
            C = W -
                DELTA[IP1] * DW -
                (D[I] - D[IP1]) * pow((Z[I] / DELTA[I]), 2);
          } else {
            C = W -
                DELTA[I] * DW -
                (D[IP1] - D[I]) * pow((Z[IP1] / DELTA[IP1]), 2);
          }
        } else {
          TEMP = Z[II] / DELTA[II];
          if (ORGATI) {
            DPSI = DPSI + TEMP * TEMP;
          } else {
            DPHI = DPHI + TEMP * TEMP;
          }
          C = W - DELTA[I] * DPSI - DELTA[IP1] * DPHI;
        }
        A = (DELTA[I] + DELTA[IP1]) * W - DELTA[I] * DELTA[IP1] * DW;
        B = DELTA[I] * DELTA[IP1] * W;
        if (C == ZERO) {
          if (A == ZERO) {
            if (!SWTCH) {
              if (ORGATI) {
                A = Z[I] * Z[I] + DELTA[IP1] * DELTA[IP1] * (DPSI + DPHI);
              } else {
                A = Z[IP1] * Z[IP1] + DELTA[I] * DELTA[I] * (DPSI + DPHI);
              }
            } else {
              A = DELTA[I] * DELTA[I] * DPSI + DELTA[IP1] * DELTA[IP1] * DPHI;
            }
          }
          ETA.value = B / A;
        } else if (A <= ZERO) {
          ETA.value = (A - sqrt((A * A - FOUR * B * C).abs())) / (TWO * C);
        } else {
          ETA.value = TWO * B / (A + sqrt((A * A - FOUR * B * C).abs()));
        }
      } else {
        // Interpolation using THREE most relevant poles

        TEMP = RHOINV + PSI + PHI;
        if (SWTCH) {
          C = TEMP - DELTA[IIM1] * DPSI - DELTA[IIP1] * DPHI;
          ZZ[1] = DELTA[IIM1] * DELTA[IIM1] * DPSI;
          ZZ[3] = DELTA[IIP1] * DELTA[IIP1] * DPHI;
        } else {
          if (ORGATI) {
            TEMP1 = Z[IIM1] / DELTA[IIM1];
            TEMP1 = TEMP1 * TEMP1;
            C = TEMP -
                DELTA[IIP1] * (DPSI + DPHI) -
                (D[IIM1] - D[IIP1]) * TEMP1;
            ZZ[1] = Z[IIM1] * Z[IIM1];
            ZZ[3] = DELTA[IIP1] * DELTA[IIP1] * ((DPSI - TEMP1) + DPHI);
          } else {
            TEMP1 = Z[IIP1] / DELTA[IIP1];
            TEMP1 = TEMP1 * TEMP1;
            C = TEMP -
                DELTA[IIM1] * (DPSI + DPHI) -
                (D[IIP1] - D[IIM1]) * TEMP1;
            ZZ[1] = DELTA[IIM1] * DELTA[IIM1] * (DPSI + (DPHI - TEMP1));
            ZZ[3] = Z[IIP1] * Z[IIP1];
          }
        }
        dlaed6(NITER, ORGATI, C, DELTA(IIM1), ZZ, W, ETA, INFO);
        if (INFO.value != 0) return;
      }

      // Note, eta should be positive if w is negative, and
      // eta should be negative otherwise. However,
      // if for some reason caused by roundoff, eta*w > 0,
      // we simply use one Newton step instead. This way
      // will guarantee eta*w < 0.

      if (W * ETA.value >= ZERO) ETA.value = -W / DW;
      TEMP = TAU + ETA.value;
      if (TEMP > DLTUB || TEMP < DLTLB) {
        if (W < ZERO) {
          ETA.value = (DLTUB - TAU) / TWO;
        } else {
          ETA.value = (DLTLB - TAU) / TWO;
        }
      }

      for (J = 1; J <= N; J++) {
        DELTA[J] = DELTA[J] - ETA.value;
      }

      TAU = TAU + ETA.value;
      PREW = W;

      // Evaluate PSI and the derivative DPSI

      DPSI = ZERO;
      PSI = ZERO;
      ERRETM = ZERO;
      for (J = 1; J <= IIM1; J++) {
        TEMP = Z[J] / DELTA[J];
        PSI = PSI + Z[J] * TEMP;
        DPSI = DPSI + TEMP * TEMP;
        ERRETM = ERRETM + PSI;
      }
      ERRETM = (ERRETM).abs();

      // Evaluate PHI and the derivative DPHI

      DPHI = ZERO;
      PHI = ZERO;
      for (J = N; J >= IIP1; J--) {
        TEMP = Z[J] / DELTA[J];
        PHI = PHI + Z[J] * TEMP;
        DPHI = DPHI + TEMP * TEMP;
        ERRETM = ERRETM + PHI;
      }

      TEMP = Z[II] / DELTA[II];
      DW = DPSI + DPHI + TEMP * TEMP;
      TEMP = Z[II] * TEMP;
      W = RHOINV + PHI + PSI + TEMP;
      ERRETM = EIGHT * (PHI - PSI) +
          ERRETM +
          TWO * RHOINV +
          THREE * (TEMP).abs() +
          (TAU).abs() * DW;
      if (W * PREW > ZERO && (W).abs() > (PREW).abs() / TEN) SWTCH = !SWTCH;
    }

    // Return with INFO.value = 1, NITER = MAXIT and not converged

    INFO.value = 1;
    if (ORGATI) {
      DLAM.value = D[I] + TAU;
    } else {
      DLAM.value = D[IP1] + TAU;
    }
  }
}
