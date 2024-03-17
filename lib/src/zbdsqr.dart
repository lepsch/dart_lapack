import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdrot.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlas2.dart';
import 'package:lapack/src/dlasq1.dart';
import 'package:lapack/src/dlasv2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlasr.dart';

void zbdsqr(
  final String UPLO,
  final int N,
  final int NCVT,
  final int NRU,
  final int NCC,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> VT_,
  final int LDVT,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  final VT = VT_.having(ld: LDVT);
  final U = U_.having(ld: LDU);
  final C = C_.having(ld: LDC);
  final RWORK = RWORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  const ONE = 1.0;
  const NEGONE = -1.0;
  const HNDRTH = 0.01;
  const TEN = 10.0;
  const HNDRD = 100.0;
  const MEIGTH = -0.125;
  const MAXITR = 6;
  bool LOWER, ROTATE;
  int I,
      IDIR,
      ISUB,
      ITER,
      ITERDIVN,
      J,
      LL = 0,
      LLL,
      M,
      MAXITDIVN,
      NM1,
      NM12,
      NM13,
      OLDLL,
      OLDM;
  double ABSE,
      ABSS,
      EPS,
      F,
      G,
      H,
      MU,
      SLL,
      SMAX,
      SMIN,
      SMINOA,
      THRESH,
      TOL,
      TOLMUL,
      UNFL;
  final CS = Box(0.0),
      SN = Box(0.0),
      R = Box(0.0),
      SIGMN = Box(0.0),
      SIGMX = Box(0.0),
      SINR = Box(0.0),
      COSR = Box(0.0),
      SINL = Box(0.0),
      COSL = Box(0.0),
      SHIFT = Box(0.0),
      OLDCS = Box(0.0),
      OLDSN = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;
  LOWER = lsame(UPLO, 'L');
  if (!lsame(UPLO, 'U') && !LOWER) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NCVT < 0) {
    INFO.value = -3;
  } else if (NRU < 0) {
    INFO.value = -4;
  } else if (NCC < 0) {
    INFO.value = -5;
  } else if ((NCVT == 0 && LDVT < 1) || (NCVT > 0 && LDVT < max(1, N))) {
    INFO.value = -9;
  } else if (LDU < max(1, NRU)) {
    INFO.value = -11;
  } else if ((NCC == 0 && LDC < 1) || (NCC > 0 && LDC < max(1, N))) {
    INFO.value = -13;
  }
  if (INFO.value != 0) {
    xerbla('ZBDSQR', -INFO.value);
    return;
  }
  if (N == 0) return;
  if (N != 1) {
    // ROTATE is true if any singular vectors desired, false otherwise

    ROTATE = (NCVT > 0) || (NRU > 0) || (NCC > 0);

    // If no singular vectors desired, use qd algorithm

    if (!ROTATE) {
      dlasq1(N, D, E, RWORK, INFO);

      // If INFO.value equals 2, dqds didn't finish, try to finish

      if (INFO.value != 2) return;
      INFO.value = 0;
    }

    NM1 = N - 1;
    NM12 = NM1 + NM1;
    NM13 = NM12 + NM1;
    IDIR = 0;

    // Get machine constants

    EPS = dlamch('Epsilon');
    UNFL = dlamch('Safe minimum');

    // If matrix lower bidiagonal, rotate to be upper bidiagonal
    // by applying Givens rotations on the left

    if (LOWER) {
      for (I = 1; I <= N - 1; I++) {
        dlartg(D[I], E[I], CS, SN, R);
        D[I] = R.value;
        E[I] = SN.value * D[I + 1];
        D[I + 1] = CS.value * D[I + 1];
        RWORK[I] = CS.value;
        RWORK[NM1 + I] = SN.value;
      }

      // Update singular vectors if desired

      if (NRU > 0) {
        zlasr('R.value', 'V', 'F', NRU, N, RWORK(1), RWORK(N), U, LDU);
      }
      if (NCC > 0) zlasr('L', 'V', 'F', N, NCC, RWORK(1), RWORK(N), C, LDC);
    }

    // Compute singular values to relative accuracy TOL
    // (By setting TOL to be negative, algorithm will compute
    // singular values to absolute accuracy ABS(TOL)*norm(input matrix))

    TOLMUL = max(TEN, min(HNDRD, pow(EPS, MEIGTH).toDouble()));
    TOL = TOLMUL * EPS;

    // Compute approximate maximum, minimum singular values

    SMAX = ZERO;
    for (I = 1; I <= N; I++) {
      SMAX = max(SMAX, (D[I]).abs());
    }
    for (I = 1; I <= N - 1; I++) {
      SMAX = max(SMAX, (E[I]).abs());
    }
    SMIN = ZERO;
    if (TOL >= ZERO) {
      // Relative accuracy desired

      SMINOA = (D[1]).abs();
      if (SMINOA != ZERO) {
        MU = SMINOA;
        for (I = 2; I <= N; I++) {
          MU = (D[I]).abs() * (MU / (MU + (E[I - 1]).abs()));
          SMINOA = min(SMINOA, MU);
          if (SMINOA == ZERO) break;
        }
      }
      SMINOA = SMINOA / sqrt(N.toDouble());
      THRESH = max(TOL * SMINOA, MAXITR * (N * (N * UNFL)));
    } else {
      // Absolute accuracy desired

      THRESH = max((TOL).abs() * SMAX, MAXITR * (N * (N * UNFL)));
    }

    // Prepare for main iteration loop for the singular values
    // (MAXIT is the maximum number of passes through the inner
    // loop permitted before nonconvergence signalled.)

    MAXITDIVN = MAXITR * N;
    ITERDIVN = 0;
    ITER = -1;
    OLDLL = -1;
    OLDM = -1;

    // M points to last element of unconverged part of matrix

    M = N;

    // Begin main iteration loop

    mainLoop:
    while (true) {
      // Check for convergence or exceeding iteration count

      if (M <= 1) break;
      if (ITER >= N) {
        ITER -= N;
        ITERDIVN++;
        if (ITERDIVN >= MAXITDIVN) {
          // Maximum number of iterations exceeded, failure to converge
          INFO.value = 0;
          for (I = 1; I <= N - 1; I++) {
            if (E[I] != ZERO) INFO.value++;
          }
          return;
        }
      }

      // Find diagonal block of matrix to work on

      if (TOL < ZERO && (D[M]).abs() <= THRESH) D[M] = ZERO;
      SMAX = (D[M]).abs();
      var diagonalBlockFound = false;
      for (LLL = 1; LLL <= M - 1; LLL++) {
        LL = M - LLL;
        ABSS = (D[LL]).abs();
        ABSE = (E[LL]).abs();
        if (TOL < ZERO && ABSS <= THRESH) D[LL] = ZERO;
        if (ABSE <= THRESH) {
          diagonalBlockFound = true;
          break;
        }
        SMAX = max(SMAX, max(ABSS, ABSE));
      }
      if (!diagonalBlockFound) {
        LL = 0;
      } else {
        E[LL] = ZERO;

        // Matrix splits since E(LL) = 0

        if (LL == M - 1) {
          // Convergence of bottom singular value, return to top of loop

          M--;
          continue;
        }
      }

      LL++;

      // E(LL) through E(M-1) are nonzero, E(LL-1) is zero

      if (LL == M - 1) {
        // 2 by 2 block, handle separately

        dlasv2(D[M - 1], E[M - 1], D[M], SIGMN, SIGMX, SINR, COSR, SINL, COSL);
        D[M - 1] = SIGMX.value;
        E[M - 1] = ZERO;
        D[M] = SIGMN.value;

        // Compute singular vectors, if desired

        if (NCVT > 0) {
          zdrot(NCVT, VT(M - 1, 1).asArray(), LDVT, VT(M, 1).asArray(), LDVT,
              COSR.value, SINR.value);
        }
        if (NRU > 0) {
          zdrot(NRU, U(1, M - 1).asArray(), 1, U(1, M).asArray(), 1, COSL.value,
              SINL.value);
        }
        if (NCC > 0) {
          zdrot(NCC, C(M - 1, 1).asArray(), LDC, C(M, 1).asArray(), LDC,
              COSL.value, SINL.value);
        }
        M -= 2;
        continue mainLoop;
      }

      // If working on new submatrix, choose shift direction
      // (from larger end diagonal element towards smaller)

      if (LL > OLDM || M < OLDLL) {
        if ((D[LL]).abs() >= (D[M]).abs()) {
          // Chase bulge from top (big end) to bottom (small end)

          IDIR = 1;
        } else {
          // Chase bulge from bottom (big end) to top (small end)

          IDIR = 2;
        }
      }

      // Apply convergence tests

      if (IDIR == 1) {
        // Run convergence test in forward direction
        // First apply standard test to bottom of matrix

        if ((E[M - 1]).abs() <= (TOL).abs() * (D[M]).abs() ||
            (TOL < ZERO && (E[M - 1]).abs() <= THRESH)) {
          E[M - 1] = ZERO;
          continue mainLoop;
        }

        if (TOL >= ZERO) {
          // If relative accuracy desired,
          // apply convergence criterion forward

          MU = (D[LL]).abs();
          SMIN = MU;
          for (LLL = LL; LLL <= M - 1; LLL++) {
            if ((E[LLL]).abs() <= TOL * MU) {
              E[LLL] = ZERO;
              continue mainLoop;
            }
            MU = (D[LLL + 1]).abs() * (MU / (MU + (E[LLL]).abs()));
            SMIN = min(SMIN, MU);
          }
        }
      } else {
        // Run convergence test in backward direction
        // First apply standard test to top of matrix

        if ((E[LL]).abs() <= (TOL).abs() * (D[LL]).abs() ||
            (TOL < ZERO && (E[LL]).abs() <= THRESH)) {
          E[LL] = ZERO;
          continue mainLoop;
        }

        if (TOL >= ZERO) {
          // If relative accuracy desired,
          // apply convergence criterion backward

          MU = (D[M]).abs();
          SMIN = MU;
          for (LLL = M - 1; LLL >= LL; LLL--) {
            if ((E[LLL]).abs() <= TOL * MU) {
              E[LLL] = ZERO;
              continue mainLoop;
            }
            MU = (D[LLL]).abs() * (MU / (MU + (E[LLL]).abs()));
            SMIN = min(SMIN, MU);
          }
        }
      }
      OLDLL = LL;
      OLDM = M;

      // Compute shift.  First, test if shifting would ruin relative
      // accuracy, and if so set the shift to zero.

      if (TOL >= ZERO && N * TOL * (SMIN / SMAX) <= max(EPS, HNDRTH * TOL)) {
        // Use a zero shift to avoid loss of relative accuracy

        SHIFT.value = ZERO;
      } else {
        // Compute the shift from 2-by-2 block at end of matrix

        if (IDIR == 1) {
          SLL = (D[LL]).abs();
          dlas2(D[M - 1], E[M - 1], D[M], SHIFT, R);
        } else {
          SLL = (D[M]).abs();
          dlas2(D[LL], E[LL], D[LL + 1], SHIFT, R);
        }

        // Test if shift negligible, and if so set to zero

        if (SLL > ZERO) {
          if (pow(SHIFT.value / SLL, 2) < EPS) SHIFT.value = ZERO;
        }
      }

      // Increment iteration count

      ITER += M - LL;

      // If SHIFT.value = 0, do simplified QR iteration

      if (SHIFT.value == ZERO) {
        if (IDIR == 1) {
          // Chase bulge from top to bottom
          // Save cosines and sines for later singular vector updates

          CS.value = ONE;
          OLDCS.value = ONE;
          for (I = LL; I <= M - 1; I++) {
            dlartg(D[I] * CS.value, E[I], CS, SN, R);
            if (I > LL) E[I - 1] = OLDSN.value * R.value;
            dlartg(OLDCS.value * R.value, D[I + 1] * SN.value, OLDCS, OLDSN,
                D.box(I));
            RWORK[I - LL + 1] = CS.value;
            RWORK[I - LL + 1 + NM1] = SN.value;
            RWORK[I - LL + 1 + NM12] = OLDCS.value;
            RWORK[I - LL + 1 + NM13] = OLDSN.value;
          }
          H = D[M] * CS.value;
          D[M] = H * OLDCS.value;
          E[M - 1] = H * OLDSN.value;

          // Update singular vectors

          if (NCVT > 0) {
            zlasr('L', 'V', 'F', M - LL + 1, NCVT, RWORK(1), RWORK(N),
                VT(LL, 1), LDVT);
          }
          if (NRU > 0) {
            zlasr('R.value', 'V', 'F', NRU, M - LL + 1, RWORK(NM12 + 1),
                RWORK(NM13 + 1), U(1, LL), LDU);
          }
          if (NCC > 0) {
            zlasr('L', 'V', 'F', M - LL + 1, NCC, RWORK(NM12 + 1),
                RWORK(NM13 + 1), C(LL, 1), LDC);
          }

          // Test convergence

          if ((E[M - 1]).abs() <= THRESH) E[M - 1] = ZERO;
        } else {
          // Chase bulge from bottom to top
          // Save cosines and sines for later singular vector updates

          CS.value = ONE;
          OLDCS.value = ONE;
          for (I = M; I >= LL + 1; I--) {
            // 130
            dlartg(D[I] * CS.value, E[I - 1], CS, SN, R);
            if (I < M) E[I] = OLDSN.value * R.value;
            dlartg(OLDCS.value * R.value, D[I - 1] * SN.value, OLDCS, OLDSN,
                D.box(I));
            RWORK[I - LL] = CS.value;
            RWORK[I - LL + NM1] = -SN.value;
            RWORK[I - LL + NM12] = OLDCS.value;
            RWORK[I - LL + NM13] = -OLDSN.value;
          } // 130
          H = D[LL] * CS.value;
          D[LL] = H * OLDCS.value;
          E[LL] = H * OLDSN.value;

          // Update singular vectors

          if (NCVT > 0) {
            zlasr('L', 'V', 'B', M - LL + 1, NCVT, RWORK(NM12 + 1),
                RWORK(NM13 + 1), VT(LL, 1), LDVT);
          }
          if (NRU > 0) {
            zlasr('R.value', 'V', 'B', NRU, M - LL + 1, RWORK(1), RWORK(N),
                U(1, LL), LDU);
          }
          if (NCC > 0) {
            zlasr('L', 'V', 'B', M - LL + 1, NCC, RWORK(1), RWORK(N), C(LL, 1),
                LDC);
          }

          // Test convergence

          if ((E[LL]).abs() <= THRESH) E[LL] = ZERO;
        }
      } else {
        // Use nonzero shift

        if (IDIR == 1) {
          // Chase bulge from top to bottom
          // Save cosines and sines for later singular vector updates

          F = ((D[LL]).abs() - SHIFT.value) *
              (sign(ONE, D[LL]) + SHIFT.value / D[LL]);
          G = E[LL];
          for (I = LL; I <= M - 1; I++) {
            // 140
            dlartg(F, G, COSR, SINR, R);
            if (I > LL) E[I - 1] = R.value;
            F = COSR.value * D[I] + SINR.value * E[I];
            E[I] = COSR.value * E[I] - SINR.value * D[I];
            G = SINR.value * D[I + 1];
            D[I + 1] = COSR.value * D[I + 1];
            dlartg(F, G, COSL, SINL, R);
            D[I] = R.value;
            F = COSL.value * E[I] + SINL.value * D[I + 1];
            D[I + 1] = COSL.value * D[I + 1] - SINL.value * E[I];
            if (I < M - 1) {
              G = SINL.value * E[I + 1];
              E[I + 1] = COSL.value * E[I + 1];
            }
            RWORK[I - LL + 1] = COSR.value;
            RWORK[I - LL + 1 + NM1] = SINR.value;
            RWORK[I - LL + 1 + NM12] = COSL.value;
            RWORK[I - LL + 1 + NM13] = SINL.value;
          } // 140
          E[M - 1] = F;

          // Update singular vectors

          if (NCVT > 0) {
            zlasr('L', 'V', 'F', M - LL + 1, NCVT, RWORK(1), RWORK(N),
                VT(LL, 1), LDVT);
          }
          if (NRU > 0) {
            zlasr('R.value', 'V', 'F', NRU, M - LL + 1, RWORK(NM12 + 1),
                RWORK(NM13 + 1), U(1, LL), LDU);
          }
          if (NCC > 0) {
            zlasr('L', 'V', 'F', M - LL + 1, NCC, RWORK(NM12 + 1),
                RWORK(NM13 + 1), C(LL, 1), LDC);
          }

          // Test convergence

          if ((E[M - 1]).abs() <= THRESH) E[M - 1] = ZERO;
        } else {
          // Chase bulge from bottom to top
          // Save cosines and sines for later singular vector updates

          F = ((D[M]).abs() - SHIFT.value) *
              (sign(ONE, D[M]) + SHIFT.value / D[M]);
          G = E[M - 1];
          for (I = M; I >= LL + 1; I--) {
            // 150
            dlartg(F, G, COSR, SINR, R);
            if (I < M) E[I] = R.value;
            F = COSR.value * D[I] + SINR.value * E[I - 1];
            E[I - 1] = COSR.value * E[I - 1] - SINR.value * D[I];
            G = SINR.value * D[I - 1];
            D[I - 1] = COSR.value * D[I - 1];
            dlartg(F, G, COSL, SINL, R);
            D[I] = R.value;
            F = COSL.value * E[I - 1] + SINL.value * D[I - 1];
            D[I - 1] = COSL.value * D[I - 1] - SINL.value * E[I - 1];
            if (I > LL + 1) {
              G = SINL.value * E[I - 2];
              E[I - 2] = COSL.value * E[I - 2];
            }
            RWORK[I - LL] = COSR.value;
            RWORK[I - LL + NM1] = -SINR.value;
            RWORK[I - LL + NM12] = COSL.value;
            RWORK[I - LL + NM13] = -SINL.value;
          } // 150
          E[LL] = F;

          // Test convergence

          if ((E[LL]).abs() <= THRESH) E[LL] = ZERO;

          // Update singular vectors if desired

          if (NCVT > 0) {
            zlasr('L', 'V', 'B', M - LL + 1, NCVT, RWORK(NM12 + 1),
                RWORK(NM13 + 1), VT(LL, 1), LDVT);
          }
          if (NRU > 0) {
            zlasr('R.value', 'V', 'B', NRU, M - LL + 1, RWORK(1), RWORK(N),
                U(1, LL), LDU);
          }
          if (NCC > 0) {
            zlasr('L', 'V', 'B', M - LL + 1, NCC, RWORK(1), RWORK(N), C(LL, 1),
                LDC);
          }
        }
      }

      // QR iteration finished, go back and check convergence
    }
  }

  // All singular values converged, so make them positive

  for (I = 1; I <= N; I++) {
    if (D[I] < ZERO) {
      D[I] = -D[I];

      // Change sign of singular vectors, if desired

      if (NCVT > 0) zdscal(NCVT, NEGONE, VT(I, 1).asArray(), LDVT);
    }
  }

  // Sort the singular values into decreasing order (insertion sort on
  // singular values, but only one transposition per singular vector)

  for (I = 1; I <= N - 1; I++) {
    // Scan for smallest D(I)

    ISUB = 1;
    SMIN = D[1];
    for (J = 2; J <= N + 1 - I; J++) {
      // 180
      if (D[J] <= SMIN) {
        ISUB = J;
        SMIN = D[J];
      }
    } // 180
    if (ISUB != N + 1 - I) {
      // Swap singular values and vectors

      D[ISUB] = D[N + 1 - I];
      D[N + 1 - I] = SMIN;
      if (NCVT > 0) {
        zswap(NCVT, VT(ISUB, 1).asArray(), LDVT, VT(N + 1 - I, 1).asArray(),
            LDVT);
      }
      if (NRU > 0) {
        zswap(NRU, U(1, ISUB).asArray(), 1, U(1, N + 1 - I).asArray(), 1);
      }
      if (NCC > 0) {
        zswap(NCC, C(ISUB, 1).asArray(), LDC, C(N + 1 - I, 1).asArray(), LDC);
      }
    }
  }
}
