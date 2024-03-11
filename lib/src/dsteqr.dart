import 'dart:math';

import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlae2.dart';
import 'package:lapack/src/dlaev2.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dlasr.dart';
import 'package:lapack/src/dlasrt.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsteqr(
  final String COMPZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0;
  const MAXIT = 30;
  int I,
      ICOMPZ,
      II,
      ISCALE,
      J,
      JTOT,
      K,
      L,
      L1,
      LEND,
      LENDM1,
      LENDP1,
      LENDSV,
      LM1,
      LSV,
      M = 0,
      MM,
      MM1,
      NM1,
      NMAXIT;
  double ANORM, B, EPS, EPS2, F, G, P, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST;
  final C = Box(0.0),
      S = Box(0.0),
      R = Box(0.0),
      RT1 = Box(0.0),
      RT2 = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;

  if (lsame(COMPZ, 'N')) {
    ICOMPZ = 0;
  } else if (lsame(COMPZ, 'V')) {
    ICOMPZ = 1;
  } else if (lsame(COMPZ, 'I')) {
    ICOMPZ = 2;
  } else {
    ICOMPZ = -1;
  }
  if (ICOMPZ < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if ((LDZ < 1) || (ICOMPZ > 0 && LDZ < max(1, N))) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DSTEQR', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    if (ICOMPZ == 2) Z[1][1] = ONE;
    return;
  }

  // Determine the unit roundoff and over/underflow thresholds.

  EPS = dlamch('E');
  EPS2 = pow(EPS, 2).toDouble();
  SAFMIN = dlamch('S');
  SAFMAX = ONE / SAFMIN;
  SSFMAX = sqrt(SAFMAX) / THREE;
  SSFMIN = sqrt(SAFMIN) / EPS2;

  // Compute the eigenvalues and eigenvectors of the tridiagonal
  // matrix.

  if (ICOMPZ == 2) dlaset('Full', N, N, ZERO, ONE, Z, LDZ);

  NMAXIT = N * MAXIT;
  JTOT = 0;

  // Determine where the matrix splits and choose QL or QR iteration
  // for each block, according to whether top or bottom diagonal
  // element is smaller.

  L1 = 1;
  NM1 = N - 1;
  var sortEigen = false;
  while (true) {
    if (L1 > N) {
      sortEigen = true;
      break;
    }
    if (L1 > 1) E[L1 - 1] = ZERO;
    var isSmall = false;
    if (L1 <= NM1) {
      for (M = L1; M <= NM1; M++) {
        TST = (E[M]).abs();
        if (TST == ZERO) {
          isSmall = true;
          break;
        }
        if (TST <= (sqrt((D[M]).abs()) * sqrt((D[M + 1]).abs())) * EPS) {
          E[M] = ZERO;
          isSmall = true;
          break;
        }
      }
    }
    if (!isSmall) {
      M = N;
    }

    L = L1;
    LSV = L;
    LEND = M;
    LENDSV = LEND;
    L1 = M + 1;
    if (LEND == L) continue;

    // Scale submatrix in rows and columns L to LEND

    ANORM = dlanst('M', LEND - L + 1, D(L), E(L));
    ISCALE = 0;
    if (ANORM == ZERO) continue;
    if (ANORM > SSFMAX) {
      ISCALE = 1;
      dlascl(
          'G', 0, 0, ANORM, SSFMAX, LEND - L + 1, 1, D(L).asMatrix(N), N, INFO);
      dlascl('G', 0, 0, ANORM, SSFMAX, LEND - L, 1, E(L).asMatrix(N), N, INFO);
    } else if (ANORM < SSFMIN) {
      ISCALE = 2;
      dlascl(
          'G', 0, 0, ANORM, SSFMIN, LEND - L + 1, 1, D(L).asMatrix(N), N, INFO);
      dlascl('G', 0, 0, ANORM, SSFMIN, LEND - L, 1, E(L).asMatrix(N), N, INFO);
    }

    // Choose between QL and QR iteration

    if ((D[LEND]).abs() < (D[L]).abs()) {
      LEND = LSV;
      L = LENDSV;
    }

    if (LEND > L) {
      // QL Iteration

      // Look for small subdiagonal element.

      while (true) {
        var hasSmallSubdiagonalElement = false;
        if (L != LEND) {
          LENDM1 = LEND - 1;
          for (M = L; M <= LENDM1; M++) {
            TST = pow((E[M]).abs(), 2).toDouble();
            if (TST <= (EPS2 * (D[M]).abs()) * (D[M + 1]).abs() + SAFMIN) {
              hasSmallSubdiagonalElement = true;
              break;
            }
          }
        }
        if (!hasSmallSubdiagonalElement) {
          M = LEND;
        }

        if (M < LEND) E[M] = ZERO;
        P = D[L];
        if (M != L) {
          // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
          // to compute its eigensystem.

          if (M == L + 1) {
            if (ICOMPZ > 0) {
              dlaev2(D[L], E[L], D[L + 1], RT1, RT2, C, S);
              WORK[L] = C.value;
              WORK[N - 1 + L] = S.value;
              dlasr(
                  'R', 'V', 'B', N, 2, WORK(L), WORK(N - 1 + L), Z(1, L), LDZ);
            } else {
              dlae2(D[L], E[L], D[L + 1], RT1, RT2);
            }
            D[L] = RT1.value;
            D[L + 1] = RT2.value;
            E[L] = ZERO;
            L += 2;
            if (L <= LEND) continue;
            break;
          }
          if (JTOT == NMAXIT) break;

          JTOT++;

          // Form shift.

          G = (D[L + 1] - P) / (TWO * E[L]);
          R.value = dlapy2(G, ONE);
          G = D[M] - P + (E[L] / (G + sign(R.value, G)));

          S.value = ONE;
          C.value = ONE;
          P = ZERO;

          // Inner loop

          MM1 = M - 1;
          for (I = MM1; I >= L; I--) {
            F = S.value * E[I];
            B = C.value * E[I];
            dlartg(G, F, C, S, R);
            if (I != M - 1) E[I + 1] = R.value;
            G = D[I + 1] - P;
            R.value = (D[I] - G) * S.value + TWO * C.value * B;
            P = S.value * R.value;
            D[I + 1] = G + P;
            G = C.value * R.value - B;

            // If eigenvectors are desired, then save rotations.

            if (ICOMPZ > 0) {
              WORK[I] = C.value;
              WORK[N - 1 + I] = -S.value;
            }
          }

          // If eigenvectors are desired, then apply saved rotations.

          if (ICOMPZ > 0) {
            MM = M - L + 1;
            dlasr('R', 'V', 'B', N, MM, WORK(L), WORK(N - 1 + L), Z(1, L), LDZ);
          }

          D[L] = D[L] - P;
          E[L] = G;
          continue;
        }
        // Eigenvalue found.
        D[L] = P;
        L++;
        if (L <= LEND) continue;
        break;
      }
    } else {
      // QR Iteration

      // Look for small superdiagonal element.

      while (true) {
        var hasSmallSuperdiagonalElement = false;
        if (L != LEND) {
          LENDP1 = LEND + 1;
          for (M = L; M >= LENDP1; M--) {
            TST = pow((E[M - 1]).abs(), 2).toDouble();
            if (TST <= (EPS2 * (D[M]).abs()) * (D[M - 1]).abs() + SAFMIN) {
              hasSmallSuperdiagonalElement = true;
              break;
            }
          }
        }
        if (!hasSmallSuperdiagonalElement) {
          M = LEND;
        }

        if (M > LEND) E[M - 1] = ZERO;
        P = D[L];
        if (M == L) {
          // Eigenvalue found.
          D[L] = P;
          L--;
          if (L >= LEND) continue;
        } else

        // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
        // to compute its eigensystem.

        if (M == L - 1) {
          if (ICOMPZ > 0) {
            dlaev2(D[L - 1], E[L - 1], D[L], RT1, RT2, C, S);
            WORK[M] = C.value;
            WORK[N - 1 + M] = S.value;
            dlasr('R', 'V', 'F', N, 2, WORK(M), WORK(N - 1 + M), Z(1, L - 1),
                LDZ);
          } else {
            dlae2(D[L - 1], E[L - 1], D[L], RT1, RT2);
          }
          D[L - 1] = RT1.value;
          D[L] = RT2.value;
          E[L - 1] = ZERO;
          L -= 2;
          if (L >= LEND) continue;
        } else if (JTOT != NMAXIT) {
          JTOT++;

          // Form shift.

          G = (D[L - 1] - P) / (TWO * E[L - 1]);
          R.value = dlapy2(G, ONE);
          G = D[M] - P + (E[L - 1] / (G + sign(R.value, G)));

          S.value = ONE;
          C.value = ONE;
          P = ZERO;

          // Inner loop

          LM1 = L - 1;
          for (I = M; I <= LM1; I++) {
            F = S.value * E[I];
            B = C.value * E[I];
            dlartg(G, F, C, S, R);
            if (I != M) E[I - 1] = R.value;
            G = D[I] - P;
            R.value = (D[I + 1] - G) * S.value + TWO * C.value * B;
            P = S.value * R.value;
            D[I] = G + P;
            G = C.value * R.value - B;

            // If eigenvectors are desired, then save rotations.

            if (ICOMPZ > 0) {
              WORK[I] = C.value;
              WORK[N - 1 + I] = S.value;
            }
          }

          // If eigenvectors are desired, then apply saved rotations.

          if (ICOMPZ > 0) {
            MM = L - M + 1;
            dlasr('R', 'V', 'F', N, MM, WORK(M), WORK(N - 1 + M), Z(1, M), LDZ);
          }

          D[L] = D[L] - P;
          E[LM1] = G;
          continue;
        }
        break;
      }
    }

    // Undo scaling if necessary

    if (ISCALE == 1) {
      dlascl('G', 0, 0, SSFMAX, ANORM, LENDSV - LSV + 1, 1, D(LSV).asMatrix(N),
          N, INFO);
      dlascl('G', 0, 0, SSFMAX, ANORM, LENDSV - LSV, 1, E(LSV).asMatrix(N), N,
          INFO);
    } else if (ISCALE == 2) {
      dlascl('G', 0, 0, SSFMIN, ANORM, LENDSV - LSV + 1, 1, D(LSV).asMatrix(N),
          N, INFO);
      dlascl('G', 0, 0, SSFMIN, ANORM, LENDSV - LSV, 1, E(LSV).asMatrix(N), N,
          INFO);
    }

    // Check for no convergence to an eigenvalue after a total
    // of N*MAXIT iterations.

    if (JTOT >= NMAXIT) break;
  }
  if (!sortEigen) {
    for (I = 1; I <= N - 1; I++) {
      if (E[I] != ZERO) INFO.value = INFO.value + 1;
    }
    return;
  }

  // Order eigenvalues and eigenvectors.

  if (ICOMPZ == 0) {
    // Use Quick Sort

    dlasrt('I', N, D, INFO);
  } else {
    // Use Selection Sort to minimize swaps of eigenvectors

    for (II = 2; II <= N; II++) {
      I = II - 1;
      K = I;
      P = D[I];
      for (J = II; J <= N; J++) {
        if (D[J] < P) {
          K = J;
          P = D[J];
        }
      }
      if (K != I) {
        D[K] = D[I];
        D[I] = P;
        dswap(N, Z(1, I).asArray(), 1, Z(1, K).asArray(), 1);
      }
    }
  }
}
