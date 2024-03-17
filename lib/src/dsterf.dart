import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlae2.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasrt.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsterf(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0;
  const MAXIT = 30;
  int I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, NMAXIT;
  double ALPHA,
      ANORM,
      BB,
      C,
      EPS,
      EPS2,
      GAMMA,
      OLDC,
      OLDGAM,
      P,
      R,
      RTE,
      S,
      SAFMAX,
      SAFMIN,
      SIGMA,
      SSFMAX,
      SSFMIN;
  // RMAX;
  final RT1 = Box(0.0), RT2 = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;

  // Quick return if possible

  if (N < 0) {
    INFO.value = -1;
    xerbla('DSTERF', -INFO.value);
    return;
  }
  if (N <= 1) return;

  // Determine the unit roundoff for this environment.

  EPS = dlamch('E');
  EPS2 = pow(EPS, 2).toDouble();
  SAFMIN = dlamch('S');
  SAFMAX = ONE / SAFMIN;
  SSFMAX = sqrt(SAFMAX) / THREE;
  SSFMIN = sqrt(SAFMIN) / EPS2;
  // RMAX = dlamch('O');

  // Compute the eigenvalues of the tridiagonal matrix.

  NMAXIT = N * MAXIT;
  SIGMA = ZERO;
  JTOT = 0;

  // Determine where the matrix splits and choose QL or QR iteration
  // for each block, according to whether top or bottom diagonal
  // element is smaller.

  L1 = 1;

  while (true) {
    if (L1 > N) break;
    if (L1 > 1) E[L1 - 1] = ZERO;
    bool isLoopExhausted = true;
    for (M = L1; M <= N - 1; M++) {
      if ((E[M]).abs() <= (sqrt((D[M]).abs()) * sqrt((D[M + 1]).abs())) * EPS) {
        E[M] = ZERO;
        isLoopExhausted = false;
        break;
      }
    }
    if (isLoopExhausted) {
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
    if ((ANORM > SSFMAX)) {
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

    for (I = L; I <= LEND - 1; I++) {
      E[I] = pow(E[I], 2).toDouble();
    }

    // Choose between QL and QR iteration

    if ((D[LEND]).abs() < (D[L]).abs()) {
      LEND = LSV;
      L = LENDSV;
    }

    if (LEND >= L) {
      // QL Iteration

      // Look for small subdiagonal element.

      while (LEND >= L) {
        var hasSubdiagonalElement = false;
        if (L != LEND) {
          for (M = L; M <= LEND - 1; M++) {
            if ((E[M]).abs() <= EPS2 * (D[M] * D[M + 1]).abs()) {
              hasSubdiagonalElement = true;
              break;
            }
          }
        }
        if (!hasSubdiagonalElement) {
          M = LEND;
        }

        if (M < LEND) E[M] = ZERO;
        P = D[L];
        if (M != L) {
          // If remaining matrix is 2 by 2, use DLAE2 to compute its
          // eigenvalues.

          if (M == L + 1) {
            RTE = sqrt(E[L]);
            dlae2(D[L], RTE, D[L + 1], RT1, RT2);
            D[L] = RT1.value;
            D[L + 1] = RT2.value;
            E[L] = ZERO;
            L += 2;
            continue;
          }

          if (JTOT == NMAXIT) break;
          JTOT++;

          // Form shift.

          RTE = sqrt(E[L]);
          SIGMA = (D[L + 1] - P) / (TWO * RTE);
          R = dlapy2(SIGMA, ONE);
          SIGMA = P - (RTE / (SIGMA + sign(R, SIGMA)));

          C = ONE;
          S = ZERO;
          GAMMA = D[M] - SIGMA;
          P = GAMMA * GAMMA;

          // Inner loop

          for (I = M - 1; I >= L; I--) {
            BB = E[I];
            R = P + BB;
            if (I != M - 1) E[I + 1] = S * R;
            OLDC = C;
            C = P / R;
            S = BB / R;
            OLDGAM = GAMMA;
            ALPHA = D[I];
            GAMMA = C * (ALPHA - SIGMA) - S * OLDGAM;
            D[I + 1] = OLDGAM + (ALPHA - GAMMA);
            if (C != ZERO) {
              P = (GAMMA * GAMMA) / C;
            } else {
              P = OLDC * BB;
            }
          }

          E[L] = S * P;
          D[L] = SIGMA + GAMMA;
          continue;
        }

        // Eigenvalue found.

        D[L] = P;
        L++;
      }
    } else {
      // QR Iteration

      // Look for small superdiagonal element.

      while (L >= LEND) {
        var hasSuperdiagonalElement = false;
        for (M = L; M >= LEND + 1; M--) {
          if ((E[M - 1]).abs() <= EPS2 * (D[M] * D[M - 1]).abs()) {
            hasSuperdiagonalElement = true;
            break;
          }
        }
        if (!hasSuperdiagonalElement) {
          M = LEND;
        }

        if (M > LEND) E[M - 1] = ZERO;
        P = D[L];
        if (M != L) {
          // If remaining matrix is 2 by 2, use DLAE2 to compute its
          // eigenvalues.

          if (M == L - 1) {
            RTE = sqrt(E[L - 1]);
            dlae2(D[L], RTE, D[L - 1], RT1, RT2);
            D[L] = RT1.value;
            D[L - 1] = RT2.value;
            E[L - 1] = ZERO;
            L -= 2;
            continue;
          }

          if (JTOT == NMAXIT) break;
          JTOT++;

          // Form shift.

          RTE = sqrt(E[L - 1]);
          SIGMA = (D[L - 1] - P) / (TWO * RTE);
          R = dlapy2(SIGMA, ONE);
          SIGMA = P - (RTE / (SIGMA + sign(R, SIGMA)));

          C = ONE;
          S = ZERO;
          GAMMA = D[M] - SIGMA;
          P = GAMMA * GAMMA;

          // Inner loop

          for (I = M; I <= L - 1; I++) {
            BB = E[I];
            R = P + BB;
            if (I != M) E[I - 1] = S * R;
            OLDC = C;
            C = P / R;
            S = BB / R;
            OLDGAM = GAMMA;
            ALPHA = D[I + 1];
            GAMMA = C * (ALPHA - SIGMA) - S * OLDGAM;
            D[I] = OLDGAM + (ALPHA - GAMMA);
            if (C != ZERO) {
              P = (GAMMA * GAMMA) / C;
            } else {
              P = OLDC * BB;
            }
          }

          E[L - 1] = S * P;
          D[L] = SIGMA + GAMMA;
          continue;
        }

        // Eigenvalue found.
        D[L] = P;
        L--;
      }
    }

    // Undo scaling if necessary
    if (ISCALE == 1) {
      dlascl('G', 0, 0, SSFMAX, ANORM, LENDSV - LSV + 1, 1, D(LSV).asMatrix(N),
          N, INFO);
    }
    if (ISCALE == 2) {
      dlascl('G', 0, 0, SSFMIN, ANORM, LENDSV - LSV + 1, 1, D(LSV).asMatrix(N),
          N, INFO);
    }

    // Check for no convergence to an eigenvalue after a total
    // of N*MAXIT iterations.

    if (JTOT >= NMAXIT) {
      for (I = 1; I <= N - 1; I++) {
        if (E[I] != ZERO) INFO.value++;
      }
      return;
    }
  }

  // Sort eigenvalues in increasing order.

  dlasrt('I', N, D, INFO);
}
