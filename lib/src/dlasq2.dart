// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlasq3.dart';
import 'package:dart_lapack/src/dlasrt.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlasq2(final int N, final Array<double> Z_, final Box<int> INFO) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.having();
  const CBIAS = 1.50;
  const ZERO = 0.0,
      HALF = 0.5,
      ONE = 1.0,
      TWO = 2.0,
      FOUR = 4.0,
      HUNDRD = 100.0;
  bool IEEE;
  int I0, I1, I4 = 0, IPN4, IWHILA, IWHILB, K, KMIN, N1, NBIG, SPLT;
  double D,
      DEE,
      DEEMIN,
      E,
      EMAX = 0,
      EMIN = 0,
      EPS,
      OLDEMN,
      QMIN = 0,
      S,
      SAFMIN,
      T,
      TEMP,
      TOL,
      TOL2,
      TRACE,
      ZMAX,
      TEMPE,
      TEMPQ;
  final IINFO = Box(0),
      N0 = Box(0),
      NFAIL = Box(0),
      ITER = Box(0),
      NDIV = Box(0),
      PP = Box(0),
      TTYPE = Box(0);

  final DMIN = Box(0.0),
      DESIG = Box(0.0),
      SIGMA = Box(0.0),
      QMAX = Box(0.0),
      DMIN1 = Box(0.0),
      DMIN2 = Box(0.0),
      DN = Box(0.0),
      DN1 = Box(0.0),
      DN2 = Box(0.0),
      G = Box(0.0),
      TAU = Box(0.0);

  // Test the input arguments.
  // (in case DLASQ2 is not called by DLASQ1)

  INFO.value = 0;
  EPS = dlamch('Precision');
  SAFMIN = dlamch('Safe minimum');
  TOL = EPS * HUNDRD;
  TOL2 = pow(TOL, 2).toDouble();

  if (N < 0) {
    INFO.value = -1;
    xerbla('DLASQ2', 1);
    return;
  } else if (N == 0) {
    return;
  } else if (N == 1) {
    // 1-by-1 case.

    if (Z[1] < ZERO) {
      INFO.value = -201;
      xerbla('DLASQ2', 2);
    }
    return;
  } else if (N == 2) {
    // 2-by-2 case.

    if (Z[1] < ZERO) {
      INFO.value = -201;
      xerbla('DLASQ2', 2);
      return;
    } else if (Z[2] < ZERO) {
      INFO.value = -202;
      xerbla('DLASQ2', 2);
      return;
    } else if (Z[3] < ZERO) {
      INFO.value = -203;
      xerbla('DLASQ2', 2);
      return;
    } else if (Z[3] > Z[1]) {
      D = Z[3];
      Z[3] = Z[1];
      Z[1] = D;
    }
    Z[5] = Z[1] + Z[2] + Z[3];
    if (Z[2] > Z[3] * TOL2) {
      T = HALF * ((Z[1] - Z[3]) + Z[2]);
      S = Z[3] * (Z[2] / T);
      if (S <= T) {
        S = Z[3] * (Z[2] / (T * (ONE + sqrt(ONE + S / T))));
      } else {
        S = Z[3] * (Z[2] / (T + sqrt(T) * sqrt(T + S)));
      }
      T = Z[1] + (S + Z[2]);
      Z[3] *= (Z[1] / T);
      Z[1] = T;
    }
    Z[2] = Z[3];
    Z[6] = Z[2] + Z[1];
    return;
  }

  // Check for negative data and compute sums of q's and e's.

  Z[2 * N] = ZERO;
  EMIN = Z[2];
  QMAX.value = ZERO;
  ZMAX = ZERO;
  D = ZERO;
  E = ZERO;

  for (K = 1; K <= 2 * (N - 1); K += 2) {
    if (Z[K] < ZERO) {
      INFO.value = -(200 + K);
      xerbla('DLASQ2', 2);
      return;
    } else if (Z[K + 1] < ZERO) {
      INFO.value = -(200 + K + 1);
      xerbla('DLASQ2', 2);
      return;
    }
    D += Z[K];
    E += Z[K + 1];
    QMAX.value = max(QMAX.value, Z[K]);
    EMIN = min(EMIN, Z[K + 1]);
    ZMAX = max(QMAX.value, max(ZMAX, Z[K + 1]));
  }
  if (Z[2 * N - 1] < ZERO) {
    INFO.value = -(200 + 2 * N - 1);
    xerbla('DLASQ2', 2);
    return;
  }
  D += Z[2 * N - 1];
  QMAX.value = max(QMAX.value, Z[2 * N - 1]);
  ZMAX = max(QMAX.value, ZMAX);

  // Check for diagonality.

  if (E == ZERO) {
    for (K = 2; K <= N; K++) {
      Z[K] = Z[2 * K - 1];
    }
    dlasrt('D', N, Z, IINFO);
    Z[2 * N - 1] = D;
    return;
  }

  TRACE = D + E;

  // Check for zero data.

  if (TRACE == ZERO) {
    Z[2 * N - 1] = ZERO;
    return;
  }

  // Check whether the machine is IEEE conformable.

  IEEE = ilaenv(10, 'DLASQ2', 'N', 1, 2, 3, 4) == 1;

  // Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).

  for (K = 2 * N; K >= 2; K -= 2) {
    Z[2 * K] = ZERO;
    Z[2 * K - 1] = Z[K];
    Z[2 * K - 2] = ZERO;
    Z[2 * K - 3] = Z[K - 1];
  }

  I0 = 1;
  N0.value = N;

  // Reverse the qd-array, if warranted.

  if (CBIAS * Z[4 * I0 - 3] < Z[4 * N0.value - 3]) {
    IPN4 = 4 * (I0 + N0.value);
    for (I4 = 4 * I0; I4 <= 2 * (I0 + N0.value - 1); I4 += 4) {
      TEMP = Z[I4 - 3];
      Z[I4 - 3] = Z[IPN4 - I4 - 3];
      Z[IPN4 - I4 - 3] = TEMP;
      TEMP = Z[I4 - 1];
      Z[I4 - 1] = Z[IPN4 - I4 - 5];
      Z[IPN4 - I4 - 5] = TEMP;
    }
  }

  // Initial split checking via dqd and Li's test.

  PP.value = 0;

  for (K = 1; K <= 2; K++) {
    D = Z[4 * N0.value + PP.value - 3];
    for (I4 = 4 * (N0.value - 1) + PP.value; I4 >= 4 * I0 + PP.value; I4 -= 4) {
      if (Z[I4 - 1] <= TOL2 * D) {
        Z[I4 - 1] = -ZERO;
        D = Z[I4 - 3];
      } else {
        D = Z[I4 - 3] * (D / (D + Z[I4 - 1]));
      }
    }

    // dqd maps Z to ZZ plus Li's test.

    EMIN = Z[4 * I0 + PP.value + 1];
    D = Z[4 * I0 + PP.value - 3];
    for (I4 = 4 * I0 + PP.value; I4 <= 4 * (N0.value - 1) + PP.value; I4 += 4) {
      Z[I4 - 2 * PP.value - 2] = D + Z[I4 - 1];
      if (Z[I4 - 1] <= TOL2 * D) {
        Z[I4 - 1] = -ZERO;
        Z[I4 - 2 * PP.value - 2] = D;
        Z[I4 - 2 * PP.value] = ZERO;
        D = Z[I4 + 1];
      } else if (SAFMIN * Z[I4 + 1] < Z[I4 - 2 * PP.value - 2] &&
          SAFMIN * Z[I4 - 2 * PP.value - 2] < Z[I4 + 1]) {
        TEMP = Z[I4 + 1] / Z[I4 - 2 * PP.value - 2];
        Z[I4 - 2 * PP.value] = Z[I4 - 1] * TEMP;
        D *= TEMP;
      } else {
        Z[I4 - 2 * PP.value] =
            Z[I4 + 1] * (Z[I4 - 1] / Z[I4 - 2 * PP.value - 2]);
        D = Z[I4 + 1] * (D / Z[I4 - 2 * PP.value - 2]);
      }
      EMIN = min(EMIN, Z[I4 - 2 * PP.value]);
    }
    Z[4 * N0.value - PP.value - 2] = D;

    // Now find qmax.

    QMAX.value = Z[4 * I0 - PP.value - 2];
    for (I4 = 4 * I0 - PP.value + 2;
        4 < 0
            ? I4 >= 4 * N0.value - PP.value - 2
            : I4 <= 4 * N0.value - PP.value - 2;
        I4 += 4) {
      QMAX.value = max(QMAX.value, Z[I4]);
    }

    // Prepare for the next iteration on K.

    PP.value = 1 - PP.value;
  }

  // Initialise variables to pass to DLASQ3.

  TTYPE.value = 0;
  DMIN1.value = ZERO;
  DMIN2.value = ZERO;
  DN.value = ZERO;
  DN1.value = ZERO;
  DN2.value = ZERO;
  G.value = ZERO;
  TAU.value = ZERO;

  ITER.value = 2;
  NFAIL.value = 0;
  NDIV.value = 2 * (N0.value - I0);
  var flag = false;
  iwhilaLoop:
  for (IWHILA = 1; IWHILA <= N + 1; IWHILA++) {
    if (N0.value < 1) {
      flag = true;
      break;
    }

    // While array unfinished do

    // E(N0) holds the value of SIGMA when submatrix in I0:N0
    // splits from the rest of the array, but is negated.

    DESIG.value = ZERO;
    if (N0.value == N) {
      SIGMA.value = ZERO;
    } else {
      SIGMA.value = -Z[4 * N0.value - 1];
    }
    if (SIGMA.value < ZERO) {
      INFO.value = 1;
      return;
    }

    // Find last unreduced submatrix's top index I0, find QMAX and
    // EMIN. Find Gershgorin-type bound if Q's much greater than E's.

    EMAX = ZERO;
    if (N0.value > I0) {
      EMIN = Z[4 * N0.value - 5].abs();
    } else {
      EMIN = ZERO;
    }
    QMIN = Z[4 * N0.value - 3];
    QMAX.value = QMIN;
    var topIndexFound = false;
    for (I4 = 4 * N0.value; I4 >= 8; I4 -= 4) {
      if (Z[I4 - 5] <= ZERO) {
        topIndexFound = true;
        break;
      }
      if (QMIN >= FOUR * EMAX) {
        QMIN = min(QMIN, Z[I4 - 3]);
        EMAX = max(EMAX, Z[I4 - 5]);
      }
      QMAX.value = max(QMAX.value, Z[I4 - 7] + Z[I4 - 5]);
      EMIN = min(EMIN, Z[I4 - 5]);
    }
    if (!topIndexFound) {
      I4 = 4;
    }

    I0 = I4 ~/ 4;
    PP.value = 0;

    if (N0.value - I0 > 1) {
      DEE = Z[4 * I0 - 3];
      DEEMIN = DEE;
      KMIN = I0;
      for (I4 = 4 * I0 + 1; I4 <= 4 * N0.value - 3; I4 += 4) {
        DEE = Z[I4] * (DEE / (DEE + Z[I4 - 2]));
        if (DEE <= DEEMIN) {
          DEEMIN = DEE;
          KMIN = (I4 + 3) ~/ 4;
        }
      }
      if ((KMIN - I0) * 2 < N0.value - KMIN &&
          DEEMIN <= HALF * Z[4 * N0.value - 3]) {
        IPN4 = 4 * (I0 + N0.value);
        PP.value = 2;
        for (I4 = 4 * I0; I4 <= 2 * (I0 + N0.value - 1); I4 += 4) {
          TEMP = Z[I4 - 3];
          Z[I4 - 3] = Z[IPN4 - I4 - 3];
          Z[IPN4 - I4 - 3] = TEMP;
          TEMP = Z[I4 - 2];
          Z[I4 - 2] = Z[IPN4 - I4 - 2];
          Z[IPN4 - I4 - 2] = TEMP;
          TEMP = Z[I4 - 1];
          Z[I4 - 1] = Z[IPN4 - I4 - 5];
          Z[IPN4 - I4 - 5] = TEMP;
          TEMP = Z[I4];
          Z[I4] = Z[IPN4 - I4 - 4];
          Z[IPN4 - I4 - 4] = TEMP;
        }
      }
    }

    // Put -(initial shift) into DMIN.

    DMIN.value = -max(ZERO, QMIN - TWO * sqrt(QMIN) * sqrt(EMAX));

    // Now I0:N0 is unreduced.
    // PP = 0 for ping, PP = 1 for pong.
    // PP = 2 indicates that flipping was applied to the Z array and
    //        and that the tests for deflation upon entry in DLASQ3
    //        should not be performed.

    NBIG = 100 * (N0.value - I0 + 1);
    for (IWHILB = 1; IWHILB <= NBIG; IWHILB++) {
      if (I0 > N0.value) continue iwhilaLoop;

      // While submatrix unfinished take a good dqds step.

      dlasq3(I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, ITER, NDIV, IEEE,
          TTYPE, DMIN1, DMIN2, DN, DN1, DN2, G, TAU);

      PP.value = 1 - PP.value;

      // When EMIN is very small check for splits.

      if (PP.value == 0 && N0.value - I0 >= 3) {
        if (Z[4 * N0.value] <= TOL2 * QMAX.value ||
            Z[4 * N0.value - 1] <= TOL2 * SIGMA.value) {
          SPLT = I0 - 1;
          QMAX.value = Z[4 * I0 - 3];
          EMIN = Z[4 * I0 - 1];
          OLDEMN = Z[4 * I0];
          for (I4 = 4 * I0; I4 <= 4 * (N0.value - 3); I4 += 4) {
            if (Z[I4] <= TOL2 * Z[I4 - 3] || Z[I4 - 1] <= TOL2 * SIGMA.value) {
              Z[I4 - 1] = -SIGMA.value;
              SPLT = I4 ~/ 4;
              QMAX.value = ZERO;
              EMIN = Z[I4 + 3];
              OLDEMN = Z[I4 + 4];
            } else {
              QMAX.value = max(QMAX.value, Z[I4 + 1]);
              EMIN = min(EMIN, Z[I4 - 1]);
              OLDEMN = min(OLDEMN, Z[I4]);
            }
          }
          Z[4 * N0.value - 1] = EMIN;
          Z[4 * N0.value] = OLDEMN;
          I0 = SPLT + 1;
        }
      }
    }

    INFO.value = 2;

    // Maximum number of iterations exceeded, restore the shift
    // SIGMA and place the new d's and e's in a qd array.
    // This might need to be done for several blocks

    I1 = I0;
    N1 = N0.value;
    while (true) {
      TEMPQ = Z[4 * I0 - 3];
      Z[4 * I0 - 3] += SIGMA.value;
      for (K = I0 + 1; K <= N0.value; K++) {
        TEMPE = Z[4 * K - 5];
        Z[4 * K - 5] *= (TEMPQ / Z[4 * K - 7]);
        TEMPQ = Z[4 * K - 3];
        Z[4 * K - 3] += SIGMA.value + TEMPE - Z[4 * K - 5];
      }

      // Prepare to do this on the previous block if there is one

      if (I1 > 1) {
        N1 = I1 - 1;
        while ((I1 >= 2) && (Z[4 * I1 - 5] >= ZERO)) {
          I1--;
        }
        SIGMA.value = -Z[4 * N1 - 1];
        continue;
      }
      break;
    }

    for (K = 1; K <= N; K++) {
      Z[2 * K - 1] = Z[4 * K - 3];

      // Only the block 1..N0 is unfinished.  The rest of the e's
      // must be essentially zero, although sometimes other data
      // has been stored in them.

      if (K < N0.value) {
        Z[2 * K] = Z[4 * K - 1];
      } else {
        Z[2 * K] = 0;
      }
    }
    return;

    // end IWHILB
  }

  if (!flag) {
    INFO.value = 3;
    return;
  }

  // end IWHILA

  // Move q's to the front.

  for (K = 2; K <= N; K++) {
    Z[K] = Z[4 * K - 3];
  }

  // Sort and compute sum of eigenvalues.

  dlasrt('D', N, Z, IINFO);

  E = ZERO;
  for (K = N; K >= 1; K--) {
    E += Z[K];
  }

  // Store trace, sum(eigenvalues) and information on performance.

  Z[2 * N + 1] = TRACE;
  Z[2 * N + 2] = E;
  Z[2 * N + 3] = ITER.value.toDouble();
  Z[2 * N + 4] = NDIV.value / pow(N, 2);
  Z[2 * N + 5] = HUNDRD * NFAIL.value / ITER.value;
}
