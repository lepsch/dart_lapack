// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

import '../matgen/dlaran.dart';
import '../matgen/dlarnd.dart';

void dlatm4(
  final int ITYPE,
  final int N,
  final int NZ1,
  final int NZ2,
  final int ISIGN,
  final double AMAGN,
  final double RCOND,
  final double TRIANG,
  final int IDIST,
  final Array<int> ISEED_,
  final Matrix<double> A_,
  final int LDA,
) {
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const HALF = ONE / TWO;
  int I, IOFF, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN;
  double ALPHA, CL, CR, SAFMIN, SL, SR, SV1, SV2, TEMP;

  if (N <= 0) return;
  dlaset('Full', N, N, ZERO, ZERO, A, LDA);

  // Insure a correct ISEED
  if ((ISEED[4] % 2) != 1) ISEED[4]++;

  // Compute diagonal and subdiagonal according to ITYPE, NZ1, NZ2,
  // and RCOND
  if (ITYPE != 0) {
    if (ITYPE.abs() >= 4) {
      KBEG = max(1, min(N, NZ1 + 1));
      KEND = max(KBEG, min(N, N - NZ2));
      KLEN = KEND + 1 - KBEG;
    } else {
      KBEG = 1;
      KEND = N;
      KLEN = N;
    }
    ISDB = 1;
    ISDE = 0;
    switch (ITYPE.abs()) {
      // abs(ITYPE) = 1: Identity
      case 1:
        for (JD = 1; JD <= N; JD++) {
          A[JD][JD] = ONE;
        }
        break;

      // abs(ITYPE) = 2: Transposed Jordan block
      case 2:
        for (JD = 1; JD <= N - 1; JD++) {
          A[JD + 1][JD] = ONE;
        }
        ISDB = 1;
        ISDE = N - 1;
        break;

      // abs(ITYPE) = 3: Transposed Jordan block, followed by the
      //                 identity.
      case 3:
        K = (N - 1) ~/ 2;
        for (JD = 1; JD <= K; JD++) {
          A[JD + 1][JD] = ONE;
        }
        ISDB = 1;
        ISDE = K;
        for (JD = K + 2; JD <= 2 * K + 1; JD++) {
          A[JD][JD] = ONE;
        }
        break;

      // abs(ITYPE) = 4: 1,...,k
      case 4:
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = (JD - NZ1).toDouble();
        }
        break;

      // abs(ITYPE) = 5: One large D value:
      case 5:
        for (JD = KBEG + 1; JD <= KEND; JD++) {
          A[JD][JD] = RCOND;
        }
        A[KBEG][KBEG] = ONE;
        break;

      // abs(ITYPE) = 6: One small D value:
      case 6:
        for (JD = KBEG; JD <= KEND - 1; JD++) {
          A[JD][JD] = ONE;
        }
        A[KEND][KEND] = RCOND;
        break;

      // abs(ITYPE) = 7: Exponentially distributed D values:
      case 7:
        A[KBEG][KBEG] = ONE;
        if (KLEN > 1) {
          ALPHA = pow(RCOND, (ONE / (KLEN - 1))).toDouble();
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] = pow(ALPHA, (I - 1)).toDouble();
          }
        }
        break;

      // abs(ITYPE) = 8: Arithmetically distributed D values:
      case 8:
        A[KBEG][KBEG] = ONE;
        if (KLEN > 1) {
          ALPHA = (ONE - RCOND) / (KLEN - 1);
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] = (KLEN - I) * ALPHA + RCOND;
          }
        }
        break;

      // abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):
      case 9:
        ALPHA = log(RCOND);
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = exp(ALPHA * dlaran(ISEED));
        }
        break;

      // abs(ITYPE) = 10: Randomly distributed D values from DIST
      case 10:
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = dlarnd(IDIST, ISEED);
        }
        break;
    }

    // Scale by AMAGN
    for (JD = KBEG; JD <= KEND; JD++) {
      A[JD][JD] = AMAGN * A[JD][JD];
    }
    for (JD = ISDB; JD <= ISDE; JD++) {
      A[JD + 1][JD] = AMAGN * A[JD + 1][JD];
    }

    // If ISIGN = 1 or 2, assign random signs to diagonal and
    // subdiagonal
    if (ISIGN > 0) {
      for (JD = KBEG; JD <= KEND; JD++) {
        if (A[JD][JD] != ZERO) {
          if (dlaran(ISEED) > HALF) A[JD][JD] = -A[JD][JD];
        }
      }
      for (JD = ISDB; JD <= ISDE; JD++) {
        if (A[JD + 1][JD] != ZERO) {
          if (dlaran(ISEED) > HALF) A[JD + 1][JD] = -A[JD + 1][JD];
        }
      }
    }

    // Reverse if ITYPE < 0
    if (ITYPE < 0) {
      for (JD = KBEG; JD <= (KBEG + KEND - 1) ~/ 2; JD++) {
        TEMP = A[JD][JD];
        A[JD][JD] = A[KBEG + KEND - JD][KBEG + KEND - JD];
        A[KBEG + KEND - JD][KBEG + KEND - JD] = TEMP;
      }
      for (JD = 1; JD <= (N - 1) ~/ 2; JD++) {
        TEMP = A[JD + 1][JD];
        A[JD + 1][JD] = A[N + 1 - JD][N - JD];
        A[N + 1 - JD][N - JD] = TEMP;
      }
    }

    // If ISIGN = 2, and no subdiagonals already, then apply
    // random rotations to make 2x2 blocks.
    if (ISIGN == 2 && ITYPE != 2 && ITYPE != 3) {
      SAFMIN = dlamch('S');
      for (JD = KBEG; JD <= KEND - 1; JD += 2) {
        if (dlaran(ISEED) > HALF) {
          // Rotation on left.
          CL = TWO * dlaran(ISEED) - ONE;
          SL = TWO * dlaran(ISEED) - ONE;
          TEMP = ONE / max(SAFMIN, sqrt(pow(CL, 2) + pow(SL, 2)));
          CL *= TEMP;
          SL *= TEMP;

          // Rotation on right.
          CR = TWO * dlaran(ISEED) - ONE;
          SR = TWO * dlaran(ISEED) - ONE;
          TEMP = ONE / max(SAFMIN, sqrt(pow(CR, 2) + pow(SR, 2)));
          CR *= TEMP;
          SR *= TEMP;

          // Apply
          SV1 = A[JD][JD];
          SV2 = A[JD + 1][JD + 1];
          A[JD][JD] = CL * CR * SV1 + SL * SR * SV2;
          A[JD + 1][JD] = -SL * CR * SV1 + CL * SR * SV2;
          A[JD][JD + 1] = -CL * SR * SV1 + SL * CR * SV2;
          A[JD + 1][JD + 1] = SL * SR * SV1 + CL * CR * SV2;
        }
      }
    }
  }

  // Fill in upper triangle (except for 2x2 blocks)
  if (TRIANG != ZERO) {
    if (ISIGN != 2 || ITYPE == 2 || ITYPE == 3) {
      IOFF = 1;
    } else {
      IOFF = 2;
      for (JR = 1; JR <= N - 1; JR++) {
        if (A[JR + 1][JR] == ZERO)
          // ignore: curly_braces_in_flow_control_structures
          A[JR][JR + 1] = TRIANG * dlarnd(IDIST, ISEED);
      }
    }

    for (JC = 2; JC <= N; JC++) {
      for (JR = 1; JR <= JC - IOFF; JR++) {
        A[JR][JC] = TRIANG * dlarnd(IDIST, ISEED);
      }
    }
  }
}
