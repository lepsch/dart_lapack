// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlaset.dart';

import '../matgen/dlaran.dart';
import '../matgen/zlarnd.dart';

void zlatm4(
  final int ITYPE,
  final int N,
  final int NZ1,
  final int NZ2,
  final bool RSIGN,
  final double AMAGN,
  final double RCOND,
  final double TRIANG,
  final int IDIST,
  final Array<int> ISEED_,
  final Matrix<Complex> A_,
  final int LDA,
) {
  final A = A_.having(ld: LDA);
  final ISEED = ISEED_.having(length: 4);
  const ZERO = 0.0, ONE = 1.0;
  int I, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN;
  double ALPHA;
  Complex CTEMP;

  if (N <= 0) return;
  zlaset('Full', N, N, Complex.zero, Complex.zero, A, LDA);

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
    //  GO TO ( 10, 30, 50, 80, 100, 120, 140, 160, 180, 200 )
    switch (ITYPE.abs()) {
      // abs(ITYPE) = 1: Identity
      case 1:
        for (JD = 1; JD <= N; JD++) {
          A[JD][JD] = Complex.one;
        }
        break;

      // abs(ITYPE) = 2: Transposed Jordan block
      case 2:
        for (JD = 1; JD <= N - 1; JD++) {
          A[JD + 1][JD] = Complex.one;
        }
        ISDB = 1;
        ISDE = N - 1;
        break;

      // abs(ITYPE) = 3: Transposed Jordan block, followed by the
      //                 identity.
      case 3:
        K = (N - 1) ~/ 2;
        for (JD = 1; JD <= K; JD++) {
          A[JD + 1][JD] = Complex.one;
        }
        ISDB = 1;
        ISDE = K;
        for (JD = K + 2; JD <= 2 * K + 1; JD++) {
          A[JD][JD] = Complex.one;
        }
        break;

      // abs(ITYPE) = 4: 1,...,k
      case 4:
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = (JD - NZ1).toComplex();
        }
        break;

      // abs(ITYPE) = 5: One large D value:
      case 5:
        for (JD = KBEG + 1; JD <= KEND; JD++) {
          A[JD][JD] = RCOND.toComplex();
        }
        A[KBEG][KBEG] = Complex.one;
        break;

      // abs(ITYPE) = 6: One small D value:
      case 6:
        for (JD = KBEG; JD <= KEND - 1; JD++) {
          A[JD][JD] = Complex.one;
        }
        A[KEND][KEND] = RCOND.toComplex();
        break;

      // abs(ITYPE) = 7: Exponentially distributed D values:
      case 7:
        A[KBEG][KBEG] = Complex.one;
        if (KLEN > 1) {
          ALPHA = pow(RCOND, ONE / (KLEN - 1)).toDouble();
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] = pow(ALPHA, I - 1).toComplex();
          }
        }
        break;

      // abs(ITYPE) = 8: Arithmetically distributed D values:
      case 8:
        A[KBEG][KBEG] = Complex.one;
        if (KLEN > 1) {
          ALPHA = (ONE - RCOND) / (KLEN - 1);
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] = ((KLEN - I) * ALPHA + RCOND).toComplex();
          }
        }
        break;

      // abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):
      case 9:
        ALPHA = log(RCOND);
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = exp(ALPHA * dlaran(ISEED)).toComplex();
        }
        break;

      // abs(ITYPE) = 10: Randomly distributed D values from DIST
      case 10:
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = zlarnd(IDIST, ISEED);
        }
        break;
    }

    // Scale by AMAGN
    for (JD = KBEG; JD <= KEND; JD++) {
      A[JD][JD] = (AMAGN * A[JD][JD].real).toComplex();
    }
    for (JD = ISDB; JD <= ISDE; JD++) {
      A[JD + 1][JD] = AMAGN.toComplex() * A[JD + 1][JD].real.toComplex();
    }

    // If RSIGN = true , assign random signs to diagonal and
    // subdiagonal
    if (RSIGN) {
      for (JD = KBEG; JD <= KEND; JD++) {
        if (A[JD][JD].real != ZERO) {
          CTEMP = zlarnd(3, ISEED);
          CTEMP /= CTEMP.abs().toComplex();
          A[JD][JD] = CTEMP * A[JD][JD].real.toComplex();
        }
      }
      for (JD = ISDB; JD <= ISDE; JD++) {
        if (A[JD + 1][JD].real != ZERO) {
          CTEMP = zlarnd(3, ISEED);
          CTEMP /= CTEMP.abs().toComplex();
          A[JD + 1][JD] = CTEMP * A[JD + 1][JD].real.toComplex();
        }
      }
    }

    // Reverse if ITYPE < 0
    if (ITYPE < 0) {
      for (JD = KBEG; JD <= (KBEG + KEND - 1) ~/ 2; JD++) {
        CTEMP = A[JD][JD];
        A[JD][JD] = A[KBEG + KEND - JD][KBEG + KEND - JD];
        A[KBEG + KEND - JD][KBEG + KEND - JD] = CTEMP;
      }
      for (JD = 1; JD <= (N - 1) ~/ 2; JD++) {
        CTEMP = A[JD + 1][JD];
        A[JD + 1][JD] = A[N + 1 - JD][N - JD];
        A[N + 1 - JD][N - JD] = CTEMP;
      }
    }
  }

  // Fill in upper triangle
  if (TRIANG != ZERO) {
    for (JC = 2; JC <= N; JC++) {
      for (JR = 1; JR <= JC - 1; JR++) {
        A[JR][JC] = TRIANG.toComplex() * zlarnd(IDIST, ISEED);
      }
    }
  }
}
