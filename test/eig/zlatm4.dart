import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlaset.dart';

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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final ISEED = ISEED_.having(length: 4);
  const ZERO = 0.0, ONE = 1.0;
  int I, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN;
  double ALPHA;
  Complex CTEMP;
  // ..
  // .. External Functions ..
  //- double             DLARAN;
  //- Complex         zlarnd;
  // EXTERNAL DLARAN, zlarnd
  // ..
  // .. External Subroutines ..
  // EXTERNAL ZLASET
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC ABS, DBLE, DCMPLX, EXP, log, MAX, MIN, MOD

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
      case 1:
        // abs(ITYPE) = 1: Identity

        for (JD = 1; JD <= N; JD++) {
          A[JD][JD] = Complex.one;
        }
        break;

      case 2:

        // abs(ITYPE) = 2: Transposed Jordan block

        for (JD = 1; JD <= N - 1; JD++) {
          A[JD + 1][JD] = Complex.one;
        }
        ISDB = 1;
        ISDE = N - 1;
        break;

      case 3:

        // abs(ITYPE) = 3: Transposed Jordan block, followed by the
        //                 identity.

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

      case 4:

        // abs(ITYPE) = 4: 1,...,k

        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = (JD - NZ1).toComplex();
        }
        break;

      case 5:

        // abs(ITYPE) = 5: One large D value:

        for (JD = KBEG + 1; JD <= KEND; JD++) {
          A[JD][JD] = RCOND.toComplex();
        }
        A[KBEG][KBEG] = Complex.one;
        break;

      case 6:

        // abs(ITYPE) = 6: One small D value:

        for (JD = KBEG; JD <= KEND - 1; JD++) {
          A[JD][JD] = Complex.one;
        }
        A[KEND][KEND] = RCOND.toComplex();
        break;

      case 7:

        // abs(ITYPE) = 7: Exponentially distributed D values:

        A[KBEG][KBEG] = Complex.one;
        if (KLEN > 1) {
          ALPHA = pow(RCOND, ONE / (KLEN - 1)).toDouble();
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] = pow(ALPHA, (I - 1).toDouble()).toComplex();
          }
        }
        break;

      case 8:

        // abs(ITYPE) = 8: Arithmetically distributed D values:

        A[KBEG][KBEG] = Complex.one;
        if (KLEN > 1) {
          ALPHA = (ONE - RCOND) / (KLEN - 1).toDouble();
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] =
                ((KLEN - I).toDouble() * ALPHA + RCOND).toComplex();
          }
        }
        break;

      case 9:

        // abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):

        ALPHA = log(RCOND);
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = exp(ALPHA * dlaran(ISEED)).toComplex();
        }
        break;

      case 10:

        // abs(ITYPE) = 10: Randomly distributed D values from DIST

        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = zlarnd(IDIST, ISEED);
        }
        break;
    }

    // Scale by AMAGN

    for (JD = KBEG; JD <= KEND; JD++) {
      A[JD][JD] = (AMAGN * A[JD][JD].toDouble()).toComplex();
    }
    for (JD = ISDB; JD <= ISDE; JD++) {
      A[JD + 1][JD] = AMAGN.toComplex() * A[JD + 1][JD].real.toComplex();
    }

    // If RSIGN = true , assign random signs to diagonal and
    // subdiagonal

    if (RSIGN) {
      for (JD = KBEG; JD <= KEND; JD++) {
        if (A[JD][JD].toDouble() != ZERO) {
          CTEMP = zlarnd(3, ISEED);
          CTEMP /= CTEMP.abs().toComplex();
          A[JD][JD] = CTEMP * A[JD][JD].real.toComplex();
        }
      }
      for (JD = ISDB; JD <= ISDE; JD++) {
        if (A[JD + 1][JD].toDouble() != ZERO) {
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
