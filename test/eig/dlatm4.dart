import 'dart:math';

import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

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
  final Array<int> ISEED,
  final Matrix<double> A,
  final int LDA,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const HALF = ONE / TWO;
  int I, IOFF, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN;
  double ALPHA, CL, CR, SAFMIN, SL, SR, SV1, SV2, TEMP;

  if (N <= 0) return;
  dlaset('Full', N, N, ZERO, ZERO, A, LDA);

  // Insure a correct ISEED

  if ((ISEED[4] % 2) != 1) ISEED[4] = ISEED[4] + 1;

  // Compute diagonal and subdiagonal according to ITYPE, NZ1, NZ2,
  // and RCOND

  if (ITYPE != 0) {
    if ((ITYPE).abs() >= 4) {
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
    //  GO TO ( 10, 30, 50, 80, 100, 120, 140, 160, 180, 200 )( ITYPE ).abs();
    switch (ITYPE.abs()) {
      case 1:
        // abs(ITYPE) = 1: Identity
        for (JD = 1; JD <= N; JD++) {
          A[JD][JD] = ONE;
        }
        break;

      case 2:
        // abs(ITYPE) = 2: Transposed Jordan block
        for (JD = 1; JD <= N - 1; JD++) {
          A[JD + 1][JD] = ONE;
        }
        ISDB = 1;
        ISDE = N - 1;
        break;

      case 3:
        // abs(ITYPE) = 3: Transposed Jordan block, followed by the
        //                 identity.

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

      case 4:
        // abs(ITYPE) = 4: 1,...,k

        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = (JD - NZ1).toDouble();
        }
        break;

      case 5:
        // abs(ITYPE) = 5: One large D value:

        for (JD = KBEG + 1; JD <= KEND; JD++) {
          A[JD][JD] = RCOND;
        }
        A[KBEG][KBEG] = ONE;
        break;

      case 6:
        // abs(ITYPE) = 6: One small D value:

        for (JD = KBEG; JD <= KEND - 1; JD++) {
          A[JD][JD] = ONE;
        }
        A[KEND][KEND] = RCOND;
        break;

      case 7:
        // abs(ITYPE) = 7: Exponentially distributed D values:

        A[KBEG][KBEG] = ONE;
        if (KLEN > 1) {
          ALPHA = pow(RCOND, (ONE / (KLEN - 1))).toDouble();
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] = pow(ALPHA, (I - 1)).toDouble();
          }
        }
        break;

      case 8:
        // abs(ITYPE) = 8: Arithmetically distributed D values:

        A[KBEG][KBEG] = ONE;
        if (KLEN > 1) {
          ALPHA = (ONE - RCOND) / (KLEN - 1).toDouble();
          for (I = 2; I <= KLEN; I++) {
            A[NZ1 + I][NZ1 + I] = (KLEN - I).toDouble() * ALPHA + RCOND;
          }
        }
        break;

      case 9:
        // abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):

        ALPHA = log(RCOND);
        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = exp(ALPHA * dlaran(ISEED));
        }
        break;

      case 10:
        // abs(ITYPE) = 10: Randomly distributed D values from DIST

        for (JD = KBEG; JD <= KEND; JD++) {
          A[JD][JD] = dlarnd(IDIST, ISEED);
        }
        break;
    }

    // Scale by AMAGN

    for (JD = KBEG; JD <= KEND; JD++) {
      A[JD][JD] = AMAGN * (A[JD][JD]).toDouble();
    }
    for (JD = ISDB; JD <= ISDE; JD++) {
      A[JD + 1][JD] = AMAGN * (A[JD + 1][JD]).toDouble();
    }

    // If ISIGN = 1 or 2, assign random signs to diagonal and
    // subdiagonal

    if (ISIGN > 0) {
      for (JD = KBEG; JD <= KEND; JD++) {
        if ((A[JD][JD]).toDouble() != ZERO) {
          if (dlaran(ISEED) > HALF) A[JD][JD] = -A[JD][JD];
        }
      }
      for (JD = ISDB; JD <= ISDE; JD++) {
        if ((A[JD + 1][JD]).toDouble() != ZERO) {
          if (dlaran(ISEED) > HALF) A[JD + 1][JD] = -A[JD + 1][JD];
        }
      }
    }

    // Reverse if ITYPE < 0

    if (ITYPE < 0) {
      for (JD = KBEG; JD <= (KBEG + KEND - 1) / 2; JD++) {
        TEMP = A[JD][JD];
        A[JD][JD] = A[KBEG + KEND - JD][KBEG + KEND - JD];
        A[KBEG + KEND - JD][KBEG + KEND - JD] = TEMP;
      }
      for (JD = 1; JD <= (N - 1) / 2; JD++) {
        TEMP = A[JD + 1][JD];
        A[JD + 1][JD] = A[N + 1 - JD][N - JD];
        A[N + 1 - JD][N - JD] = TEMP;
      }
    }

    // If ISIGN = 2, and no subdiagonals already, then apply
    // random rotations to make 2x2 blocks.

    if (ISIGN == 2 && ITYPE != 2 && ITYPE != 3) {
      SAFMIN = dlamch('S');
      for (JD = KBEG; 2 < 0 ? JD >= KEND - 1 : JD <= KEND - 1; JD += 2) {
        if (dlaran(ISEED) > HALF) {
          // Rotation on left.

          CL = TWO * dlaran(ISEED) - ONE;
          SL = TWO * dlaran(ISEED) - ONE;
          TEMP = ONE / max(SAFMIN, sqrt(pow(CL, 2) + pow(SL, 2)));
          CL = CL * TEMP;
          SL = SL * TEMP;

          // Rotation on right.

          CR = TWO * dlaran(ISEED) - ONE;
          SR = TWO * dlaran(ISEED) - ONE;
          TEMP = ONE / max(SAFMIN, sqrt(pow(CR, 2) + pow(SR, 2)));
          CR = CR * TEMP;
          SR = SR * TEMP;

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
