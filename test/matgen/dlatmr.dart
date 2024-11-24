// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlangb.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlansb.dart';
import 'package:dart_lapack/src/dlansp.dart';
import 'package:dart_lapack/src/dlansy.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

import 'dlatm1.dart';
import 'dlatm2.dart';
import 'dlatm3.dart';

void dlatmr(
  final int M,
  final int N,
  final String DIST,
  final Array<int> ISEED_,
  final String SYM,
  final Array<double> D_,
  final int MODE,
  final double COND,
  final double DMAX,
  final String RSIGN,
  final String GRADE,
  final Array<double> DL_,
  final int MODEL,
  final double CONDL,
  final Array<double> DR_,
  final int MODER,
  final double CONDR,
  final String PIVTNG,
  final Array<int> IPIVOT_,
  final int KL,
  final int KU,
  final double SPARSE,
  final double ANORM,
  final String PACK,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having();
  final D = D_.having();
  final DL = DL_.having();
  final DR = DR_.having();
  final IPIVOT = IPIVOT_.having();
  final A = A_.having(ld: LDA);
  final IWORK = IWORK_.having();
  const ZERO = 0.0;
  const ONE = 1.0;
  bool BADPVT, DZERO, FULBND;
  int I,
      IDIST,
      IGRADE,
      IISUB = 0,
      IPACK,
      IPVTNG,
      IRSIGN,
      ISYM,
      J,
      JJSUB,
      K,
      KLL,
      KUU,
      MNMIN,
      MNSUB,
      MXSUB,
      NPVTS = 0;
  double ALPHA, ONORM = 0, TEMP;
  final TEMPA = Array<double>(1);
  final ISUB = Box(0), JSUB = Box(0);

  // 1)      Decode and Test the input parameters.
  // Initialize flags & seed.

  INFO.value = 0;

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Decode DIST

  if (lsame(DIST, 'U')) {
    IDIST = 1;
  } else if (lsame(DIST, 'S')) {
    IDIST = 2;
  } else if (lsame(DIST, 'N')) {
    IDIST = 3;
  } else {
    IDIST = -1;
  }

  // Decode SYM

  if (lsame(SYM, 'S')) {
    ISYM = 0;
  } else if (lsame(SYM, 'N')) {
    ISYM = 1;
  } else if (lsame(SYM, 'H')) {
    ISYM = 0;
  } else {
    ISYM = -1;
  }

  // Decode RSIGN

  if (lsame(RSIGN, 'F')) {
    IRSIGN = 0;
  } else if (lsame(RSIGN, 'T')) {
    IRSIGN = 1;
  } else {
    IRSIGN = -1;
  }

  // Decode PIVTNG

  if (lsame(PIVTNG, 'N')) {
    IPVTNG = 0;
  } else if (lsame(PIVTNG, ' ')) {
    IPVTNG = 0;
  } else if (lsame(PIVTNG, 'L')) {
    IPVTNG = 1;
    NPVTS = M;
  } else if (lsame(PIVTNG, 'R')) {
    IPVTNG = 2;
    NPVTS = N;
  } else if (lsame(PIVTNG, 'B')) {
    IPVTNG = 3;
    NPVTS = min(N, M);
  } else if (lsame(PIVTNG, 'F')) {
    IPVTNG = 3;
    NPVTS = min(N, M);
  } else {
    IPVTNG = -1;
  }

  // Decode GRADE

  if (lsame(GRADE, 'N')) {
    IGRADE = 0;
  } else if (lsame(GRADE, 'L')) {
    IGRADE = 1;
  } else if (lsame(GRADE, 'R')) {
    IGRADE = 2;
  } else if (lsame(GRADE, 'B')) {
    IGRADE = 3;
  } else if (lsame(GRADE, 'E')) {
    IGRADE = 4;
  } else if (lsame(GRADE, 'H') || lsame(GRADE, 'S')) {
    IGRADE = 5;
  } else {
    IGRADE = -1;
  }

  // Decode PACK

  if (lsame(PACK, 'N')) {
    IPACK = 0;
  } else if (lsame(PACK, 'U')) {
    IPACK = 1;
  } else if (lsame(PACK, 'L')) {
    IPACK = 2;
  } else if (lsame(PACK, 'C')) {
    IPACK = 3;
  } else if (lsame(PACK, 'R')) {
    IPACK = 4;
  } else if (lsame(PACK, 'B')) {
    IPACK = 5;
  } else if (lsame(PACK, 'Q')) {
    IPACK = 6;
  } else if (lsame(PACK, 'Z')) {
    IPACK = 7;
  } else {
    IPACK = -1;
  }

  // Set certain internal parameters

  MNMIN = min(M, N);
  KLL = min(KL, M - 1);
  KUU = min(KU, N - 1);

  // If inv(DL) is used, check to see if DL has a zero entry.

  DZERO = false;
  if (IGRADE == 4 && MODEL == 0) {
    for (I = 1; I <= M; I++) {
      if (DL[I] == ZERO) DZERO = true;
    }
  }

  // Check values in IPIVOT

  BADPVT = false;
  if (IPVTNG > 0) {
    for (J = 1; J <= NPVTS; J++) {
      if (IPIVOT[J] <= 0 || IPIVOT[J] > NPVTS) BADPVT = true;
    }
  }

  // Set INFO if an error

  if (M < 0) {
    INFO.value = -1;
  } else if (M != N && ISYM == 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (IDIST == -1) {
    INFO.value = -3;
  } else if (ISYM == -1) {
    INFO.value = -5;
  } else if (MODE < -6 || MODE > 6) {
    INFO.value = -7;
  } else if ((MODE != -6 && MODE != 0 && MODE != 6) && COND < ONE) {
    INFO.value = -8;
  } else if ((MODE != -6 && MODE != 0 && MODE != 6) && IRSIGN == -1) {
    INFO.value = -10;
  } else if (IGRADE == -1 ||
      (IGRADE == 4 && M != N) ||
      ((IGRADE >= 1 && IGRADE <= 4) && ISYM == 0)) {
    INFO.value = -11;
  } else if (IGRADE == 4 && DZERO) {
    INFO.value = -12;
  } else if ((IGRADE == 1 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5) &&
      (MODEL < -6 || MODEL > 6)) {
    INFO.value = -13;
  } else if ((IGRADE == 1 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5) &&
      (MODEL != -6 && MODEL != 0 && MODEL != 6) &&
      CONDL < ONE) {
    INFO.value = -14;
  } else if ((IGRADE == 2 || IGRADE == 3) && (MODER < -6 || MODER > 6)) {
    INFO.value = -16;
  } else if ((IGRADE == 2 || IGRADE == 3) &&
      (MODER != -6 && MODER != 0 && MODER != 6) &&
      CONDR < ONE) {
    INFO.value = -17;
  } else if (IPVTNG == -1 ||
      (IPVTNG == 3 && M != N) ||
      ((IPVTNG == 1 || IPVTNG == 2) && ISYM == 0)) {
    INFO.value = -18;
  } else if (IPVTNG != 0 && BADPVT) {
    INFO.value = -19;
  } else if (KL < 0) {
    INFO.value = -20;
  } else if (KU < 0 || (ISYM == 0 && KL != KU)) {
    INFO.value = -21;
  } else if (SPARSE < ZERO || SPARSE > ONE) {
    INFO.value = -22;
  } else if (IPACK == -1 ||
      ((IPACK == 1 || IPACK == 2 || IPACK == 5 || IPACK == 6) && ISYM == 1) ||
      (IPACK == 3 && ISYM == 1 && (KL != 0 || M != N)) ||
      (IPACK == 4 && ISYM == 1 && (KU != 0 || M != N))) {
    INFO.value = -24;
  } else if (((IPACK == 0 || IPACK == 1 || IPACK == 2) && LDA < max(1, M)) ||
      ((IPACK == 3 || IPACK == 4) && LDA < 1) ||
      ((IPACK == 5 || IPACK == 6) && LDA < KUU + 1) ||
      (IPACK == 7 && LDA < KLL + KUU + 1)) {
    INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('DLATMR', -INFO.value);
    return;
  }

  // Decide if we can pivot consistently

  FULBND = false;
  if (KUU == N - 1 && KLL == M - 1) FULBND = true;

  // Initialize random number generator

  for (I = 1; I <= 4; I++) {
    ISEED[I] = (ISEED[I].abs() % 4096);
  }

  ISEED[4] = 2 * (ISEED[4] ~/ 2) + 1;

  // 2)      Set up D, DL, and DR, if indicated.

  // Compute D according to COND and MODE

  dlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO);
  if (INFO.value != 0) {
    INFO.value = 1;
    return;
  }
  if (MODE != 0 && MODE != -6 && MODE != 6) {
    // Scale by DMAX

    TEMP = D[1].abs();
    for (I = 2; I <= MNMIN; I++) {
      TEMP = max(TEMP, D[I].abs());
    }
    if (TEMP == ZERO && DMAX != ZERO) {
      INFO.value = 2;
      return;
    }
    if (TEMP != ZERO) {
      ALPHA = DMAX / TEMP;
    } else {
      ALPHA = ONE;
    }
    for (I = 1; I <= MNMIN; I++) {
      D[I] = ALPHA * D[I];
    }
  }

  // Compute DL if grading set

  if (IGRADE == 1 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5) {
    dlatm1(MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO);
    if (INFO.value != 0) {
      INFO.value = 3;
      return;
    }
  }

  // Compute DR if grading set

  if (IGRADE == 2 || IGRADE == 3) {
    dlatm1(MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO);
    if (INFO.value != 0) {
      INFO.value = 4;
      return;
    }
  }

  // 3)     Generate IWORK if pivoting

  if (IPVTNG > 0) {
    for (I = 1; I <= NPVTS; I++) {
      IWORK[I] = I;
    }
    if (FULBND) {
      for (I = 1; I <= NPVTS; I++) {
        K = IPIVOT[I];
        J = IWORK[I];
        IWORK[I] = IWORK[K];
        IWORK[K] = J;
      }
    } else {
      for (I = NPVTS; I >= 1; I--) {
        K = IPIVOT[I];
        J = IWORK[I];
        IWORK[I] = IWORK[K];
        IWORK[K] = J;
      }
    }
  }

  // 4)      Generate matrices for each kind of PACKing
  // Always sweep matrix columnwise (if symmetric, upper
  // half only) so that matrix generated does not depend
  // on PACK

  if (FULBND) {
    // Use DLATM3 so matrices generated with differing PIVOTing only
    // differ only in the order of their rows and/or columns.

    if (IPACK == 0) {
      if (ISYM == 0) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            A[ISUB.value][JSUB.value] = TEMP;
            A[JSUB.value][ISUB.value] = TEMP;
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= M; I++) {
            TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            A[ISUB.value][JSUB.value] = TEMP;
          }
        }
      }
    } else if (IPACK == 1) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE,
              DL, DR, IPVTNG, IWORK, SPARSE);
          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          A[MNSUB][MXSUB] = TEMP;
          if (MNSUB != MXSUB) A[MXSUB][MNSUB] = ZERO;
        }
      }
    } else if (IPACK == 2) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE,
              DL, DR, IPVTNG, IWORK, SPARSE);
          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          A[MXSUB][MNSUB] = TEMP;
          if (MNSUB != MXSUB) A[MNSUB][MXSUB] = ZERO;
        }
      }
    } else if (IPACK == 3) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE,
              DL, DR, IPVTNG, IWORK, SPARSE);

          // Compute K = location of (ISUB,JSUB) entry in packed
          // array

          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          K = MXSUB * (MXSUB - 1) ~/ 2 + MNSUB;

          // Convert K to (IISUB,JJSUB) location

          JJSUB = (K - 1) ~/ LDA + 1;
          IISUB = K - LDA * (JJSUB - 1);

          A[IISUB][JJSUB] = TEMP;
        }
      }
    } else if (IPACK == 4) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE,
              DL, DR, IPVTNG, IWORK, SPARSE);

          // Compute K = location of (I,J) entry in packed array

          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          if (MNSUB == 1) {
            K = MXSUB;
          } else {
            K = N * (N + 1) ~/ 2 -
                (N - MNSUB + 1) * (N - MNSUB + 2) ~/ 2 +
                MXSUB -
                MNSUB +
                1;
          }

          // Convert K to (IISUB,JJSUB) location

          JJSUB = (K - 1) ~/ LDA + 1;
          IISUB = K - LDA * (JJSUB - 1);

          A[IISUB][JJSUB] = TEMP;
        }
      }
    } else if (IPACK == 5) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          if (I < 1) {
            A[J - I + 1][I + N] = ZERO;
          } else {
            TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            MNSUB = min(ISUB.value, JSUB.value);
            MXSUB = max(ISUB.value, JSUB.value);
            A[MXSUB - MNSUB + 1][MNSUB] = TEMP;
          }
        }
      }
    } else if (IPACK == 6) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE,
              DL, DR, IPVTNG, IWORK, SPARSE);
          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          A[MNSUB - MXSUB + KUU + 1][MXSUB] = TEMP;
        }
      }
    } else if (IPACK == 7) {
      if (ISYM == 0) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J; I++) {
            TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            MNSUB = min(ISUB.value, JSUB.value);
            MXSUB = max(ISUB.value, JSUB.value);
            A[MNSUB - MXSUB + KUU + 1][MXSUB] = TEMP;
            if (I < 1) A[J - I + 1 + KUU][I + N] = ZERO;
            if (I >= 1 && MNSUB != MXSUB) {
              A[MXSUB - MNSUB + 1 + KUU][MNSUB] = TEMP;
            }
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J + KLL; I++) {
            TEMP = dlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            A[ISUB.value - JSUB.value + KUU + 1][JSUB.value] = TEMP;
          }
        }
      }
    }
  } else {
    // Use DLATM2

    if (IPACK == 0) {
      if (ISYM == 0) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            A[I][J] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL,
                DR, IPVTNG, IWORK, SPARSE);
            A[J][I] = A[I][J];
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= M; I++) {
            A[I][J] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL,
                DR, IPVTNG, IWORK, SPARSE);
          }
        }
      }
    } else if (IPACK == 1) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          A[I][J] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR,
              IPVTNG, IWORK, SPARSE);
          if (I != J) A[J][I] = ZERO;
        }
      }
    } else if (IPACK == 2) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          A[J][I] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR,
              IPVTNG, IWORK, SPARSE);
          if (I != J) A[I][J] = ZERO;
        }
      }
    } else if (IPACK == 3) {
      ISUB.value = 0;
      JSUB.value = 1;
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          ISUB.value++;
          if (ISUB.value > LDA) {
            ISUB.value = 1;
            JSUB.value++;
          }
          A[ISUB.value][JSUB.value] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED,
              D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
        }
      }
    } else if (IPACK == 4) {
      if (ISYM == 0) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            // Compute K = location of (I,J) entry in packed array

            if (I == 1) {
              K = J;
            } else {
              K = N * (N + 1) ~/ 2 - (N - I + 1) * (N - I + 2) ~/ 2 + J - I + 1;
            }

            // Convert K to (ISUB,JSUB) location

            JSUB.value = (K - 1) ~/ LDA + 1;
            ISUB.value = K - LDA * (JSUB.value - 1);

            A[ISUB.value][JSUB.value] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED,
                D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          }
        }
      } else {
        ISUB.value = 0;
        JSUB.value = 1;
        for (J = 1; J <= N; J++) {
          for (I = J; I <= M; I++) {
            ISUB.value++;
            if (ISUB.value > LDA) {
              ISUB.value = 1;
              JSUB.value++;
            }
            A[ISUB.value][JSUB.value] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED,
                D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          }
        }
      }
    } else if (IPACK == 5) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          if (I < 1) {
            A[J - I + 1][I + N] = ZERO;
          } else {
            A[J - I + 1][I] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          }
        }
      }
    } else if (IPACK == 6) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          A[I - J + KUU + 1][J] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
              IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
        }
      }
    } else if (IPACK == 7) {
      if (ISYM == 0) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J; I++) {
            A[I - J + KUU + 1][J] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            if (I < 1) A[J - I + 1 + KUU][I + N] = ZERO;
            if (I >= 1 && I != J) A[J - I + 1 + KUU][I] = A[I - J + KUU + 1][J];
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J + KLL; I++) {
            A[I - J + KUU + 1][J] = dlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          }
        }
      }
    }
  }

  // 5)      Scaling the norm

  if (IPACK == 0) {
    ONORM = dlange('M', M, N, A, LDA, TEMPA);
  } else if (IPACK == 1) {
    ONORM = dlansy('M', 'U', N, A, LDA, TEMPA);
  } else if (IPACK == 2) {
    ONORM = dlansy('M', 'L', N, A, LDA, TEMPA);
  } else if (IPACK == 3) {
    ONORM = dlansp('M', 'U', N, A.asArray(), TEMPA);
  } else if (IPACK == 4) {
    ONORM = dlansp('M', 'L', N, A.asArray(), TEMPA);
  } else if (IPACK == 5) {
    ONORM = dlansb('M', 'L', N, KLL, A, LDA, TEMPA);
  } else if (IPACK == 6) {
    ONORM = dlansb('M', 'U', N, KUU, A, LDA, TEMPA);
  } else if (IPACK == 7) {
    ONORM = dlangb('M', N, KLL, KUU, A, LDA, TEMPA);
  }

  if (ANORM >= ZERO) {
    if (ANORM > ZERO && ONORM == ZERO) {
      // Desired scaling impossible

      INFO.value = 5;
      return;
    } else if ((ANORM > ONE && ONORM < ONE) || (ANORM < ONE && ONORM > ONE)) {
      // Scale carefully to avoid over / underflow

      if (IPACK <= 2) {
        for (J = 1; J <= N; J++) {
          dscal(M, ONE / ONORM, A(1, J).asArray(), 1);
          dscal(M, ANORM, A(1, J).asArray(), 1);
        }
      } else if (IPACK == 3 || IPACK == 4) {
        dscal(N * (N + 1) ~/ 2, ONE / ONORM, A.asArray(), 1);
        dscal(N * (N + 1) ~/ 2, ANORM, A.asArray(), 1);
      } else if (IPACK >= 5) {
        for (J = 1; J <= N; J++) {
          dscal(KLL + KUU + 1, ONE / ONORM, A(1, J).asArray(), 1);
          dscal(KLL + KUU + 1, ANORM, A(1, J).asArray(), 1);
        }
      }
    } else {
      // Scale straightforwardly

      if (IPACK <= 2) {
        for (J = 1; J <= N; J++) {
          dscal(M, ANORM / ONORM, A(1, J).asArray(), 1);
        }
      } else if (IPACK == 3 || IPACK == 4) {
        dscal(N * (N + 1) ~/ 2, ANORM / ONORM, A.asArray(), 1);
      } else if (IPACK >= 5) {
        for (J = 1; J <= N; J++) {
          dscal(KLL + KUU + 1, ANORM / ONORM, A(1, J).asArray(), 1);
        }
      }
    }
  }
}
