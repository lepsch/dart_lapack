// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlangb.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlansb.dart';
import 'package:dart_lapack/src/zlansp.dart';
import 'package:dart_lapack/src/zlansy.dart';

import 'zlatm1.dart';
import 'zlatm2.dart';
import 'zlatm3.dart';

void zlatmr(
  final int M,
  final int N,
  final String DIST,
  final Array<int> ISEED_,
  final String SYM,
  final Array<Complex> D_,
  final int MODE,
  final double COND,
  final Complex DMAX,
  final String RSIGN,
  final String GRADE,
  final Array<Complex> DL_,
  final int MODEL,
  final double CONDL,
  final Array<Complex> DR_,
  final int MODER,
  final double CONDR,
  final String PIVTNG,
  final Array<int> IPIVOT_,
  final int KL,
  final int KU,
  final double SPARSE,
  final double ANORM,
  final String PACK,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final ISEED = ISEED_.having();
  final D = D_.having();
  final DL = DL_.having();
  final DR = DR_.having();
  final A = A_.having(ld: LDA);
  final IPIVOT = IPIVOT_.having();
  final IWORK = IWORK_.having();

  const ZERO = 0.0;
  const ONE = 1.0;
  bool BADPVT, DZERO, FULBND;
  int I,
      IDIST,
      IGRADE,
      IISUB,
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
  double ONORM = 0, TEMP;
  Complex CALPHA, CTEMP;
  final TEMPA = Array<double>(1);
  final ISUB = Box(0), JSUB = Box(0);

  // 1)      Decode and Test the input parameters.
  //         Initialize flags & seed.

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
  } else if (lsame(DIST, 'D')) {
    IDIST = 4;
  } else {
    IDIST = -1;
  }

  // Decode SYM

  if (lsame(SYM, 'H')) {
    ISYM = 0;
  } else if (lsame(SYM, 'N')) {
    ISYM = 1;
  } else if (lsame(SYM, 'S')) {
    ISYM = 2;
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
  } else if (lsame(GRADE, 'H')) {
    IGRADE = 5;
  } else if (lsame(GRADE, 'S')) {
    IGRADE = 6;
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
      if (DL[I] == Complex.zero) DZERO = true;
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
  } else if (M != N && (ISYM == 0 || ISYM == 2)) {
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
      ((IGRADE == 1 ||
              IGRADE == 2 ||
              IGRADE == 3 ||
              IGRADE == 4 ||
              IGRADE == 6) &&
          ISYM == 0) ||
      ((IGRADE == 1 ||
              IGRADE == 2 ||
              IGRADE == 3 ||
              IGRADE == 4 ||
              IGRADE == 5) &&
          ISYM == 2)) {
    INFO.value = -11;
  } else if (IGRADE == 4 && DZERO) {
    INFO.value = -12;
  } else if ((IGRADE == 1 ||
          IGRADE == 3 ||
          IGRADE == 4 ||
          IGRADE == 5 ||
          IGRADE == 6) &&
      (MODEL < -6 || MODEL > 6)) {
    INFO.value = -13;
  } else if ((IGRADE == 1 ||
          IGRADE == 3 ||
          IGRADE == 4 ||
          IGRADE == 5 ||
          IGRADE == 6) &&
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
      ((IPVTNG == 1 || IPVTNG == 2) && (ISYM == 0 || ISYM == 2))) {
    INFO.value = -18;
  } else if (IPVTNG != 0 && BADPVT) {
    INFO.value = -19;
  } else if (KL < 0) {
    INFO.value = -20;
  } else if (KU < 0 || ((ISYM == 0 || ISYM == 2) && KL != KU)) {
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
    xerbla('ZLATMR', -INFO.value);
    return;
  }

  // Decide if we can pivot consistently

  FULBND = false;
  if (KUU == N - 1 && KLL == M - 1) FULBND = true;

  // Initialize random number generator

  for (I = 1; I <= 4; I++) {
    ISEED[I] = ISEED[I].abs() % 4096;
  }

  ISEED[4] = 2 * (ISEED[4] ~/ 2) + 1;

  // 2)      Set up D, DL, and DR, if indicated.

  // Compute D according to COND and MODE

  zlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO);
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
    if (TEMP == ZERO && DMAX != Complex.zero) {
      INFO.value = 2;
      return;
    }
    if (TEMP != ZERO) {
      CALPHA = DMAX / TEMP.toComplex();
    } else {
      CALPHA = Complex.one;
    }
    for (I = 1; I <= MNMIN; I++) {
      D[I] = CALPHA * D[I];
    }
  }

  // If matrix Hermitian, make D real

  if (ISYM == 0) {
    for (I = 1; I <= MNMIN; I++) {
      D[I] = D[I].real.toComplex();
    }
  }

  // Compute DL if grading set

  if (IGRADE == 1 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5 || IGRADE == 6) {
    zlatm1(MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO);
    if (INFO.value != 0) {
      INFO.value = 3;
      return;
    }
  }

  // Compute DR if grading set

  if (IGRADE == 2 || IGRADE == 3) {
    zlatm1(MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO);
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
  //         Always sweep matrix columnwise (if symmetric, upper
  //         half only) so that matrix generated does not depend
  //         on PACK

  if (FULBND) {
    // Use ZLATM3 so matrices generated with differing PIVOTing only
    // differ only in the order of their rows and/or columns.

    if (IPACK == 0) {
      if (ISYM == 0) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            A[ISUB.value][JSUB.value] = CTEMP;
            A[JSUB.value][ISUB.value] = CTEMP.conjugate();
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= M; I++) {
            CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            A[ISUB.value][JSUB.value] = CTEMP;
          }
        }
      } else if (ISYM == 2) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            A[ISUB.value][JSUB.value] = CTEMP;
            A[JSUB.value][ISUB.value] = CTEMP;
          }
        }
      }
    } else if (IPACK == 1) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
              IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          if (MXSUB == ISUB.value && ISYM == 0) {
            A[MNSUB][MXSUB] = CTEMP.conjugate();
          } else {
            A[MNSUB][MXSUB] = CTEMP;
          }
          if (MNSUB != MXSUB) A[MXSUB][MNSUB] = Complex.zero;
        }
      }
    } else if (IPACK == 2) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
              IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          if (MXSUB == JSUB.value && ISYM == 0) {
            A[MXSUB][MNSUB] = CTEMP.conjugate();
          } else {
            A[MXSUB][MNSUB] = CTEMP;
          }
          if (MNSUB != MXSUB) A[MNSUB][MXSUB] = Complex.zero;
        }
      }
    } else if (IPACK == 3) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
              IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);

          // Compute K = location of (ISUB,JSUB) entry in packed
          // array

          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          K = MXSUB * (MXSUB - 1) ~/ 2 + MNSUB;

          // Convert K to (IISUB,JJSUB) location

          JJSUB = (K - 1) ~/ LDA + 1;
          IISUB = K - LDA * (JJSUB - 1);

          if (MXSUB == ISUB.value && ISYM == 0) {
            A[IISUB][JJSUB] = CTEMP.conjugate();
          } else {
            A[IISUB][JJSUB] = CTEMP;
          }
        }
      }
    } else if (IPACK == 4) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
              IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);

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

          if (MXSUB == JSUB.value && ISYM == 0) {
            A[IISUB][JJSUB] = CTEMP.conjugate();
          } else {
            A[IISUB][JJSUB] = CTEMP;
          }
        }
      }
    } else if (IPACK == 5) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          if (I < 1) {
            A[J - I + 1][I + N] = Complex.zero;
          } else {
            CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            MNSUB = min(ISUB.value, JSUB.value);
            MXSUB = max(ISUB.value, JSUB.value);
            if (MXSUB == JSUB.value && ISYM == 0) {
              A[MXSUB - MNSUB + 1][MNSUB] = CTEMP.conjugate();
            } else {
              A[MXSUB - MNSUB + 1][MNSUB] = CTEMP;
            }
          }
        }
      }
    } else if (IPACK == 6) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
              IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          MNSUB = min(ISUB.value, JSUB.value);
          MXSUB = max(ISUB.value, JSUB.value);
          if (MXSUB == ISUB.value && ISYM == 0) {
            A[MNSUB - MXSUB + KUU + 1][MXSUB] = CTEMP.conjugate();
          } else {
            A[MNSUB - MXSUB + KUU + 1][MXSUB] = CTEMP;
          }
        }
      }
    } else if (IPACK == 7) {
      if (ISYM != 1) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J; I++) {
            CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            MNSUB = min(ISUB.value, JSUB.value);
            MXSUB = max(ISUB.value, JSUB.value);
            if (I < 1) A[J - I + 1 + KUU][I + N] = Complex.zero;
            if (MXSUB == ISUB.value && ISYM == 0) {
              A[MNSUB - MXSUB + KUU + 1][MXSUB] = CTEMP.conjugate();
            } else {
              A[MNSUB - MXSUB + KUU + 1][MXSUB] = CTEMP;
            }
            if (I >= 1 && MNSUB != MXSUB) {
              if (MNSUB == ISUB.value && ISYM == 0) {
                A[MXSUB - MNSUB + 1 + KUU][MNSUB] = CTEMP.conjugate();
              } else {
                A[MXSUB - MNSUB + 1 + KUU][MNSUB] = CTEMP;
              }
            }
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J + KLL; I++) {
            CTEMP = zlatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            A[ISUB.value - JSUB.value + KUU + 1][JSUB.value] = CTEMP;
          }
        }
      }
    }
  } else {
    // Use ZLATM2

    if (IPACK == 0) {
      if (ISYM == 0) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            A[I][J] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL,
                DR, IPVTNG, IWORK, SPARSE);
            A[J][I] = A[I][J].conjugate();
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= M; I++) {
            A[I][J] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL,
                DR, IPVTNG, IWORK, SPARSE);
          }
        }
      } else if (ISYM == 2) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            A[I][J] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL,
                DR, IPVTNG, IWORK, SPARSE);
            A[J][I] = A[I][J];
          }
        }
      }
    } else if (IPACK == 1) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          A[I][J] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR,
              IPVTNG, IWORK, SPARSE);
          if (I != J) A[J][I] = Complex.zero;
        }
      }
    } else if (IPACK == 2) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J; I++) {
          if (ISYM == 0) {
            A[J][I] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL,
                    DR, IPVTNG, IWORK, SPARSE)
                .conjugate();
          } else {
            A[J][I] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL,
                DR, IPVTNG, IWORK, SPARSE);
          }
          if (I != J) A[I][J] = Complex.zero;
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
          A[ISUB.value][JSUB.value] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED,
              D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
        }
      }
    } else if (IPACK == 4) {
      if (ISYM == 0 || ISYM == 2) {
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

            A[ISUB.value][JSUB.value] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED,
                D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            if (ISYM == 0) {
              A[ISUB.value][JSUB.value] = A[ISUB.value][JSUB.value].conjugate();
            }
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
            A[ISUB.value][JSUB.value] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED,
                D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          }
        }
      }
    } else if (IPACK == 5) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          if (I < 1) {
            A[J - I + 1][I + N] = Complex.zero;
          } else {
            if (ISYM == 0) {
              A[J - I + 1][I] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
                      IGRADE, DL, DR, IPVTNG, IWORK, SPARSE)
                  .conjugate();
            } else {
              A[J - I + 1][I] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
                  IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            }
          }
        }
      }
    } else if (IPACK == 6) {
      for (J = 1; J <= N; J++) {
        for (I = J - KUU; I <= J; I++) {
          A[I - J + KUU + 1][J] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
              IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
        }
      }
    } else if (IPACK == 7) {
      if (ISYM != 1) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J; I++) {
            A[I - J + KUU + 1][J] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
            if (I < 1) A[J - I + 1 + KUU][I + N] = Complex.zero;
            if (I >= 1 && I != J) {
              if (ISYM == 0) {
                A[J - I + 1 + KUU][I] = A[I - J + KUU + 1][J].conjugate();
              } else {
                A[J - I + 1 + KUU][I] = A[I - J + KUU + 1][J];
              }
            }
          }
        }
      } else if (ISYM == 1) {
        for (J = 1; J <= N; J++) {
          for (I = J - KUU; I <= J + KLL; I++) {
            A[I - J + KUU + 1][J] = zlatm2(M, N, I, J, KL, KU, IDIST, ISEED, D,
                IGRADE, DL, DR, IPVTNG, IWORK, SPARSE);
          }
        }
      }
    }
  }

  // 5)      Scaling the norm

  if (IPACK == 0) {
    ONORM = zlange('M', M, N, A, LDA, TEMPA);
  } else if (IPACK == 1) {
    ONORM = zlansy('M', 'U', N, A, LDA, TEMPA);
  } else if (IPACK == 2) {
    ONORM = zlansy('M', 'L', N, A, LDA, TEMPA);
  } else if (IPACK == 3) {
    ONORM = zlansp('M', 'U', N, A.asArray(), TEMPA);
  } else if (IPACK == 4) {
    ONORM = zlansp('M', 'L', N, A.asArray(), TEMPA);
  } else if (IPACK == 5) {
    ONORM = zlansb('M', 'L', N, KLL, A, LDA, TEMPA);
  } else if (IPACK == 6) {
    ONORM = zlansb('M', 'U', N, KUU, A, LDA, TEMPA);
  } else if (IPACK == 7) {
    ONORM = zlangb('M', N, KLL, KUU, A, LDA, TEMPA);
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
          zdscal(M, ONE / ONORM, A(1, J).asArray(), 1);
          zdscal(M, ANORM, A(1, J).asArray(), 1);
        }
      } else if (IPACK == 3 || IPACK == 4) {
        zdscal(N * (N + 1) ~/ 2, ONE / ONORM, A.asArray(), 1);
        zdscal(N * (N + 1) ~/ 2, ANORM, A.asArray(), 1);
      } else if (IPACK >= 5) {
        for (J = 1; J <= N; J++) {
          zdscal(KLL + KUU + 1, ONE / ONORM, A(1, J).asArray(), 1);
          zdscal(KLL + KUU + 1, ANORM, A(1, J).asArray(), 1);
        }
      }
    } else {
      // Scale straightforwardly

      if (IPACK <= 2) {
        for (J = 1; J <= N; J++) {
          zdscal(M, ANORM / ONORM, A(1, J).asArray(), 1);
        }
      } else if (IPACK == 3 || IPACK == 4) {
        zdscal(N * (N + 1) ~/ 2, ANORM / ONORM, A.asArray(), 1);
      } else if (IPACK >= 5) {
        for (J = 1; J <= N; J++) {
          zdscal(KLL + KUU + 1, ANORM / ONORM, A(1, J).asArray(), 1);
        }
      }
    }
  }
}
