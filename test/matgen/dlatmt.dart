// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlartg.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

import 'dlagge.dart';
import 'dlagsy.dart';
import 'dlarnd.dart';
import 'dlarot.dart';
import 'dlatm7.dart';

void dlatmt(
  final int M,
  final int N,
  final String DIST,
  final Array<int> ISEED_,
  final String SYM,
  final Array<double> D_,
  final int MODE,
  final double COND,
  final double DMAX,
  final int RANK,
  final int KL,
  final int KU,
  final String PACK,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  const ONE = 1.0;
  const TWOPI = 6.28318530717958647692528676655900576839;
  double ALPHA, ANGLE;
  int I,
      IC,
      ICOL = 0,
      IDIST,
      IENDCH,
      IL,
      ILDA,
      IOFFG,
      IOFFST,
      IPACK,
      IPACKG,
      IR,
      IR1,
      IR2,
      IROW = 0,
      IRSIGN = 0,
      ISKEW,
      ISYM,
      ISYMPK,
      J,
      JC,
      JCH,
      JKL,
      JKU,
      JR,
      K,
      LLB,
      MINLDA,
      MNMIN,
      MR,
      NC,
      UUB;
  bool GIVENS, ILEXTR, ILTEMP, TOPDWN;
  final IINFO = Box(0);
  final DUMMY = Box(0.0),
      EXTRA = Box(0.0),
      C = Box(0.0),
      S = Box(0.0),
      TEMP = Box(0.0);

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
  } else {
    IDIST = -1;
  }

  // Decode SYM

  if (lsame(SYM, 'N')) {
    ISYM = 1;
    IRSIGN = 0;
  } else if (lsame(SYM, 'P')) {
    ISYM = 2;
    IRSIGN = 0;
  } else if (lsame(SYM, 'S')) {
    ISYM = 2;
    IRSIGN = 1;
  } else if (lsame(SYM, 'H')) {
    ISYM = 2;
    IRSIGN = 1;
  } else {
    ISYM = -1;
  }

  // Decode PACK

  ISYMPK = 0;
  if (lsame(PACK, 'N')) {
    IPACK = 0;
  } else if (lsame(PACK, 'U')) {
    IPACK = 1;
    ISYMPK = 1;
  } else if (lsame(PACK, 'L')) {
    IPACK = 2;
    ISYMPK = 1;
  } else if (lsame(PACK, 'C')) {
    IPACK = 3;
    ISYMPK = 2;
  } else if (lsame(PACK, 'R')) {
    IPACK = 4;
    ISYMPK = 3;
  } else if (lsame(PACK, 'B')) {
    IPACK = 5;
    ISYMPK = 3;
  } else if (lsame(PACK, 'Q')) {
    IPACK = 6;
    ISYMPK = 2;
  } else if (lsame(PACK, 'Z')) {
    IPACK = 7;
  } else {
    IPACK = -1;
  }

  // Set certain internal parameters

  MNMIN = min(M, N);
  LLB = min(KL, M - 1);
  UUB = min(KU, N - 1);
  MR = min(M, N + LLB);
  NC = min(N, M + UUB);

  if (IPACK == 5 || IPACK == 6) {
    MINLDA = UUB + 1;
  } else if (IPACK == 7) {
    MINLDA = LLB + UUB + 1;
  } else {
    MINLDA = M;
  }

  // Use Givens rotation method if bandwidth small enough,
  // or if LDA is too small to store the matrix unpacked.

  GIVENS = false;
  if (ISYM == 1) {
    if ((LLB + UUB) < 0.3 * max(1, MR + NC)) {
      GIVENS = true;
    }
  } else {
    if (2 * LLB < M) GIVENS = true;
  }
  if (LDA < M && LDA >= MINLDA) GIVENS = true;

  // Set INFO if an error

  if (M < 0) {
    INFO.value = -1;
  } else if (M != N && ISYM != 1) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (IDIST == -1) {
    INFO.value = -3;
  } else if (ISYM == -1) {
    INFO.value = -5;
  } else if (MODE.abs() > 6) {
    INFO.value = -7;
  } else if ((MODE != 0 && MODE.abs() != 6) && COND < ONE) {
    INFO.value = -8;
  } else if (KL < 0) {
    INFO.value = -10;
  } else if (KU < 0 || (ISYM != 1 && KL != KU)) {
    INFO.value = -11;
  } else if (IPACK == -1 ||
      (ISYMPK == 1 && ISYM == 1) ||
      (ISYMPK == 2 && ISYM == 1 && KL > 0) ||
      (ISYMPK == 3 && ISYM == 1 && KU > 0) ||
      (ISYMPK != 0 && M != N)) {
    INFO.value = -12;
  } else if (LDA < max(1, MINLDA)) {
    INFO.value = -14;
  }

  if (INFO.value != 0) {
    xerbla('DLATMT', -INFO.value);
    return;
  }

  // Initialize random number generator

  for (I = 1; I <= 4; I++) {
    ISEED[I] = ISEED[I].abs() % 4096;
  }

  if ((ISEED[4] % 2) != 1) ISEED[4]++;

  // 2)      Set up D  if indicated.

  // Compute D according to COND and MODE

  dlatm7(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, RANK, IINFO);
  if (IINFO.value != 0) {
    INFO.value = 1;
    return;
  }

  // Choose Top-Down if D is (apparently) increasing,
  // Bottom-Up if D is (apparently) decreasing.

  if (D[1].abs() <= D[RANK].abs()) {
    TOPDWN = true;
  } else {
    TOPDWN = false;
  }

  if (MODE != 0 && MODE.abs() != 6) {
    // Scale by DMAX

    TEMP.value = D[1].abs();
    for (I = 2; I <= RANK; I++) {
      TEMP.value = max(TEMP.value, D[I].abs());
    }

    if (TEMP.value > ZERO) {
      ALPHA = DMAX / TEMP.value;
    } else {
      INFO.value = 2;
      return;
    }

    dscal(RANK, ALPHA, D, 1);
  }

  // 3)      Generate Banded Matrix using Givens rotations.
  //         Also the special case of UUB=LLB=0

  // Compute Addressing constants to cover all
  // storage formats.  Whether GE, SY, GB, or SB,
  // upper or lower triangle or both,
  // the (i,j)-th element is in
  // A( i - ISKEW*j + IOFFST, j )

  if (IPACK > 4) {
    ILDA = LDA - 1;
    ISKEW = 1;
    if (IPACK > 5) {
      IOFFST = UUB + 1;
    } else {
      IOFFST = 1;
    }
  } else {
    ILDA = LDA;
    ISKEW = 0;
    IOFFST = 0;
  }

  // IPACKG is the format that the matrix is generated in. If this is
  // different from IPACK, then the matrix must be repacked at the
  // end.  It also signals how to compute the norm, for scaling.

  IPACKG = 0;
  dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);

  // Diagonal Matrix -- We are done, unless it
  // is to be stored SP/PP/TP (PACK='R' or 'C')

  if (LLB == 0 && UUB == 0) {
    dcopy(MNMIN, D, 1, A(1 - ISKEW + IOFFST, 1).asArray(), ILDA + 1);
    if (IPACK <= 2 || IPACK >= 5) IPACKG = IPACK;
  } else if (GIVENS) {
    // Check whether to use Givens rotations,
    // Householder transformations, or nothing.

    if (ISYM == 1) {
      // Non-symmetric -- A = U D V

      if (IPACK > 4) {
        IPACKG = IPACK;
      } else {
        IPACKG = 0;
      }

      dcopy(MNMIN, D, 1, A(1 - ISKEW + IOFFST, 1).asArray(), ILDA + 1);

      if (TOPDWN) {
        JKL = 0;
        for (JKU = 1; JKU <= UUB; JKU++) {
          // Transform from bandwidth JKL, JKU-1 to JKL, JKU

          // Last row actually rotated is M
          // Last column actually rotated is min( M+JKU, N )

          for (JR = 1; JR <= min(M + JKU, N) + JKL - 1; JR++) {
            EXTRA.value = ZERO;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C.value = cos(ANGLE);
            S.value = sin(ANGLE);
            ICOL = max(1, JR - JKL);
            if (JR < M) {
              IL = min(N, JR + JKU) + 1 - ICOL;
              dlarot(
                  true,
                  JR > JKL,
                  false,
                  IL,
                  C.value,
                  S.value,
                  A(JR - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                  ILDA,
                  EXTRA,
                  DUMMY);
            }

            // Chase "EXTRA" back up

            IR = JR;
            IC = ICOL;
            for (JCH = JR - JKL;
                -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1;
                JCH += -JKL - JKU) {
              if (IR < M) {
                dlartg(A[IR + 1 - ISKEW * (IC + 1) + IOFFST][IC + 1],
                    EXTRA.value, C, S, DUMMY);
              }
              IROW = max(1, JCH - JKU);
              IL = IR + 2 - IROW;
              TEMP.value = ZERO;
              ILTEMP = JCH > JKU;
              dlarot(
                  false,
                  ILTEMP,
                  true,
                  IL,
                  C.value,
                  -S.value,
                  A(IROW - ISKEW * IC + IOFFST, IC).asArray(),
                  ILDA,
                  TEMP,
                  EXTRA);
              if (ILTEMP) {
                dlartg(A[IROW + 1 - ISKEW * (IC + 1) + IOFFST][IC + 1],
                    TEMP.value, C, S, DUMMY);
                ICOL = max(1, JCH - JKU - JKL);
                IL = IC + 2 - ICOL;
                EXTRA.value = ZERO;
                dlarot(
                    true,
                    JCH > JKU + JKL,
                    true,
                    IL,
                    C.value,
                    -S.value,
                    A(IROW - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                    ILDA,
                    EXTRA,
                    TEMP);
                IC = ICOL;
                IR = IROW;
              }
            }
          }
        }

        JKU = UUB;
        for (JKL = 1; JKL <= LLB; JKL++) {
          // Transform from bandwidth JKL-1, JKU to JKL, JKU

          for (JC = 1; JC <= min(N + JKL, M) + JKU - 1; JC++) {
            EXTRA.value = ZERO;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C.value = cos(ANGLE);
            S.value = sin(ANGLE);
            IROW = max(1, JC - JKU);
            if (JC < N) {
              IL = min(M, JC + JKL) + 1 - IROW;
              dlarot(
                  false,
                  JC > JKU,
                  false,
                  IL,
                  C.value,
                  S.value,
                  A(IROW - ISKEW * JC + IOFFST, JC).asArray(),
                  ILDA,
                  EXTRA,
                  DUMMY);
            }

            // Chase "EXTRA" back up

            IC = JC;
            IR = IROW;
            for (JCH = JC - JKU;
                -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1;
                JCH += -JKL - JKU) {
              if (IC < N) {
                dlartg(A[IR + 1 - ISKEW * (IC + 1) + IOFFST][IC + 1],
                    EXTRA.value, C, S, DUMMY);
              }
              ICOL = max(1, JCH - JKL);
              IL = IC + 2 - ICOL;
              TEMP.value = ZERO;
              ILTEMP = JCH > JKL;
              dlarot(
                  true,
                  ILTEMP,
                  true,
                  IL,
                  C.value,
                  -S.value,
                  A(IR - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                  ILDA,
                  TEMP,
                  EXTRA);
              if (ILTEMP) {
                dlartg(A[IR + 1 - ISKEW * (ICOL + 1) + IOFFST][ICOL + 1],
                    TEMP.value, C, S, DUMMY);
                IROW = max(1, JCH - JKL - JKU);
                IL = IR + 2 - IROW;
                EXTRA.value = ZERO;
                dlarot(
                    false,
                    JCH > JKL + JKU,
                    true,
                    IL,
                    C.value,
                    -S.value,
                    A(IROW - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                    ILDA,
                    EXTRA,
                    TEMP);
                IC = ICOL;
                IR = IROW;
              }
            }
          }
        }
      } else {
        // Bottom-Up -- Start at the bottom right.

        JKL = 0;
        for (JKU = 1; JKU <= UUB; JKU++) {
          // Transform from bandwidth JKL, JKU-1 to JKL, JKU

          // First row actually rotated is M
          // First column actually rotated is min( M+JKU, N )

          IENDCH = min(M, N + JKL) - 1;
          for (JC = min(M + JKU, N) - 1; JC >= 1 - JKL; JC--) {
            EXTRA.value = ZERO;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C.value = cos(ANGLE);
            S.value = sin(ANGLE);
            IROW = max(1, JC - JKU + 1);
            if (JC > 0) {
              IL = min(M, JC + JKL + 1) + 1 - IROW;
              dlarot(
                  false,
                  false,
                  JC + JKL < M,
                  IL,
                  C.value,
                  S.value,
                  A(IROW - ISKEW * JC + IOFFST, JC).asArray(),
                  ILDA,
                  DUMMY,
                  EXTRA);
            }

            // Chase "EXTRA" back down

            IC = JC;
            for (JCH = JC + JKL;
                JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH;
                JCH += JKL + JKU) {
              ILEXTR = IC > 0;
              if (ILEXTR) {
                dlartg(
                    A[JCH - ISKEW * IC + IOFFST][IC], EXTRA.value, C, S, DUMMY);
              }
              IC = max(1, IC);
              ICOL = min(N - 1, JCH + JKU);
              ILTEMP = JCH + JKU < N;
              TEMP.value = ZERO;
              dlarot(
                  true,
                  ILEXTR,
                  ILTEMP,
                  ICOL + 2 - IC,
                  C.value,
                  S.value,
                  A(JCH - ISKEW * IC + IOFFST, IC).asArray(),
                  ILDA,
                  EXTRA,
                  TEMP);
              if (ILTEMP) {
                dlartg(A[JCH - ISKEW * ICOL + IOFFST][ICOL], TEMP.value, C, S,
                    DUMMY);
                IL = min(IENDCH, JCH + JKL + JKU) + 2 - JCH;
                EXTRA.value = ZERO;
                dlarot(
                    false,
                    true,
                    JCH + JKL + JKU <= IENDCH,
                    IL,
                    C.value,
                    S.value,
                    A(JCH - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                    ILDA,
                    TEMP,
                    EXTRA);
                IC = ICOL;
              }
            }
          }
        }

        JKU = UUB;
        for (JKL = 1; JKL <= LLB; JKL++) {
          // Transform from bandwidth JKL-1, JKU to JKL, JKU

          // First row actually rotated is min( N+JKL, M )
          // First column actually rotated is N

          IENDCH = min(N, M + JKU) - 1;
          for (JR = min(N + JKL, M) - 1; JR >= 1 - JKU; JR--) {
            EXTRA.value = ZERO;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C.value = cos(ANGLE);
            S.value = sin(ANGLE);
            ICOL = max(1, JR - JKL + 1);
            if (JR > 0) {
              IL = min(N, JR + JKU + 1) + 1 - ICOL;
              dlarot(
                  true,
                  false,
                  JR + JKU < N,
                  IL,
                  C.value,
                  S.value,
                  A(JR - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                  ILDA,
                  DUMMY,
                  EXTRA);
            }

            // Chase "EXTRA" back down

            IR = JR;
            for (JCH = JR + JKU;
                JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH;
                JCH += JKL + JKU) {
              ILEXTR = IR > 0;
              if (ILEXTR) {
                dlartg(A[IR - ISKEW * JCH + IOFFST][JCH], EXTRA.value, C, S,
                    DUMMY);
              }
              IR = max(1, IR);
              IROW = min(M - 1, JCH + JKL);
              ILTEMP = JCH + JKL < M;
              TEMP.value = ZERO;
              dlarot(
                  false,
                  ILEXTR,
                  ILTEMP,
                  IROW + 2 - IR,
                  C.value,
                  S.value,
                  A(IR - ISKEW * JCH + IOFFST, JCH).asArray(),
                  ILDA,
                  EXTRA,
                  TEMP);
              if (ILTEMP) {
                dlartg(A[IROW - ISKEW * JCH + IOFFST][JCH], TEMP.value, C, S,
                    DUMMY);
                IL = min(IENDCH, JCH + JKL + JKU) + 2 - JCH;
                EXTRA.value = ZERO;
                dlarot(
                    true,
                    true,
                    JCH + JKL + JKU <= IENDCH,
                    IL,
                    C.value,
                    S.value,
                    A(IROW - ISKEW * JCH + IOFFST, JCH).asArray(),
                    ILDA,
                    TEMP,
                    EXTRA);
                IR = IROW;
              }
            }
          }
        }
      }
    } else {
      // Symmetric -- A = U D U'

      IPACKG = IPACK;
      IOFFG = IOFFST;

      if (TOPDWN) {
        // Top-Down -- Generate Upper triangle only

        if (IPACK >= 5) {
          IPACKG = 6;
          IOFFG = UUB + 1;
        } else {
          IPACKG = 1;
        }
        dcopy(MNMIN, D, 1, A(1 - ISKEW + IOFFG, 1).asArray(), ILDA + 1);

        for (K = 1; K <= UUB; K++) {
          for (JC = 1; JC <= N - 1; JC++) {
            IROW = max(1, JC - K);
            IL = min(JC + 1, K + 2);
            EXTRA.value = ZERO;
            TEMP.value = A[JC - ISKEW * (JC + 1) + IOFFG][JC + 1];
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C.value = cos(ANGLE);
            S.value = sin(ANGLE);
            dlarot(false, JC > K, true, IL, C.value, S.value,
                A(IROW - ISKEW * JC + IOFFG, JC).asArray(), ILDA, EXTRA, TEMP);
            dlarot(true, true, false, min(K, N - JC) + 1, C.value, S.value,
                A((1 - ISKEW) * JC + IOFFG, JC).asArray(), ILDA, TEMP, DUMMY);

            // Chase EXTRA back up the matrix

            ICOL = JC;
            for (JCH = JC - K; JCH >= 1; JCH -= K) {
              dlartg(A[JCH + 1 - ISKEW * (ICOL + 1) + IOFFG][ICOL + 1],
                  EXTRA.value, C, S, DUMMY);
              TEMP.value = A[JCH - ISKEW * (JCH + 1) + IOFFG][JCH + 1];
              dlarot(
                  true,
                  true,
                  true,
                  K + 2,
                  C.value,
                  -S.value,
                  A((1 - ISKEW) * JCH + IOFFG, JCH).asArray(),
                  ILDA,
                  TEMP,
                  EXTRA);
              IROW = max(1, JCH - K);
              IL = min(JCH + 1, K + 2);
              EXTRA.value = ZERO;
              dlarot(
                  false,
                  JCH > K,
                  true,
                  IL,
                  C.value,
                  -S.value,
                  A(IROW - ISKEW * JCH + IOFFG, JCH).asArray(),
                  ILDA,
                  EXTRA,
                  TEMP);
              ICOL = JCH;
            }
          }
        }

        // If we need lower triangle, copy from upper. Note that
        // the order of copying is chosen to work for 'q' -> 'b'

        if (IPACK != IPACKG && IPACK != 3) {
          for (JC = 1; JC <= N; JC++) {
            IROW = IOFFST - ISKEW * JC;
            for (JR = JC; JR <= min(N, JC + UUB); JR++) {
              A[JR + IROW][JC] = A[JC - ISKEW * JR + IOFFG][JR];
            }
          }
          if (IPACK == 5) {
            for (JC = N - UUB + 1; JC <= N; JC++) {
              for (JR = N + 2 - JC; JR <= UUB + 1; JR++) {
                A[JR][JC] = ZERO;
              }
            }
          }
          if (IPACKG == 6) {
            IPACKG = IPACK;
          } else {
            IPACKG = 0;
          }
        }
      } else {
        // Bottom-Up -- Generate Lower triangle only

        if (IPACK >= 5) {
          IPACKG = 5;
          if (IPACK == 6) IOFFG = 1;
        } else {
          IPACKG = 2;
        }
        dcopy(MNMIN, D, 1, A(1 - ISKEW + IOFFG, 1).asArray(), ILDA + 1);

        for (K = 1; K <= UUB; K++) {
          for (JC = N - 1; JC >= 1; JC--) {
            IL = min(N + 1 - JC, K + 2);
            EXTRA.value = ZERO;
            TEMP.value = A[1 + (1 - ISKEW) * JC + IOFFG][JC];
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C.value = cos(ANGLE);
            S.value = -sin(ANGLE);
            dlarot(false, true, N - JC > K, IL, C.value, S.value,
                A((1 - ISKEW) * JC + IOFFG, JC).asArray(), ILDA, TEMP, EXTRA);
            ICOL = max(1, JC - K + 1);
            dlarot(
                true,
                false,
                true,
                JC + 2 - ICOL,
                C.value,
                S.value,
                A(JC - ISKEW * ICOL + IOFFG, ICOL).asArray(),
                ILDA,
                DUMMY,
                TEMP);

            // Chase EXTRA back down the matrix

            ICOL = JC;
            for (JCH = JC + K; JCH <= N - 1; JCH += K) {
              dlartg(A[JCH - ISKEW * ICOL + IOFFG][ICOL], EXTRA.value, C, S,
                  DUMMY);
              TEMP.value = A[1 + (1 - ISKEW) * JCH + IOFFG][JCH];
              dlarot(
                  true,
                  true,
                  true,
                  K + 2,
                  C.value,
                  S.value,
                  A(JCH - ISKEW * ICOL + IOFFG, ICOL).asArray(),
                  ILDA,
                  EXTRA,
                  TEMP);
              IL = min(N + 1 - JCH, K + 2);
              EXTRA.value = ZERO;
              dlarot(
                  false,
                  true,
                  N - JCH > K,
                  IL,
                  C.value,
                  S.value,
                  A((1 - ISKEW) * JCH + IOFFG, JCH).asArray(),
                  ILDA,
                  TEMP,
                  EXTRA);
              ICOL = JCH;
            }
          }
        }

        // If we need upper triangle, copy from lower. Note that
        // the order of copying is chosen to work for 'b' -> 'q'

        if (IPACK != IPACKG && IPACK != 4) {
          for (JC = N; JC >= 1; JC--) {
            IROW = IOFFST - ISKEW * JC;
            for (JR = JC; JR >= max(1, JC - UUB); JR--) {
              A[JR + IROW][JC] = A[JC - ISKEW * JR + IOFFG][JR];
            }
          }
          if (IPACK == 6) {
            for (JC = 1; JC <= UUB; JC++) {
              for (JR = 1; JR <= UUB + 1 - JC; JR++) {
                A[JR][JC] = ZERO;
              }
            }
          }
          if (IPACKG == 5) {
            IPACKG = IPACK;
          } else {
            IPACKG = 0;
          }
        }
      }
    }
  } else {
    // 4)      Generate Banded Matrix by first
    //         Rotating by random Unitary matrices,
    //         then reducing the bandwidth using Householder
    //         transformations.

    // Note: we should get here only if LDA >= N

    if (ISYM == 1) {
      // Non-symmetric -- A = U D V

      dlagge(MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO);
    } else {
      // Symmetric -- A = U D U'

      dlagsy(M, LLB, D, A, LDA, ISEED, WORK, IINFO);
    }
    if (IINFO.value != 0) {
      INFO.value = 3;
      return;
    }
  }

  // 5)      Pack the matrix

  if (IPACK != IPACKG) {
    if (IPACK == 1) {
      // 'U' -- Upper triangular, not packed

      for (J = 1; J <= M; J++) {
        for (I = J + 1; I <= M; I++) {
          A[I][J] = ZERO;
        }
      }
    } else if (IPACK == 2) {
      // 'L' -- Lower triangular, not packed

      for (J = 2; J <= M; J++) {
        for (I = 1; I <= J - 1; I++) {
          A[I][J] = ZERO;
        }
      }
    } else if (IPACK == 3) {
      // 'C' -- Upper triangle packed Columnwise.

      ICOL = 1;
      IROW = 0;
      for (J = 1; J <= M; J++) {
        for (I = 1; I <= J; I++) {
          IROW++;
          if (IROW > LDA) {
            IROW = 1;
            ICOL++;
          }
          A[IROW][ICOL] = A[I][J];
        }
      }
    } else if (IPACK == 4) {
      // 'R' -- Lower triangle packed Columnwise.

      ICOL = 1;
      IROW = 0;
      for (J = 1; J <= M; J++) {
        for (I = J; I <= M; I++) {
          IROW++;
          if (IROW > LDA) {
            IROW = 1;
            ICOL++;
          }
          A[IROW][ICOL] = A[I][J];
        }
      }
    } else if (IPACK >= 5) {
      // 'B' -- The lower triangle is packed as a band matrix.
      // 'Q' -- The upper triangle is packed as a band matrix.
      // 'Z' -- The whole matrix is packed as a band matrix.

      if (IPACK == 5) UUB = 0;
      if (IPACK == 6) LLB = 0;

      for (J = 1; J <= UUB; J++) {
        for (I = min(J + LLB, M); I >= 1; I--) {
          A[I - J + UUB + 1][J] = A[I][J];
        }
      }

      for (J = UUB + 2; J <= N; J++) {
        for (I = J - UUB; I <= min(J + LLB, M); I++) {
          A[I - J + UUB + 1][J] = A[I][J];
        }
      }
    }

    // If packed, zero out extraneous elements.

    // Symmetric/Triangular Packed --
    // zero out everything after A(IROW,ICOL)

    if (IPACK == 3 || IPACK == 4) {
      for (JC = ICOL; JC <= M; JC++) {
        for (JR = IROW + 1; JR <= LDA; JR++) {
          A[JR][JC] = ZERO;
        }
        IROW = 0;
      }
    } else if (IPACK >= 5) {
      // Packed Band --
      //    1st row is now in A( UUB+2-j, j), zero above it
      //    m-th row is now in A( M+UUB-j,j), zero below it
      //    last non-zero diagonal is now in A( UUB+LLB+1,j ),
      //       zero below it, too.

      IR1 = UUB + LLB + 2;
      IR2 = UUB + M + 2;
      for (JC = 1; JC <= N; JC++) {
        for (JR = 1; JR <= UUB + 1 - JC; JR++) {
          A[JR][JC] = ZERO;
        }
        for (JR = max(1, min(IR1, IR2 - JC)); JR <= LDA; JR++) {
          A[JR][JC] = ZERO;
        }
      }
    }
  }
}
