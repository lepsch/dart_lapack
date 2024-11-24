// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/drot.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlag2.dart';
import 'package:dart_lapack/src/dlanhs.dart';
import 'package:dart_lapack/src/dlapy2.dart';
import 'package:dart_lapack/src/dlapy3.dart';
import 'package:dart_lapack/src/dlarfg.dart';
import 'package:dart_lapack/src/dlartg.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dlasv2.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dhgeqz(
  final String JOB,
  final String COMPQ,
  final String COMPZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<double> H_,
  final int LDH,
  final Matrix<double> T_,
  final int LDT,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final T = T_.having(ld: LDT);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const HALF = 0.5, ZERO = 0.0, ONE = 1.0, SAFETY = 1.0e+2;
  bool ILAZR2,
      ILAZRO,
      ILPIVT = false,
      ILQ = false,
      ILSCHR = false,
      ILZ = false,
      LQUERY;
  int ICOMPQ,
      ICOMPZ,
      IFIRST = 0,
      IFRSTM = 0,
      IITER,
      ILAST = 0,
      ILASTM,
      IN,
      ISCHUR,
      ISTART = 0,
      J,
      JC,
      JCH,
      JITER,
      JR,
      MAXIT;
  double A11,
      A12,
      A1I,
      A1R,
      A21,
      A22,
      A2I,
      A2R,
      AD11,
      AD11L,
      AD12,
      AD12L,
      AD21,
      AD21L = 0,
      AD22,
      AD22L,
      AD32L,
      AN,
      ANORM,
      ASCALE,
      ATOL,
      B1A,
      B1I,
      B1R,
      B2A,
      B2I,
      B2R,
      BN,
      BNORM,
      BSCALE,
      BTOL,
      C11I,
      C11R,
      C12,
      C21,
      C22I,
      C22R,
      CQ,
      CZ,
      ESHIFT,
      S1INV,
      SAFMAX,
      SAFMIN,
      SCALE = 0,
      SQI,
      SQR,
      SZI,
      SZR,
      T1,
      T2,
      T3,
      TEMPI,
      U1 = 0,
      U12,
      U12L,
      U2 = 0,
      ULP,
      VS,
      W11,
      W12,
      W21,
      W22,
      WABS;
  final V = Array<double>(3);
  final C = Box(0.0),
      S = Box(0.0),
      TEMPR = Box(0.0),
      TAU = Box(0.0),
      S1 = Box(0.0),
      S2 = Box(0.0),
      WI = Box(0.0),
      WR = Box(0.0),
      WR2 = Box(0.0),
      TEMP = Box(0.0),
      TEMP2 = Box(0.0),
      B11 = Box(0.0),
      B22 = Box(0.0),
      SR = Box(0.0),
      CR = Box(0.0),
      SL = Box(0.0),
      CL = Box(0.0);

  if (lsame(JOB, 'E')) {
    ILSCHR = false;
    ISCHUR = 1;
  } else if (lsame(JOB, 'S')) {
    ILSCHR = true;
    ISCHUR = 2;
  } else {
    ISCHUR = 0;
  }

  if (lsame(COMPQ, 'N')) {
    ILQ = false;
    ICOMPQ = 1;
  } else if (lsame(COMPQ, 'V')) {
    ILQ = true;
    ICOMPQ = 2;
  } else if (lsame(COMPQ, 'I')) {
    ILQ = true;
    ICOMPQ = 3;
  } else {
    ICOMPQ = 0;
  }

  if (lsame(COMPZ, 'N')) {
    ILZ = false;
    ICOMPZ = 1;
  } else if (lsame(COMPZ, 'V')) {
    ILZ = true;
    ICOMPZ = 2;
  } else if (lsame(COMPZ, 'I')) {
    ILZ = true;
    ICOMPZ = 3;
  } else {
    ICOMPZ = 0;
  }

  // Check Argument Values

  INFO.value = 0;
  WORK[1] = max(1, N.toDouble());
  LQUERY = (LWORK == -1);
  if (ISCHUR == 0) {
    INFO.value = -1;
  } else if (ICOMPQ == 0) {
    INFO.value = -2;
  } else if (ICOMPZ == 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (ILO < 1) {
    INFO.value = -5;
  } else if (IHI > N || IHI < ILO - 1) {
    INFO.value = -6;
  } else if (LDH < N) {
    INFO.value = -8;
  } else if (LDT < N) {
    INFO.value = -10;
  } else if (LDQ < 1 || (ILQ && LDQ < N)) {
    INFO.value = -15;
  } else if (LDZ < 1 || (ILZ && LDZ < N)) {
    INFO.value = -17;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -19;
  }
  if (INFO.value != 0) {
    xerbla('DHGEQZ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N <= 0) {
    WORK[1] = ONE;
    return;
  }

  // Initialize Q and Z

  if (ICOMPQ == 3) dlaset('Full', N, N, ZERO, ONE, Q, LDQ);
  if (ICOMPZ == 3) dlaset('Full', N, N, ZERO, ONE, Z, LDZ);

  // Machine Constants

  IN = IHI + 1 - ILO;
  SAFMIN = dlamch('S');
  SAFMAX = ONE / SAFMIN;
  ULP = dlamch('E') * dlamch('B');
  ANORM = dlanhs('F', IN, H(ILO, ILO), LDH, WORK);
  BNORM = dlanhs('F', IN, T(ILO, ILO), LDT, WORK);
  ATOL = max(SAFMIN, ULP * ANORM);
  BTOL = max(SAFMIN, ULP * BNORM);
  ASCALE = ONE / max(SAFMIN, ANORM);
  BSCALE = ONE / max(SAFMIN, BNORM);

  // Set Eigenvalues IHI+1:N

  for (J = IHI + 1; J <= N; J++) {
    if (T[J][J] < ZERO) {
      if (ILSCHR) {
        for (JR = 1; JR <= J; JR++) {
          H[JR][J] = -H[JR][J];
          T[JR][J] = -T[JR][J];
        }
      } else {
        H[J][J] = -H[J][J];
        T[J][J] = -T[J][J];
      }
      if (ILZ) {
        for (JR = 1; JR <= N; JR++) {
          Z[JR][J] = -Z[JR][J];
        }
      }
    }
    ALPHAR[J] = H[J][J];
    ALPHAI[J] = ZERO;
    BETA[J] = T[J][J];
  }

  // If IHI < ILO, skip QZ steps
  if (IHI >= ILO) {
    // MAIN QZ ITERATION LOOP

    // Initialize dynamic indices
    //
    // Eigenvalues ILAST+1:N have been found.
    //    Column operations modify rows IFRSTM:whatever.
    //    Row operations modify columns whatever:ILASTM.
    //
    // If only eigenvalues are being computed, then
    //    IFRSTM is the row of the last splitting row above row ILAST;
    //    this is always at least ILO.
    // IITER counts iterations since the last eigenvalue was found,
    //    to tell when to use an extraordinary shift.
    // MAXIT is the maximum number of QZ sweeps allowed.

    ILAST = IHI;
    if (ILSCHR) {
      IFRSTM = 1;
      ILASTM = N;
    } else {
      IFRSTM = ILO;
      ILASTM = IHI;
    }
    IITER = 0;
    ESHIFT = ZERO;
    MAXIT = 30 * (IHI - ILO + 1);

    var dropThroughNonConvergence = true;
    mainQzLoop:
    for (JITER = 1; JITER <= MAXIT; JITER++) {
      var standardizeB = true;

      // Split the matrix if possible.
      //
      // Two tests:
      //   1: H[j][j-1]=0  or  j=ILO
      //   2: T[j][j]=0

      if (ILAST == ILO) {
        // Special case: j=ILAST
      } else if (H[ILAST][ILAST - 1].abs() <=
          max(SAFMIN,
              ULP * (H[ILAST][ILAST].abs() + H[ILAST - 1][ILAST - 1].abs()))) {
        H[ILAST][ILAST - 1] = ZERO;
      } else {
        var splitOff1x1Block = true;
        if (T[ILAST][ILAST].abs() <= BTOL) {
          T[ILAST][ILAST] = ZERO;
        } else {
          // General case: j<ILAST

          var isDropThroughImpossible = true;
          dropThrough:
          for (J = ILAST - 1; J >= ILO; J--) {
            // Test 1: for H[j][j-1]=0 or j=ILO

            if (J == ILO) {
              ILAZRO = true;
            } else {
              if (H[J][J - 1].abs() <=
                  max(SAFMIN, ULP * (H[J][J].abs() + H[J - 1][J - 1].abs()))) {
                H[J][J - 1] = ZERO;
                ILAZRO = true;
              } else {
                ILAZRO = false;
              }
            }

            // Test 2: for T[j][j]=0

            if (T[J][J].abs() < BTOL) {
              T[J][J] = ZERO;

              // Test 1a: Check for 2 consecutive small subdiagonals in A

              ILAZR2 = false;
              if (!ILAZRO) {
                TEMP.value = H[J][J - 1].abs();
                TEMP2.value = H[J][J].abs();
                TEMPR.value = max(TEMP.value, TEMP2.value);
                if (TEMPR.value < ONE && TEMPR.value != ZERO) {
                  TEMP.value /= TEMPR.value;
                  TEMP2.value /= TEMPR.value;
                }
                if (TEMP.value * (ASCALE * H[J + 1][J].abs()) <=
                    TEMP2.value * (ASCALE * ATOL)) ILAZR2 = true;
              }

              // If both tests pass (1 & 2), i.e., the leading diagonal
              // element of B in the block is zero, split a 1x1 block off
              // at the top. (I.e., at the J-th row/column) The leading
              // diagonal element of the remainder can also be zero, so
              // this may have to be done repeatedly.

              if (ILAZRO || ILAZR2) {
                for (JCH = J; JCH <= ILAST - 1; JCH++) {
                  TEMP.value = H[JCH][JCH];
                  dlartg(TEMP.value, H[JCH + 1][JCH], C, S, H.box(JCH, JCH));
                  H[JCH + 1][JCH] = ZERO;
                  drot(ILASTM - JCH, H(JCH, JCH + 1).asArray(), LDH,
                      H(JCH + 1, JCH + 1).asArray(), LDH, C.value, S.value);
                  drot(ILASTM - JCH, T(JCH, JCH + 1).asArray(), LDT,
                      T(JCH + 1, JCH + 1).asArray(), LDT, C.value, S.value);
                  if (ILQ) {
                    drot(N, Q(1, JCH).asArray(), 1, Q(1, JCH + 1).asArray(), 1,
                        C.value, S.value);
                  }
                  if (ILAZR2) H[JCH][JCH - 1] *= C.value;
                  ILAZR2 = false;
                  if (T[JCH + 1][JCH + 1].abs() >= BTOL) {
                    isDropThroughImpossible = false;
                    splitOff1x1Block = false;
                    if (JCH + 1 < ILAST) {
                      IFIRST = JCH + 1;
                      standardizeB = false;
                    }
                    break dropThrough;
                  }
                  T[JCH + 1][JCH + 1] = ZERO;
                }
                isDropThroughImpossible = false;
                break;
              }

              // Only test 2 passed -- chase the zero to T[ILAST][ILAST]
              // Then process as in the case T[ILAST][ILAST]=0

              for (JCH = J; JCH <= ILAST - 1; JCH++) {
                TEMP.value = T[JCH][JCH + 1];
                dlartg(
                    TEMP.value, T[JCH + 1][JCH + 1], C, S, T.box(JCH, JCH + 1));
                T[JCH + 1][JCH + 1] = ZERO;
                if (JCH < ILASTM - 1) {
                  drot(ILASTM - JCH - 1, T(JCH, JCH + 2).asArray(), LDT,
                      T(JCH + 1, JCH + 2).asArray(), LDT, C.value, S.value);
                }
                drot(ILASTM - JCH + 2, H(JCH, JCH - 1).asArray(), LDH,
                    H(JCH + 1, JCH - 1).asArray(), LDH, C.value, S.value);
                if (ILQ) {
                  drot(N, Q(1, JCH).asArray(), 1, Q(1, JCH + 1).asArray(), 1,
                      C.value, S.value);
                }
                TEMP.value = H[JCH + 1][JCH];
                dlartg(
                    TEMP.value, H[JCH + 1][JCH - 1], C, S, H.box(JCH + 1, JCH));
                H[JCH + 1][JCH - 1] = ZERO;
                drot(JCH + 1 - IFRSTM, H(IFRSTM, JCH).asArray(), 1,
                    H(IFRSTM, JCH - 1).asArray(), 1, C.value, S.value);
                drot(JCH - IFRSTM, T(IFRSTM, JCH).asArray(), 1,
                    T(IFRSTM, JCH - 1).asArray(), 1, C.value, S.value);
                if (ILZ) {
                  drot(N, Z(1, JCH).asArray(), 1, Z(1, JCH - 1).asArray(), 1,
                      C.value, S.value);
                }
              }
              isDropThroughImpossible = false;
              break;
            } else if (ILAZRO) {
              // Only test 1 passed -- work on J:ILAST

              IFIRST = J;
              isDropThroughImpossible = false;
              splitOff1x1Block = false;
              standardizeB = false;
              break;
            }

            // Neither test passed -- try next J
          }

          if (isDropThroughImpossible) {
            // (Drop-through is "impossible")

            INFO.value = N + 1;
            WORK[1] = N.toDouble();
            return;
          }
        }

        if (splitOff1x1Block) {
          // T[ILAST][ILAST]=0 -- clear H[ILAST][ILAST-1] to split off a
          // 1x1 block.

          TEMP.value = H[ILAST][ILAST];
          dlartg(TEMP.value, H[ILAST][ILAST - 1], C, S, H.box(ILAST, ILAST));
          H[ILAST][ILAST - 1] = ZERO;
          drot(ILAST - IFRSTM, H(IFRSTM, ILAST).asArray(), 1,
              H(IFRSTM, ILAST - 1).asArray(), 1, C.value, S.value);
          drot(ILAST - IFRSTM, T(IFRSTM, ILAST).asArray(), 1,
              T(IFRSTM, ILAST - 1).asArray(), 1, C.value, S.value);
          if (ILZ) {
            drot(N, Z(1, ILAST).asArray(), 1, Z(1, ILAST - 1).asArray(), 1,
                C.value, S.value);
          }
        }
      }

      if (standardizeB) {
        // H[ILAST][ILAST-1]=0 -- Standardize B, set ALPHAR, ALPHAI, and BETA

        if (T[ILAST][ILAST] < ZERO) {
          if (ILSCHR) {
            for (J = IFRSTM; J <= ILAST; J++) {
              H[J][ILAST] = -H[J][ILAST];
              T[J][ILAST] = -T[J][ILAST];
            }
          } else {
            H[ILAST][ILAST] = -H[ILAST][ILAST];
            T[ILAST][ILAST] = -T[ILAST][ILAST];
          }
          if (ILZ) {
            for (J = 1; J <= N; J++) {
              Z[J][ILAST] = -Z[J][ILAST];
            }
          }
        }
        ALPHAR[ILAST] = H[ILAST][ILAST];
        ALPHAI[ILAST] = ZERO;
        BETA[ILAST] = T[ILAST][ILAST];

        // Gotonext block -- exit if finished.

        ILAST--;
        if (ILAST < ILO) {
          dropThroughNonConvergence = false;
          break mainQzLoop;
        }

        // Reset counters

        IITER = 0;
        ESHIFT = ZERO;
        if (!ILSCHR) {
          ILASTM = ILAST;
          if (IFRSTM > ILAST) IFRSTM = ILO;
        }
        continue mainQzLoop;
      }

      // QZ step

      // This iteration only involves rows/columns IFIRST:ILAST. We
      // assume IFIRST < ILAST, and that the diagonal of B is non-zero.

      IITER++;
      if (!ILSCHR) {
        IFRSTM = IFIRST;
      }

      // Compute single shifts.

      // At this point, IFIRST < ILAST, and the diagonal elements of
      // T[IFIRST:ILAST,IFIRST][ILAST] are larger than BTOL (in
      // magnitude)

      var useFrancisDoubleShift = false;
      if ((IITER ~/ 10) * 10 == IITER) {
        // Exceptional shift.  Chosen for no particularly good reason.
        // (Single shift only.)

        if ((MAXIT * SAFMIN) * H[ILAST][ILAST - 1].abs() <
            T[ILAST - 1][ILAST - 1].abs()) {
          ESHIFT = H[ILAST][ILAST - 1] / T[ILAST - 1][ILAST - 1];
        } else {
          ESHIFT += ONE / (SAFMIN * MAXIT);
        }
        S1.value = ONE;
        WR.value = ESHIFT;
      } else {
        // Shifts based on the generalized eigenvalues of the
        // bottom-right 2x2 block of A and B. The first eigenvalue
        // returned by DLAG2 is the Wilkinson shift (AEP p.512),

        dlag2(H(ILAST - 1, ILAST - 1), LDH, T(ILAST - 1, ILAST - 1), LDT,
            SAFMIN * SAFETY, S1, S2, WR, WR2, WI);

        if (((WR.value / S1.value) * T[ILAST][ILAST] - H[ILAST][ILAST]).abs() >
            ((WR2.value / S2.value) * T[ILAST][ILAST] - H[ILAST][ILAST])
                .abs()) {
          TEMP.value = WR.value;
          WR.value = WR2.value;
          WR2.value = TEMP.value;
          TEMP.value = S1.value;
          S1.value = S2.value;
          S2.value = TEMP.value;
        }
        TEMP.value = max(
            S1.value, SAFMIN * max(ONE, max(WR.value.abs(), WI.value.abs())));
        if (WI.value != ZERO) {
          useFrancisDoubleShift = true;
        }
      }

      if (!useFrancisDoubleShift) {
        // Fiddle with shift to avoid overflow

        TEMP.value = min(ASCALE, ONE) * (HALF * SAFMAX);
        if (S1.value > TEMP.value) {
          SCALE = TEMP.value / S1.value;
        } else {
          SCALE = ONE;
        }

        TEMP.value = min(BSCALE, ONE) * (HALF * SAFMAX);
        if (WR.value.abs() > TEMP.value) {
          SCALE = min(SCALE, TEMP.value / WR.value.abs());
        }
        S1.value = SCALE * S1.value;
        WR.value = SCALE * WR.value;

        // Now check for two consecutive small subdiagonals.
        var hasConsecutiveSmallSubdiagonals = false;
        for (J = ILAST - 1; J >= IFIRST + 1; J--) {
          ISTART = J;
          TEMP.value = (S1.value * H[J][J - 1]).abs();
          TEMP2.value = (S1.value * H[J][J] - WR.value * T[J][J]).abs();
          TEMPR.value = max(TEMP.value, TEMP2.value);
          if (TEMPR.value < ONE && TEMPR.value != ZERO) {
            TEMP.value /= TEMPR.value;
            TEMP2.value /= TEMPR.value;
          }
          if (((ASCALE * H[J + 1][J]) * TEMP.value).abs() <=
              (ASCALE * ATOL) * TEMP2.value) {
            hasConsecutiveSmallSubdiagonals = true;
            break;
          }
        }

        if (!hasConsecutiveSmallSubdiagonals) {
          ISTART = IFIRST;
        }

        // Do an implicit single-shift QZ sweep.

        // Initial Q

        TEMP.value =
            S1.value * H[ISTART][ISTART] - WR.value * T[ISTART][ISTART];
        TEMP2.value = S1.value * H[ISTART + 1][ISTART];
        dlartg(TEMP.value, TEMP2.value, C, S, TEMPR);

        // Sweep

        for (J = ISTART; J <= ILAST - 1; J++) {
          if (J > ISTART) {
            TEMP.value = H[J][J - 1];
            dlartg(TEMP.value, H[J + 1][J - 1], C, S, H.box(J, J - 1));
            H[J + 1][J - 1] = ZERO;
          }

          for (JC = J; JC <= ILASTM; JC++) {
            TEMP.value = C.value * H[J][JC] + S.value * H[J + 1][JC];
            H[J + 1][JC] = -S.value * H[J][JC] + C.value * H[J + 1][JC];
            H[J][JC] = TEMP.value;
            TEMP2.value = C.value * T[J][JC] + S.value * T[J + 1][JC];
            T[J + 1][JC] = -S.value * T[J][JC] + C.value * T[J + 1][JC];
            T[J][JC] = TEMP2.value;
          }
          if (ILQ) {
            for (JR = 1; JR <= N; JR++) {
              TEMP.value = C.value * Q[JR][J] + S.value * Q[JR][J + 1];
              Q[JR][J + 1] = -S.value * Q[JR][J] + C.value * Q[JR][J + 1];
              Q[JR][J] = TEMP.value;
            }
          }

          TEMP.value = T[J + 1][J + 1];
          dlartg(TEMP.value, T[J + 1][J], C, S, T.box(J + 1, J + 1));
          T[J + 1][J] = ZERO;

          for (JR = IFRSTM; JR <= min(J + 2, ILAST); JR++) {
            TEMP.value = C.value * H[JR][J + 1] + S.value * H[JR][J];
            H[JR][J] = -S.value * H[JR][J + 1] + C.value * H[JR][J];
            H[JR][J + 1] = TEMP.value;
          }
          for (JR = IFRSTM; JR <= J; JR++) {
            TEMP.value = C.value * T[JR][J + 1] + S.value * T[JR][J];
            T[JR][J] = -S.value * T[JR][J + 1] + C.value * T[JR][J];
            T[JR][J + 1] = TEMP.value;
          }
          if (ILZ) {
            for (JR = 1; JR <= N; JR++) {
              TEMP.value = C.value * Z[JR][J + 1] + S.value * Z[JR][J];
              Z[JR][J] = -S.value * Z[JR][J + 1] + C.value * Z[JR][J];
              Z[JR][J + 1] = TEMP.value;
            }
          }
        }

        continue mainQzLoop;
      }

      // Use Francis double-shift

      // Note: the Francis double-shift should work with real shifts,
      // but only if the block is at least 3x3.
      // This code may break if this point is reached with
      // a 2x2 block with real eigenvalues.

      if (IFIRST + 1 == ILAST) {
        // Special case -- 2x2 block with complex eigenvectors

        // Step 1: Standardize, that is, rotate so that

        //     ( B11  0  )
        // B = (         )  with B11 non-negative.
        //     (  0  B22 )

        dlasv2(T[ILAST - 1][ILAST - 1], T[ILAST - 1][ILAST], T[ILAST][ILAST],
            B22, B11, SR, CR, SL, CL);

        if (B11.value < ZERO) {
          CR.value = -CR.value;
          SR.value = -SR.value;
          B11.value = -B11.value;
          B22.value = -B22.value;
        }

        drot(ILASTM + 1 - IFIRST, H(ILAST - 1, ILAST - 1).asArray(), LDH,
            H(ILAST, ILAST - 1).asArray(), LDH, CL.value, SL.value);
        drot(ILAST + 1 - IFRSTM, H(IFRSTM, ILAST - 1).asArray(), 1,
            H(IFRSTM, ILAST).asArray(), 1, CR.value, SR.value);

        if (ILAST < ILASTM) {
          drot(ILASTM - ILAST, T(ILAST - 1, ILAST + 1).asArray(), LDT,
              T(ILAST, ILAST + 1).asArray(), LDT, CL.value, SL.value);
        }
        if (IFRSTM < ILAST - 1) {
          drot(IFIRST - IFRSTM, T(IFRSTM, ILAST - 1).asArray(), 1,
              T(IFRSTM, ILAST).asArray(), 1, CR.value, SR.value);
        }

        if (ILQ) {
          drot(N, Q(1, ILAST - 1).asArray(), 1, Q(1, ILAST).asArray(), 1,
              CL.value, SL.value);
        }
        if (ILZ) {
          drot(N, Z(1, ILAST - 1).asArray(), 1, Z(1, ILAST).asArray(), 1,
              CR.value, SR.value);
        }

        T[ILAST - 1][ILAST - 1] = B11.value;
        T[ILAST - 1][ILAST] = ZERO;
        T[ILAST][ILAST - 1] = ZERO;
        T[ILAST][ILAST] = B22.value;

        // If B22 is negative, negate column ILAST

        if (B22.value < ZERO) {
          for (J = IFRSTM; J <= ILAST; J++) {
            H[J][ILAST] = -H[J][ILAST];
            T[J][ILAST] = -T[J][ILAST];
          }

          if (ILZ) {
            for (J = 1; J <= N; J++) {
              Z[J][ILAST] = -Z[J][ILAST];
            }
          }
          B22.value = -B22.value;
        }

        // Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)

        // Recompute shift

        dlag2(H(ILAST - 1, ILAST - 1), LDH, T(ILAST - 1, ILAST - 1), LDT,
            SAFMIN * SAFETY, S1, TEMP, WR, TEMP2, WI);

        // If standardization has perturbed the shift onto real line,
        // do another (real single-shift) QR step.

        if (WI.value == ZERO) continue mainQzLoop;
        S1INV = ONE / S1.value;

        // Do EISPACK (QZVAL) computation of alpha and beta

        A11 = H[ILAST - 1][ILAST - 1];
        A21 = H[ILAST][ILAST - 1];
        A12 = H[ILAST - 1][ILAST];
        A22 = H[ILAST][ILAST];

        // Compute complex Givens rotation on right
        // (Assume some element of C = (sA - wB) > unfl )
        // __
        // (sA - wB) ( CZ   -SZ )
        // ( SZ    CZ )

        C11R = S1.value * A11 - WR.value * B11.value;
        C11I = -WI.value * B11.value;
        C12 = S1.value * A12;
        C21 = S1.value * A21;
        C22R = S1.value * A22 - WR.value * B22.value;
        C22I = -WI.value * B22.value;

        if (C11R.abs() + C11I.abs() + C12.abs() >
            C21.abs() + C22R.abs() + C22I.abs()) {
          T1 = dlapy3(C12, C11R, C11I);
          CZ = C12 / T1;
          SZR = -C11R / T1;
          SZI = -C11I / T1;
        } else {
          CZ = dlapy2(C22R, C22I);
          if (CZ <= SAFMIN) {
            CZ = ZERO;
            SZR = ONE;
            SZI = ZERO;
          } else {
            TEMPR.value = C22R / CZ;
            TEMPI = C22I / CZ;
            T1 = dlapy2(CZ, C21);
            CZ /= T1;
            SZR = -C21 * TEMPR.value / T1;
            SZI = C21 * TEMPI / T1;
          }
        }

        // Compute Givens rotation on left

        // (  CQ   SQ )
        // (  __      )  A or B
        // ( -SQ   CQ )

        AN = A11.abs() + A12.abs() + A21.abs() + A22.abs();
        BN = B11.value.abs() + B22.value.abs();
        WABS = WR.value.abs() + WI.value.abs();
        if (S1.value * AN > WABS * BN) {
          CQ = CZ * B11.value;
          SQR = SZR * B22.value;
          SQI = -SZI * B22.value;
        } else {
          A1R = CZ * A11 + SZR * A12;
          A1I = SZI * A12;
          A2R = CZ * A21 + SZR * A22;
          A2I = SZI * A22;
          CQ = dlapy2(A1R, A1I);
          if (CQ <= SAFMIN) {
            CQ = ZERO;
            SQR = ONE;
            SQI = ZERO;
          } else {
            TEMPR.value = A1R / CQ;
            TEMPI = A1I / CQ;
            SQR = TEMPR.value * A2R + TEMPI * A2I;
            SQI = TEMPI * A2R - TEMPR.value * A2I;
          }
        }
        T1 = dlapy3(CQ, SQR, SQI);
        CQ /= T1;
        SQR /= T1;
        SQI /= T1;

        // Compute diagonal elements of QBZ

        TEMPR.value = SQR * SZR - SQI * SZI;
        TEMPI = SQR * SZI + SQI * SZR;
        B1R = CQ * CZ * B11.value + TEMPR.value * B22.value;
        B1I = TEMPI * B22.value;
        B1A = dlapy2(B1R, B1I);
        B2R = CQ * CZ * B22.value + TEMPR.value * B11.value;
        B2I = -TEMPI * B11.value;
        B2A = dlapy2(B2R, B2I);

        // Normalize so beta > 0, and Im( alpha1 ) > 0

        BETA[ILAST - 1] = B1A;
        BETA[ILAST] = B2A;
        ALPHAR[ILAST - 1] = (WR.value * B1A) * S1INV;
        ALPHAI[ILAST - 1] = (WI.value * B1A) * S1INV;
        ALPHAR[ILAST] = (WR.value * B2A) * S1INV;
        ALPHAI[ILAST] = -(WI.value * B2A) * S1INV;

        // Step 3: Go to next block -- exit if finished.

        ILAST = IFIRST - 1;
        if (ILAST < ILO) {
          dropThroughNonConvergence = false;
          break mainQzLoop;
        }

        // Reset counters

        IITER = 0;
        ESHIFT = ZERO;
        if (!ILSCHR) {
          ILASTM = ILAST;
          if (IFRSTM > ILAST) IFRSTM = ILO;
        }
        continue mainQzLoop;
      } else {
        // Usual case: 3x3 or larger block, using Francis implicit
        //             double-shift
        //
        //                          2
        // Eigenvalue equation is  w  - c w + d = 0,
        //
        //                               -1 2        -1
        // so compute 1st column of  (A B  )  - c A B   + d
        // using the formula in QZIT (from EISPACK)
        //
        // We assume that the block is at least 3x3

        AD11 = (ASCALE * H[ILAST - 1][ILAST - 1]) /
            (BSCALE * T[ILAST - 1][ILAST - 1]);
        AD21 =
            (ASCALE * H[ILAST][ILAST - 1]) / (BSCALE * T[ILAST - 1][ILAST - 1]);
        AD12 = (ASCALE * H[ILAST - 1][ILAST]) / (BSCALE * T[ILAST][ILAST]);
        AD22 = (ASCALE * H[ILAST][ILAST]) / (BSCALE * T[ILAST][ILAST]);
        U12 = T[ILAST - 1][ILAST] / T[ILAST][ILAST];
        AD11L = (ASCALE * H[IFIRST][IFIRST]) / (BSCALE * T[IFIRST][IFIRST]);
        AD21L = (ASCALE * H[IFIRST + 1][IFIRST]) / (BSCALE * T[IFIRST][IFIRST]);
        AD12L = (ASCALE * H[IFIRST][IFIRST + 1]) /
            (BSCALE * T[IFIRST + 1][IFIRST + 1]);
        AD22L = (ASCALE * H[IFIRST + 1][IFIRST + 1]) /
            (BSCALE * T[IFIRST + 1][IFIRST + 1]);
        AD32L = (ASCALE * H[IFIRST + 2][IFIRST + 1]) /
            (BSCALE * T[IFIRST + 1][IFIRST + 1]);
        U12L = T[IFIRST][IFIRST + 1] / T[IFIRST + 1][IFIRST + 1];

        V[1] = (AD11 - AD11L) * (AD22 - AD11L) -
            AD12 * AD21 +
            AD21 * U12 * AD11L +
            (AD12L - AD11L * U12L) * AD21L;
        V[2] = ((AD22L - AD11L) -
                AD21L * U12L -
                (AD11 - AD11L) -
                (AD22 - AD11L) +
                AD21 * U12) *
            AD21L;
        V[3] = AD32L * AD21L;

        ISTART = IFIRST;

        dlarfg(3, V.box(1), V(2), 1, TAU);
        V[1] = ONE;

        // Sweep

        for (J = ISTART; J <= ILAST - 2; J++) {
          // All but last elements: use 3x3 Householder transforms.

          // Zero (j-1)st column of A

          if (J > ISTART) {
            V[1] = H[J][J - 1];
            V[2] = H[J + 1][J - 1];
            V[3] = H[J + 2][J - 1];

            dlarfg(3, H.box(J, J - 1), V(2), 1, TAU);
            V[1] = ONE;
            H[J + 1][J - 1] = ZERO;
            H[J + 2][J - 1] = ZERO;
          }

          T2 = TAU.value * V[2];
          T3 = TAU.value * V[3];
          for (JC = J; JC <= ILASTM; JC++) {
            TEMP.value = H[J][JC] + V[2] * H[J + 1][JC] + V[3] * H[J + 2][JC];
            H[J][JC] -= TEMP.value * TAU.value;
            H[J + 1][JC] -= TEMP.value * T2;
            H[J + 2][JC] -= TEMP.value * T3;
            TEMP2.value = T[J][JC] + V[2] * T[J + 1][JC] + V[3] * T[J + 2][JC];
            T[J][JC] -= TEMP2.value * TAU.value;
            T[J + 1][JC] -= TEMP2.value * T2;
            T[J + 2][JC] -= TEMP2.value * T3;
          }
          if (ILQ) {
            for (JR = 1; JR <= N; JR++) {
              TEMP.value = Q[JR][J] + V[2] * Q[JR][J + 1] + V[3] * Q[JR][J + 2];
              Q[JR][J] -= TEMP.value * TAU.value;
              Q[JR][J + 1] -= TEMP.value * T2;
              Q[JR][J + 2] -= TEMP.value * T3;
            }
          }

          // Zero j-th column of B (see DLAGBC for details)

          // Swap rows to pivot

          ILPIVT = false;
          TEMP.value = max(T[J + 1][J + 1].abs(), T[J + 1][J + 2].abs());
          TEMP2.value = max(T[J + 2][J + 1].abs(), T[J + 2][J + 2].abs());
          if (max(TEMP.value, TEMP2.value) < SAFMIN) {
            SCALE = ZERO;
            U1 = ONE;
            U2 = ZERO;
            continue;
          } else if (TEMP.value >= TEMP2.value) {
            W11 = T[J + 1][J + 1];
            W21 = T[J + 2][J + 1];
            W12 = T[J + 1][J + 2];
            W22 = T[J + 2][J + 2];
            U1 = T[J + 1][J];
            U2 = T[J + 2][J];
          } else {
            W21 = T[J + 1][J + 1];
            W11 = T[J + 2][J + 1];
            W22 = T[J + 1][J + 2];
            W12 = T[J + 2][J + 2];
            U2 = T[J + 1][J];
            U1 = T[J + 2][J];
          }

          // Swap columns if nec.

          if (W12.abs() > W11.abs()) {
            ILPIVT = true;
            TEMP.value = W12;
            TEMP2.value = W22;
            W12 = W11;
            W22 = W21;
            W11 = TEMP.value;
            W21 = TEMP2.value;
          }

          // LU-factor

          TEMP.value = W21 / W11;
          U2 -= TEMP.value * U1;
          W22 -= TEMP.value * W12;
          W21 = ZERO;

          // Compute SCALE

          SCALE = ONE;
          if (W22.abs() < SAFMIN) {
            SCALE = ZERO;
            U2 = ONE;
            U1 = -W12 / W11;
          } else {
            if (W22.abs() < U2.abs()) SCALE = (W22 / U2).abs();
            if (W11.abs() < U1.abs()) SCALE = min(SCALE, (W11 / U1).abs());

            // Solve

            U2 = (SCALE * U2) / W22;
            U1 = (SCALE * U1 - W12 * U2) / W11;
          }

          if (ILPIVT) {
            TEMP.value = U2;
            U2 = U1;
            U1 = TEMP.value;
          }

          // Compute Householder Vector

          T1 = sqrt(pow(SCALE, 2) + pow(U1, 2) + pow(U2, 2));
          TAU.value = ONE + SCALE / T1;
          VS = -ONE / (SCALE + T1);
          V[1] = ONE;
          V[2] = VS * U1;
          V[3] = VS * U2;

          // Apply transformations from the right.

          T2 = TAU.value * V[2];
          T3 = TAU.value * V[3];
          for (JR = IFRSTM; JR <= min(J + 3, ILAST); JR++) {
            TEMP.value = H[JR][J] + V[2] * H[JR][J + 1] + V[3] * H[JR][J + 2];
            H[JR][J] -= TEMP.value * TAU.value;
            H[JR][J + 1] -= TEMP.value * T2;
            H[JR][J + 2] -= TEMP.value * T3;
          }
          for (JR = IFRSTM; JR <= J + 2; JR++) {
            TEMP.value = T[JR][J] + V[2] * T[JR][J + 1] + V[3] * T[JR][J + 2];
            T[JR][J] -= TEMP.value * TAU.value;
            T[JR][J + 1] -= TEMP.value * T2;
            T[JR][J + 2] -= TEMP.value * T3;
          }
          if (ILZ) {
            for (JR = 1; JR <= N; JR++) {
              TEMP.value = Z[JR][J] + V[2] * Z[JR][J + 1] + V[3] * Z[JR][J + 2];
              Z[JR][J] -= TEMP.value * TAU.value;
              Z[JR][J + 1] -= TEMP.value * T2;
              Z[JR][J + 2] -= TEMP.value * T3;
            }
          }
          T[J + 1][J] = ZERO;
          T[J + 2][J] = ZERO;
        }

        // Last elements: Use Givens rotations

        // Rotations from the left

        J = ILAST - 1;
        TEMP.value = H[J][J - 1];
        dlartg(TEMP.value, H[J + 1][J - 1], C, S, H.box(J, J - 1));
        H[J + 1][J - 1] = ZERO;

        for (JC = J; JC <= ILASTM; JC++) {
          TEMP.value = C.value * H[J][JC] + S.value * H[J + 1][JC];
          H[J + 1][JC] = -S.value * H[J][JC] + C.value * H[J + 1][JC];
          H[J][JC] = TEMP.value;
          TEMP2.value = C.value * T[J][JC] + S.value * T[J + 1][JC];
          T[J + 1][JC] = -S.value * T[J][JC] + C.value * T[J + 1][JC];
          T[J][JC] = TEMP2.value;
        }
        if (ILQ) {
          for (JR = 1; JR <= N; JR++) {
            TEMP.value = C.value * Q[JR][J] + S.value * Q[JR][J + 1];
            Q[JR][J + 1] = -S.value * Q[JR][J] + C.value * Q[JR][J + 1];
            Q[JR][J] = TEMP.value;
          }
        }

        // Rotations from the right.

        TEMP.value = T[J + 1][J + 1];
        dlartg(TEMP.value, T[J + 1][J], C, S, T.box(J + 1, J + 1));
        T[J + 1][J] = ZERO;

        for (JR = IFRSTM; JR <= ILAST; JR++) {
          TEMP.value = C.value * H[JR][J + 1] + S.value * H[JR][J];
          H[JR][J] = -S.value * H[JR][J + 1] + C.value * H[JR][J];
          H[JR][J + 1] = TEMP.value;
        }
        for (JR = IFRSTM; JR <= ILAST - 1; JR++) {
          TEMP.value = C.value * T[JR][J + 1] + S.value * T[JR][J];
          T[JR][J] = -S.value * T[JR][J + 1] + C.value * T[JR][J];
          T[JR][J + 1] = TEMP.value;
        }
        if (ILZ) {
          for (JR = 1; JR <= N; JR++) {
            TEMP.value = C.value * Z[JR][J + 1] + S.value * Z[JR][J];
            Z[JR][J] = -S.value * Z[JR][J + 1] + C.value * Z[JR][J];
            Z[JR][J + 1] = TEMP.value;
          }
        }

        // End of Double-Shift code
      }

      // End of iteration loop
    }

    if (dropThroughNonConvergence) {
      // Drop-through = non-convergence

      INFO.value = ILAST;
      WORK[1] = N.toDouble();
      return;
    }

    // Successful completion of all QZ steps
  }

  // Set Eigenvalues 1:ILO-1

  for (J = 1; J <= ILO - 1; J++) {
    if (T[J][J] < ZERO) {
      if (ILSCHR) {
        for (JR = 1; JR <= J; JR++) {
          H[JR][J] = -H[JR][J];
          T[JR][J] = -T[JR][J];
        }
      } else {
        H[J][J] = -H[J][J];
        T[J][J] = -T[J][J];
      }
      if (ILZ) {
        for (JR = 1; JR <= N; JR++) {
          Z[JR][J] = -Z[JR][J];
        }
      }
    }
    ALPHAR[J] = H[J][J];
    ALPHAI[J] = ZERO;
    BETA[J] = T[J][J];
  }

  // Normal Termination

  INFO.value = 0;

  // Exit (other than argument error) -- return optimal workspace size

  WORK[1] = N.toDouble();
}
