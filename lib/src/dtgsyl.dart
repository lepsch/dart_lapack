// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dtgsy2.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtgsyl(
  final String TRANS,
  final int IJOB,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> C_,
  final int LDC,
  final Matrix<double> D_,
  final int LDD,
  final Matrix<double> E_,
  final int LDE,
  final Matrix<double> F_,
  final int LDF,
  final Box<double> SCALE,
  final Box<double> DIF,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final D = D_.having(ld: LDD);
  final E = E_.having(ld: LDE);
  final F = F_.having(ld: LDF);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, NOTRAN;
  int I, IE, IFUNC, IROUND, IS, ISOLVE, J, JE, JS, K, LWMIN = 0, MB, NB, P, Q;
  double SCALE2 = 0;
  final PQ = Box(0), PPQQ = Box(0), LINFO = Box(0);
  final SCALOC = Box(0.0), DSCALE = Box(0.0), DSUM = Box(0.0);

  // Decode and test input parameters

  INFO.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  LQUERY = (LWORK == -1);

  if (!NOTRAN && !lsame(TRANS, 'T')) {
    INFO.value = -1;
  } else if (NOTRAN) {
    if ((IJOB < 0) || (IJOB > 4)) {
      INFO.value = -2;
    }
  }
  if (INFO.value == 0) {
    if (M <= 0) {
      INFO.value = -3;
    } else if (N <= 0) {
      INFO.value = -4;
    } else if (LDA < max(1, M)) {
      INFO.value = -6;
    } else if (LDB < max(1, N)) {
      INFO.value = -8;
    } else if (LDC < max(1, M)) {
      INFO.value = -10;
    } else if (LDD < max(1, M)) {
      INFO.value = -12;
    } else if (LDE < max(1, N)) {
      INFO.value = -14;
    } else if (LDF < max(1, M)) {
      INFO.value = -16;
    }
  }

  if (INFO.value == 0) {
    if (NOTRAN) {
      if (IJOB == 1 || IJOB == 2) {
        LWMIN = max(1, 2 * M * N);
      } else {
        LWMIN = 1;
      }
    } else {
      LWMIN = 1;
    }
    WORK[1] = LWMIN.toDouble();

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -20;
    }
  }

  if (INFO.value != 0) {
    xerbla('DTGSYL', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    SCALE.value = 1;
    if (NOTRAN) {
      if (IJOB != 0) {
        DIF.value = 0;
      }
    }
    return;
  }

  // Determine optimal block sizes MB and NB

  MB = ilaenv(2, 'DTGSYL', TRANS, M, N, -1, -1);
  NB = ilaenv(5, 'DTGSYL', TRANS, M, N, -1, -1);

  ISOLVE = 1;
  IFUNC = 0;
  if (NOTRAN) {
    if (IJOB >= 3) {
      IFUNC = IJOB - 2;
      dlaset('F', M, N, ZERO, ZERO, C, LDC);
      dlaset('F', M, N, ZERO, ZERO, F, LDF);
    } else if (IJOB >= 1) {
      ISOLVE = 2;
    }
  }

  if ((MB <= 1 && NB <= 1) || (MB >= M && NB >= N)) {
    for (IROUND = 1; IROUND <= ISOLVE; IROUND++) {
      // Use unblocked Level 2 solver

      DSCALE.value = ZERO;
      DSUM.value = ONE;
      PQ.value = 0;
      dtgsy2(TRANS, IFUNC, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF,
          SCALE, DSUM, DSCALE, IWORK, PQ, INFO);
      if (DSCALE.value != ZERO) {
        if (IJOB == 1 || IJOB == 3) {
          DIF.value = sqrt(2 * M * N) / (DSCALE.value * sqrt(DSUM.value));
        } else {
          DIF.value = sqrt(PQ.value) / (DSCALE.value * sqrt(DSUM.value));
        }
      }

      if (ISOLVE == 2 && IROUND == 1) {
        if (NOTRAN) {
          IFUNC = IJOB;
        }
        SCALE2 = SCALE.value;
        dlacpy('F', M, N, C, LDC, WORK.asMatrix(M), M);
        dlacpy('F', M, N, F, LDF, WORK(M * N + 1).asMatrix(M), M);
        dlaset('F', M, N, ZERO, ZERO, C, LDC);
        dlaset('F', M, N, ZERO, ZERO, F, LDF);
      } else if (ISOLVE == 2 && IROUND == 2) {
        dlacpy('F', M, N, WORK.asMatrix(M), M, C, LDC);
        dlacpy('F', M, N, WORK(M * N + 1).asMatrix(M), M, F, LDF);
        SCALE.value = SCALE2;
      }
    }

    return;
  }

  // Determine block structure of A

  P = 0;
  I = 1;
  while (true) {
    if (I > M) break;
    P++;
    IWORK[P] = I;
    I += MB;
    if (I >= M) break;
    if (A[I][I - 1] != ZERO) I++;
  }

  IWORK[P + 1] = M + 1;
  if (IWORK[P] == IWORK[P + 1]) P--;

  // Determine block structure of B

  Q = P + 1;
  J = 1;
  while (true) {
    if (J > N) break;
    Q++;
    IWORK[Q] = J;
    J += NB;
    if (J >= N) break;
    if (B[J][J - 1] != ZERO) J++;
  }

  IWORK[Q + 1] = N + 1;
  if (IWORK[Q] == IWORK[Q + 1]) Q--;

  if (NOTRAN) {
    for (IROUND = 1; IROUND <= ISOLVE; IROUND++) {
      // Solve (I, J)-subsystem
      // A[I][ I] * R(I, J) - L(I, J) * B[J][ J] = C[I][ J]
      // D[I][ I] * R(I, J) - L(I, J) * E[J][ J] = F[I][ J]
      // for I = P, P - 1,..., 1; J = 1, 2,..., Q

      DSCALE.value = ZERO;
      DSUM.value = ONE;
      PQ.value = 0;
      SCALE.value = ONE;
      for (J = P + 2; J <= Q; J++) {
        JS = IWORK[J];
        JE = IWORK[J + 1] - 1;
        NB = JE - JS + 1;
        for (I = P; I >= 1; I--) {
          IS = IWORK[I];
          IE = IWORK[I + 1] - 1;
          MB = IE - IS + 1;
          PPQQ.value = 0;
          dtgsy2(
              TRANS,
              IFUNC,
              MB,
              NB,
              A(IS, IS),
              LDA,
              B(JS, JS),
              LDB,
              C(IS, JS),
              LDC,
              D(IS, IS),
              LDD,
              E(JS, JS),
              LDE,
              F(IS, JS),
              LDF,
              SCALOC,
              DSUM,
              DSCALE,
              IWORK(Q + 2),
              PPQQ,
              LINFO);
          if (LINFO.value > 0) INFO.value = LINFO.value;

          PQ.value += PPQQ.value;
          if (SCALOC.value != ONE) {
            for (K = 1; K <= JS - 1; K++) {
              dscal(M, SCALOC.value, C(1, K).asArray(), 1);
              dscal(M, SCALOC.value, F(1, K).asArray(), 1);
            }
            for (K = JS; K <= JE; K++) {
              dscal(IS - 1, SCALOC.value, C(1, K).asArray(), 1);
              dscal(IS - 1, SCALOC.value, F(1, K).asArray(), 1);
            }
            for (K = JS; K <= JE; K++) {
              dscal(M - IE, SCALOC.value, C(IE + 1, K).asArray(), 1);
              dscal(M - IE, SCALOC.value, F(IE + 1, K).asArray(), 1);
            }
            for (K = JE + 1; K <= N; K++) {
              dscal(M, SCALOC.value, C(1, K).asArray(), 1);
              dscal(M, SCALOC.value, F(1, K).asArray(), 1);
            }
            SCALE.value *= SCALOC.value;
          }

          // Substitute R(I, J) and L(I, J) into remaining
          // equation.

          if (I > 1) {
            dgemm('N', 'N', IS - 1, NB, MB, -ONE, A(1, IS), LDA, C(IS, JS), LDC,
                ONE, C(1, JS), LDC);
            dgemm('N', 'N', IS - 1, NB, MB, -ONE, D(1, IS), LDD, C(IS, JS), LDC,
                ONE, F(1, JS), LDF);
          }
          if (J < Q) {
            dgemm('N', 'N', MB, N - JE, NB, ONE, F(IS, JS), LDF, B(JS, JE + 1),
                LDB, ONE, C(IS, JE + 1), LDC);
            dgemm('N', 'N', MB, N - JE, NB, ONE, F(IS, JS), LDF, E(JS, JE + 1),
                LDE, ONE, F(IS, JE + 1), LDF);
          }
        }
      }
      if (DSCALE.value != ZERO) {
        if (IJOB == 1 || IJOB == 3) {
          DIF.value = sqrt(2 * M * N) / (DSCALE.value * sqrt(DSUM.value));
        } else {
          DIF.value = sqrt(PQ.value) / (DSCALE.value * sqrt(DSUM.value));
        }
      }
      if (ISOLVE == 2 && IROUND == 1) {
        if (NOTRAN) {
          IFUNC = IJOB;
        }
        SCALE2 = SCALE.value;
        dlacpy('F', M, N, C, LDC, WORK.asMatrix(M), M);
        dlacpy('F', M, N, F, LDF, WORK(M * N + 1).asMatrix(M), M);
        dlaset('F', M, N, ZERO, ZERO, C, LDC);
        dlaset('F', M, N, ZERO, ZERO, F, LDF);
      } else if (ISOLVE == 2 && IROUND == 2) {
        dlacpy('F', M, N, WORK.asMatrix(M), M, C, LDC);
        dlacpy('F', M, N, WORK(M * N + 1).asMatrix(M), M, F, LDF);
        SCALE.value = SCALE2;
      }
    }
  } else {
    // Solve transposed (I, J)-subsystem
    // A[I][ I]**T * R(I, J)  + D[I][ I]**T * L(I, J)  =  C[I][ J]
    // R(I, J)  * B[J][ J]**T + L(I, J)  * E[J][ J]**T = -F[I][ J]
    // for I = 1,2,..., P; J = Q, Q-1,..., 1

    SCALE.value = ONE;
    for (I = 1; I <= P; I++) {
      IS = IWORK[I];
      IE = IWORK[I + 1] - 1;
      MB = IE - IS + 1;
      for (J = Q; J >= P + 2; J--) {
        JS = IWORK[J];
        JE = IWORK[J + 1] - 1;
        NB = JE - JS + 1;
        dtgsy2(
            TRANS,
            IFUNC,
            MB,
            NB,
            A(IS, IS),
            LDA,
            B(JS, JS),
            LDB,
            C(IS, JS),
            LDC,
            D(IS, IS),
            LDD,
            E(JS, JS),
            LDE,
            F(IS, JS),
            LDF,
            SCALOC,
            DSUM,
            DSCALE,
            IWORK(Q + 2),
            PPQQ,
            LINFO);
        if (LINFO.value > 0) INFO.value = LINFO.value;
        if (SCALOC.value != ONE) {
          for (K = 1; K <= JS - 1; K++) {
            dscal(M, SCALOC.value, C(1, K).asArray(), 1);
            dscal(M, SCALOC.value, F(1, K).asArray(), 1);
          }
          for (K = JS; K <= JE; K++) {
            dscal(IS - 1, SCALOC.value, C(1, K).asArray(), 1);
            dscal(IS - 1, SCALOC.value, F(1, K).asArray(), 1);
          }
          for (K = JS; K <= JE; K++) {
            dscal(M - IE, SCALOC.value, C(IE + 1, K).asArray(), 1);
            dscal(M - IE, SCALOC.value, F(IE + 1, K).asArray(), 1);
          }
          for (K = JE + 1; K <= N; K++) {
            dscal(M, SCALOC.value, C(1, K).asArray(), 1);
            dscal(M, SCALOC.value, F(1, K).asArray(), 1);
          }
          SCALE.value *= SCALOC.value;
        }

        // Substitute R(I, J) and L(I, J) into remaining equation.

        if (J > P + 2) {
          dgemm('N', 'T', MB, JS - 1, NB, ONE, C(IS, JS), LDC, B(1, JS), LDB,
              ONE, F(IS, 1), LDF);
          dgemm('N', 'T', MB, JS - 1, NB, ONE, F(IS, JS), LDF, E(1, JS), LDE,
              ONE, F(IS, 1), LDF);
        }
        if (I < P) {
          dgemm('T', 'N', M - IE, NB, MB, -ONE, A(IS, IE + 1), LDA, C(IS, JS),
              LDC, ONE, C(IE + 1, JS), LDC);
          dgemm('T', 'N', M - IE, NB, MB, -ONE, D(IS, IE + 1), LDD, F(IS, JS),
              LDF, ONE, C(IE + 1, JS), LDC);
        }
      }
    }
  }

  WORK[1] = LWMIN.toDouble();
}
