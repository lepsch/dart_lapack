// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/drot.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgeqr2.dart';
import 'package:dart_lapack/src/dgerq2.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlagv2.dart';
import 'package:dart_lapack/src/dlartg.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dlassq.dart';
import 'package:dart_lapack/src/dorg2r.dart';
import 'package:dart_lapack/src/dorgr2.dart';
import 'package:dart_lapack/src/dorm2r.dart';
import 'package:dart_lapack/src/dormr2.dart';
import 'package:dart_lapack/src/dtgsy2.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void dtgex2(
  final bool WANTQ,
  final bool WANTZ,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final int J1,
  final int N1,
  final int N2,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const TWENTY = 2.0e+01;
  const LDST = 4;
  const WANDS = true;
  bool STRONG, WEAK;
  int I, M;
  double BQRA21,
      BRQA21,
      DNORMA,
      DNORMB,
      EPS,
      F,
      G,
      SA,
      SB,
      SMLNUM,
      THRESHA,
      THRESHB;
  final DSCALE = Box(0.0), DSUM = Box(0.0), DDUM = Box(0.0), SCALE = Box(0.0);
  final LINFO = Box(0), IDUM = Box(0);
  final IWORK = Array<int>(LDST + 2);
  final AI = Array<double>(2),
      AR = Array<double>(2),
      TAUL = Array<double>(LDST),
      TAUR = Array<double>(LDST),
      BE = Array<double>(2);
  final IR = Matrix<double>(LDST, LDST),
      IRCOP = Matrix<double>(LDST, LDST),
      LI = Matrix<double>(LDST, LDST),
      LICOP = Matrix<double>(LDST, LDST),
      S = Matrix<double>(LDST, LDST),
      SCPY = Matrix<double>(LDST, LDST),
      T = Matrix<double>(LDST, LDST),
      TCPY = Matrix<double>(LDST, LDST);

  INFO.value = 0;

  // Quick return if possible
  if (N <= 1 || N1 <= 0 || N2 <= 0) return;
  if (N1 > N || (J1 + N1) > N) return;
  M = N1 + N2;
  if (LWORK < max(1, max(N * M, M * M * 2))) {
    INFO.value = -16;
    WORK[1] = max(1, max(N * M, M * M * 2)).toDouble();
    return;
  }

  WEAK = false;
  STRONG = false;

  // Make a local copy of selected block
  dlaset('Full', LDST, LDST, ZERO, ZERO, LI, LDST);
  dlaset('Full', LDST, LDST, ZERO, ZERO, IR, LDST);
  dlacpy('Full', M, M, A(J1, J1), LDA, S, LDST);
  dlacpy('Full', M, M, B(J1, J1), LDB, T, LDST);

  // Compute threshold for testing acceptance of swapping.
  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  DSCALE.value = ZERO;
  DSUM.value = ONE;
  dlacpy('Full', M, M, S, LDST, WORK.asMatrix(M), M);
  dlassq(M * M, WORK, 1, DSCALE, DSUM);
  DNORMA = DSCALE.value * sqrt(DSUM.value);
  DSCALE.value = ZERO;
  DSUM.value = ONE;
  dlacpy('Full', M, M, T, LDST, WORK.asMatrix(M), M);
  dlassq(M * M, WORK, 1, DSCALE, DSUM);
  DNORMB = DSCALE.value * sqrt(DSUM.value);
  THRESHA = max(TWENTY * EPS * DNORMA, SMLNUM);
  THRESHB = max(TWENTY * EPS * DNORMB, SMLNUM);

  if (M == 2) {
    // CASE 1: Swap 1-by-1 and 1-by-1 blocks.

    // Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks
    // using Givens rotations and perform the swap tentatively.
    F = S[2][2] * T[1][1] - T[2][2] * S[1][1];
    G = S[2][2] * T[1][2] - T[2][2] * S[1][2];
    SA = S[2][2].abs() * T[1][1].abs();
    SB = S[1][1].abs() * T[2][2].abs();
    dlartg(F, G, IR.box(1, 2), IR.box(1, 1), DDUM);
    IR[2][1] = -IR[1][2];
    IR[2][2] = IR[1][1];
    drot(2, S(1, 1).asArray(), 1, S(1, 2).asArray(), 1, IR[1][1], IR[2][1]);
    drot(2, T(1, 1).asArray(), 1, T(1, 2).asArray(), 1, IR[1][1], IR[2][1]);
    if (SA >= SB) {
      dlartg(S[1][1], S[2][1], LI.box(1, 1), LI.box(2, 1), DDUM);
    } else {
      dlartg(T[1][1], T[2][1], LI.box(1, 1), LI.box(2, 1), DDUM);
    }
    drot(2, S(1, 1).asArray(), LDST, S(2, 1).asArray(), LDST, LI[1][1],
        LI[2][1]);
    drot(2, T(1, 1).asArray(), LDST, T(2, 1).asArray(), LDST, LI[1][1],
        LI[2][1]);
    LI[2][2] = LI[1][1];
    LI[1][2] = -LI[2][1];

    // Weak stability test: |S21| <= O(EPS F-norm((A)))
    //                 and  |T21| <= O(EPS F-norm((B)))
    WEAK = S[2][1].abs() <= THRESHA && T[2][1].abs() <= THRESHB;
    if (!WEAK) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }

    if (WANDS) {
      // Strong stability test:
      //     F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
      //     and
      //     F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
      dlacpy('Full', M, M, A(J1, J1), LDA, WORK(M * M + 1).asMatrix(M), M);
      dgemm(
          'N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK.asMatrix(M), M);
      dgemm('N', 'T', M, M, M, -ONE, WORK.asMatrix(M), M, IR, LDST, ONE,
          WORK(M * M + 1).asMatrix(M), M);
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      dlassq(M * M, WORK(M * M + 1), 1, DSCALE, DSUM);
      SA = DSCALE.value * sqrt(DSUM.value);

      dlacpy('Full', M, M, B(J1, J1), LDB, WORK(M * M + 1).asMatrix(M), M);
      dgemm(
          'N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK.asMatrix(M), M);
      dgemm('N', 'T', M, M, M, -ONE, WORK.asMatrix(M), M, IR, LDST, ONE,
          WORK(M * M + 1).asMatrix(M), M);
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      dlassq(M * M, WORK(M * M + 1), 1, DSCALE, DSUM);
      SB = DSCALE.value * sqrt(DSUM.value);
      STRONG = SA <= THRESHA && SB <= THRESHB;
      if (!STRONG) {
        // Exit with INFO = 1 if swap was rejected.
        INFO.value = 1;
        return;
      }
    }

    // Update (A[J1:J1+M-1][ M+J1:N], B[J1:J1+M-1][ M+J1:N]) and
    //        (A[1:J1-1][ J1:J1+M], B[1:J1-1][ J1:J1+M]).
    drot(J1 + 1, A(1, J1).asArray(), 1, A(1, J1 + 1).asArray(), 1, IR[1][1],
        IR[2][1]);
    drot(J1 + 1, B(1, J1).asArray(), 1, B(1, J1 + 1).asArray(), 1, IR[1][1],
        IR[2][1]);
    drot(N - J1 + 1, A(J1, J1).asArray(), LDA, A(J1 + 1, J1).asArray(), LDA,
        LI[1][1], LI[2][1]);
    drot(N - J1 + 1, B(J1, J1).asArray(), LDB, B(J1 + 1, J1).asArray(), LDB,
        LI[1][1], LI[2][1]);

    // Set  N1-by-N2 (2,1) - blocks to ZERO.
    A[J1 + 1][J1] = ZERO;
    B[J1 + 1][J1] = ZERO;

    // Accumulate transformations into Q and Z if requested.
    if (WANTZ) {
      drot(N, Z(1, J1).asArray(), 1, Z(1, J1 + 1).asArray(), 1, IR[1][1],
          IR[2][1]);
    }
    if (WANTQ) {
      drot(N, Q(1, J1).asArray(), 1, Q(1, J1 + 1).asArray(), 1, LI[1][1],
          LI[2][1]);
    }

    // Exit with INFO = 0 if swap was successfully performed.
    return;
  } else {
    // CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
    //         and 2-by-2 blocks.

    // Solve the generalized Sylvester equation
    //          S11 * R - L * S22 = SCALE * S12
    //          T11 * R - L * T22 = SCALE * T12
    // for R and L. Solutions in LI and IR.
    dlacpy('Full', N1, N2, T(1, N1 + 1), LDST, LI, LDST);
    dlacpy('Full', N1, N2, S(1, N1 + 1), LDST, IR(N2 + 1, N1 + 1), LDST);
    dtgsy2(
        'N',
        0,
        N1,
        N2,
        S,
        LDST,
        S(N1 + 1, N1 + 1),
        LDST,
        IR(N2 + 1, N1 + 1),
        LDST,
        T,
        LDST,
        T(N1 + 1, N1 + 1),
        LDST,
        LI,
        LDST,
        SCALE,
        DSUM,
        DSCALE,
        IWORK,
        IDUM,
        LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }

    // Compute orthogonal matrix QL:
    //
    //             QL**T * LI = [ TL ]
    //                          [ 0  ]
    // where
    //             LI =  [      -L              ]
    //                   [ SCALE * identity(N2) ]
    for (I = 1; I <= N2; I++) {
      dscal(N1, -ONE, LI(1, I).asArray(), 1);
      LI[N1 + I][I] = SCALE.value;
    }
    dgeqr2(M, N2, LI, LDST, TAUL, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }
    dorg2r(M, M, N2, LI, LDST, TAUL, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }

    // Compute orthogonal matrix RQ:
    // IR * RQ**T =   [ 0  TR],
    // where IR = [ SCALE * identity(N1), R ]

    for (I = 1; I <= N1; I++) {
      IR[N2 + I][I] = SCALE.value;
    }
    dgerq2(N1, M, IR(N2 + 1, 1), LDST, TAUR, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }
    dorgr2(M, M, N1, IR, LDST, TAUR, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }

    // Perform the swapping tentatively:
    dgemm('T', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK.asMatrix(M), M);
    dgemm('N', 'T', M, M, M, ONE, WORK.asMatrix(M), M, IR, LDST, ZERO, S, LDST);
    dgemm('T', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK.asMatrix(M), M);
    dgemm('N', 'T', M, M, M, ONE, WORK.asMatrix(M), M, IR, LDST, ZERO, T, LDST);
    dlacpy('F', M, M, S, LDST, SCPY, LDST);
    dlacpy('F', M, M, T, LDST, TCPY, LDST);
    dlacpy('F', M, M, IR, LDST, IRCOP, LDST);
    dlacpy('F', M, M, LI, LDST, LICOP, LDST);

    // Triangularize the B-part by an RQ factorization.
    // Apply transformation (from left) to A-part, giving S.
    dgerq2(M, M, T, LDST, TAUR, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }
    dormr2('R', 'T', M, M, M, T, LDST, TAUR, S, LDST, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }
    dormr2('L', 'N', M, M, M, T, LDST, TAUR, IR, LDST, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }

    // Compute F-norm(S21) in BRQA21. (T21 is 0.)
    DSCALE.value = ZERO;
    DSUM.value = ONE;
    for (I = 1; I <= N2; I++) {
      dlassq(N1, S(N2 + 1, I).asArray(), 1, DSCALE, DSUM);
    }
    BRQA21 = DSCALE.value * sqrt(DSUM.value);

    // Triangularize the B-part by a QR factorization.
    // Apply transformation (from right) to A-part, giving S.
    dgeqr2(M, M, TCPY, LDST, TAUL, WORK, LINFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }
    dorm2r('L', 'T', M, M, M, TCPY, LDST, TAUL, SCPY, LDST, WORK, INFO);
    dorm2r('R', 'N', M, M, M, TCPY, LDST, TAUL, LICOP, LDST, WORK, INFO);
    if (LINFO.value != 0) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }

    // Compute F-norm(S21) in BQRA21. (T21 is 0.)
    DSCALE.value = ZERO;
    DSUM.value = ONE;
    for (I = 1; I <= N2; I++) {
      dlassq(N1, SCPY(N2 + 1, I).asArray(), 1, DSCALE, DSUM);
    }
    BQRA21 = DSCALE.value * sqrt(DSUM.value);

    // Decide which method to use.
    //   Weak stability test:
    //      F-norm(S21) <= O(EPS * F-norm((S)))
    if (BQRA21 <= BRQA21 && BQRA21 <= THRESHA) {
      dlacpy('F', M, M, SCPY, LDST, S, LDST);
      dlacpy('F', M, M, TCPY, LDST, T, LDST);
      dlacpy('F', M, M, IRCOP, LDST, IR, LDST);
      dlacpy('F', M, M, LICOP, LDST, LI, LDST);
    } else if (BRQA21 >= THRESHA) {
      // Exit with INFO = 1 if swap was rejected.
      INFO.value = 1;
      return;
    }

    // Set lower triangle of B-part to zero
    dlaset('Lower', M - 1, M - 1, ZERO, ZERO, T(2, 1), LDST);

    if (WANDS) {
      // Strong stability test:
      //     F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
      //     and
      //     F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
      dlacpy('Full', M, M, A(J1, J1), LDA, WORK(M * M + 1).asMatrix(M), M);
      dgemm(
          'N', 'N', M, M, M, ONE, LI, LDST, S, LDST, ZERO, WORK.asMatrix(M), M);
      dgemm('N', 'N', M, M, M, -ONE, WORK.asMatrix(M), M, IR, LDST, ONE,
          WORK(M * M + 1).asMatrix(M), M);
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      dlassq(M * M, WORK(M * M + 1), 1, DSCALE, DSUM);
      SA = DSCALE.value * sqrt(DSUM.value);

      dlacpy('Full', M, M, B(J1, J1), LDB, WORK(M * M + 1).asMatrix(M), M);
      dgemm(
          'N', 'N', M, M, M, ONE, LI, LDST, T, LDST, ZERO, WORK.asMatrix(M), M);
      dgemm('N', 'N', M, M, M, -ONE, WORK.asMatrix(M), M, IR, LDST, ONE,
          WORK(M * M + 1).asMatrix(M), M);
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      dlassq(M * M, WORK(M * M + 1), 1, DSCALE, DSUM);
      SB = DSCALE.value * sqrt(DSUM.value);
      STRONG = SA <= THRESHA && SB <= THRESHB;
      if (!STRONG) {
        // Exit with INFO = 1 if swap was rejected.
        INFO.value = 1;
        return;
      }
    }

    // If the swap is accepted ("weakly" and "strongly"), apply the
    // transformations and set N1-by-N2 (2,1)-block to zero.
    dlaset('Full', N1, N2, ZERO, ZERO, S(N2 + 1, 1), LDST);

    // copy back M-by-M diagonal block starting at index J1 of (A, B)
    dlacpy('F', M, M, S, LDST, A(J1, J1), LDA);
    dlacpy('F', M, M, T, LDST, B(J1, J1), LDB);
    dlaset('Full', LDST, LDST, ZERO, ZERO, T, LDST);

    // Standardize existing 2-by-2 blocks.
    dlaset('Full', M, M, ZERO, ZERO, WORK.asMatrix(M), M);
    WORK[1] = ONE;
    T[1][1] = ONE;
    IDUM.value = LWORK - M * M - 2;
    if (N2 > 1) {
      dlagv2(A(J1, J1), LDA, B(J1, J1), LDB, AR, AI, BE, WORK.box(1),
          WORK.box(2), T.box(1, 1), T.box(2, 1));
      WORK[M + 1] = -WORK[2];
      WORK[M + 2] = WORK[1];
      T[N2][N2] = T[1][1];
      T[1][2] = -T[2][1];
    }
    WORK[M * M] = ONE;
    T[M][M] = ONE;

    if (N1 > 1) {
      dlagv2(
        A(J1 + N2, J1 + N2),
        LDA,
        B(J1 + N2, J1 + N2),
        LDB,
        TAUR,
        TAUL,
        WORK(M * M + 1),
        WORK.box(N2 * M + N2 + 1),
        WORK.box(N2 * M + N2 + 2),
        T.box(N2 + 1, N2 + 1),
        T.box(M, M - 1),
      );
      WORK[M * M] = WORK[N2 * M + N2 + 1];
      WORK[M * M - 1] = -WORK[N2 * M + N2 + 2];
      T[M][M] = T[N2 + 1][N2 + 1];
      T[M - 1][M] = -T[M][M - 1];
    }
    dgemm('T', 'N', N2, N1, N2, ONE, WORK.asMatrix(M), M, A(J1, J1 + N2), LDA,
        ZERO, WORK(M * M + 1).asMatrix(N2), N2);
    dlacpy(
        'Full', N2, N1, WORK(M * M + 1).asMatrix(N2), N2, A(J1, J1 + N2), LDA);
    dgemm('T', 'N', N2, N1, N2, ONE, WORK.asMatrix(M), M, B(J1, J1 + N2), LDB,
        ZERO, WORK(M * M + 1).asMatrix(N2), N2);
    dlacpy(
        'Full', N2, N1, WORK(M * M + 1).asMatrix(N2), N2, B(J1, J1 + N2), LDB);
    dgemm('N', 'N', M, M, M, ONE, LI, LDST, WORK.asMatrix(M), M, ZERO,
        WORK(M * M + 1).asMatrix(M), M);
    dlacpy('Full', M, M, WORK(M * M + 1).asMatrix(M), M, LI, LDST);
    dgemm('N', 'N', N2, N1, N1, ONE, A(J1, J1 + N2), LDA, T(N2 + 1, N2 + 1),
        LDST, ZERO, WORK.asMatrix(N2), N2);
    dlacpy('Full', N2, N1, WORK.asMatrix(N2), N2, A(J1, J1 + N2), LDA);
    dgemm('N', 'N', N2, N1, N1, ONE, B(J1, J1 + N2), LDB, T(N2 + 1, N2 + 1),
        LDST, ZERO, WORK.asMatrix(N2), N2);
    dlacpy('Full', N2, N1, WORK.asMatrix(N2), N2, B(J1, J1 + N2), LDB);
    dgemm('T', 'N', M, M, M, ONE, IR, LDST, T, LDST, ZERO, WORK.asMatrix(M), M);
    dlacpy('Full', M, M, WORK.asMatrix(M), M, IR, LDST);

    // Accumulate transformations into Q and Z if requested.
    if (WANTQ) {
      dgemm('N', 'N', N, M, M, ONE, Q(1, J1), LDQ, LI, LDST, ZERO,
          WORK.asMatrix(N), N);
      dlacpy('Full', N, M, WORK.asMatrix(N), N, Q(1, J1), LDQ);
    }

    if (WANTZ) {
      dgemm('N', 'N', N, M, M, ONE, Z(1, J1), LDZ, IR, LDST, ZERO,
          WORK.asMatrix(N), N);
      dlacpy('Full', N, M, WORK.asMatrix(N), N, Z(1, J1), LDZ);
    }

    // Update (A[J1:J1+M-1][ M+J1:N], B[J1:J1+M-1][ M+J1:N]) and
    // (A[1:J1-1][ J1:J1+M], B[1:J1-1][ J1:J1+M]).
    I = J1 + M;
    if (I <= N) {
      dgemm('T', 'N', M, N - I + 1, M, ONE, LI, LDST, A(J1, I), LDA, ZERO,
          WORK.asMatrix(M), M);
      dlacpy('Full', M, N - I + 1, WORK.asMatrix(M), M, A(J1, I), LDA);
      dgemm('T', 'N', M, N - I + 1, M, ONE, LI, LDST, B(J1, I), LDB, ZERO,
          WORK.asMatrix(M), M);
      dlacpy('Full', M, N - I + 1, WORK.asMatrix(M), M, B(J1, I), LDB);
    }
    I = J1 - 1;
    if (I > 0) {
      dgemm('N', 'N', I, M, M, ONE, A(1, J1), LDA, IR, LDST, ZERO,
          WORK.asMatrix(I), I);
      dlacpy('Full', I, M, WORK.asMatrix(I), I, A(1, J1), LDA);
      dgemm('N', 'N', I, M, M, ONE, B(1, J1), LDB, IR, LDST, ZERO,
          WORK.asMatrix(I), I);
      dlacpy('Full', I, M, WORK.asMatrix(I), I, B(1, J1), LDB);
    }

    // Exit with INFO = 0 if swap was successfully performed.
  }
}
