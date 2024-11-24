// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlalsa.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasda.dart';
import 'package:lapack/src/dlasdq.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dlasrt.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlalsd(
  final String UPLO,
  final int SMLSIZ,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> B_,
  final int LDB,
  final double RCOND,
  final Box<int> RANK,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having(length: LDB);
  final E = E_.having(length: LDB);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  int BX,
      BXST,
      C,
      DIFL,
      DIFR,
      GIVCOL,
      GIVNUM,
      GIVPTR,
      I,
      ICMPQ1,
      ICMPQ2,
      IWK,
      J,
      K,
      NLVL,
      NM1,
      NSIZE,
      NSUB,
      NWORK,
      PERM,
      POLES,
      S,
      SIZEI,
      SMLSZP,
      SQRE,
      ST,
      ST1,
      U,
      VT,
      Z;
  double EPS, ORGNRM, RCND, TOL;
  final CS = Box(0.0), R = Box(0.0), SN = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;

  if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 1) {
    INFO.value = -4;
  } else if ((LDB < 1) || (LDB < N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('DLALSD', -INFO.value);
    return;
  }

  EPS = dlamch('Epsilon');

  // Set up the tolerance.

  if ((RCOND <= ZERO) || (RCOND >= ONE)) {
    RCND = EPS;
  } else {
    RCND = RCOND;
  }

  RANK.value = 0;

  // Quick return if possible.

  if (N == 0) {
    return;
  } else if (N == 1) {
    if (D[1] == ZERO) {
      dlaset('A', 1, NRHS, ZERO, ZERO, B, LDB);
    } else {
      RANK.value = 1;
      dlascl('G', 0, 0, D[1], ONE, 1, NRHS, B, LDB, INFO);
      D[1] = D[1].abs();
    }
    return;
  }

  // Rotate the matrix if it is lower bidiagonal.

  if (UPLO == 'L') {
    for (I = 1; I <= N - 1; I++) {
      dlartg(D[I], E[I], CS, SN, R);
      D[I] = R.value;
      E[I] = SN.value * D[I + 1];
      D[I + 1] = CS.value * D[I + 1];
      if (NRHS == 1) {
        drot(1, B(I, 1).asArray(), 1, B(I + 1, 1).asArray(), 1, CS.value,
            SN.value);
      } else {
        WORK[I * 2 - 1] = CS.value;
        WORK[I * 2] = SN.value;
      }
    }
    if (NRHS > 1) {
      for (I = 1; I <= NRHS; I++) {
        for (J = 1; J <= N - 1; J++) {
          CS.value = WORK[J * 2 - 1];
          SN.value = WORK[J * 2];
          drot(1, B(J, I).asArray(), 1, B(J + 1, I).asArray(), 1, CS.value,
              SN.value);
        }
      }
    }
  }

  // Scale.

  NM1 = N - 1;
  ORGNRM = dlanst('M', N, D, E);
  if (ORGNRM == ZERO) {
    dlaset('A', N, NRHS, ZERO, ZERO, B, LDB);
    return;
  }

  dlascl('G', 0, 0, ORGNRM, ONE, N, 1, D.asMatrix(N), N, INFO);
  dlascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E.asMatrix(NM1), NM1, INFO);

  // If N is smaller than the minimum divide size SMLSIZ, then solve
  // the problem with another solver.

  if (N <= SMLSIZ) {
    NWORK = 1 + N * N;
    dlaset('A', N, N, ZERO, ONE, WORK.asMatrix(N), N);
    dlasdq('U', 0, N, N, 0, NRHS, D, E, WORK.asMatrix(N), N, WORK.asMatrix(N),
        N, B, LDB, WORK(NWORK), INFO);
    if (INFO.value != 0) {
      return;
    }
    TOL = RCND * D[idamax(N, D, 1)].abs();
    for (I = 1; I <= N; I++) {
      if (D[I] <= TOL) {
        dlaset('A', 1, NRHS, ZERO, ZERO, B(I, 1), LDB);
      } else {
        dlascl('G', 0, 0, D[I], ONE, 1, NRHS, B(I, 1), LDB, INFO);
        RANK.value++;
      }
    }
    dgemm('T', 'N', N, NRHS, N, ONE, WORK.asMatrix(N), N, B, LDB, ZERO,
        WORK(NWORK).asMatrix(N), N);
    dlacpy('A', N, NRHS, WORK(NWORK).asMatrix(N), N, B, LDB);

    // Unscale.

    dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D.asMatrix(N), N, INFO);
    dlasrt('D', N, D, INFO);
    dlascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO);

    return;
  }

  // Book-keeping and setting up some constants.

  NLVL = log(N / (SMLSIZ + 1)) ~/ log(TWO) + 1;

  SMLSZP = SMLSIZ + 1;

  U = 1;
  VT = 1 + SMLSIZ * N;
  DIFL = VT + SMLSZP * N;
  DIFR = DIFL + NLVL * N;
  Z = DIFR + NLVL * N * 2;
  C = Z + NLVL * N;
  S = C + N;
  POLES = S + N;
  GIVNUM = POLES + 2 * NLVL * N;
  BX = GIVNUM + 2 * NLVL * N;
  NWORK = BX + N * NRHS;

  SIZEI = 1 + N;
  K = SIZEI + N;
  GIVPTR = K + N;
  PERM = GIVPTR + N;
  GIVCOL = PERM + NLVL * N;
  IWK = GIVCOL + NLVL * N * 2;

  ST = 1;
  SQRE = 0;
  ICMPQ1 = 1;
  ICMPQ2 = 0;
  NSUB = 0;

  for (I = 1; I <= N; I++) {
    if (D[I].abs() < EPS) {
      D[I] = sign(EPS, D[I]);
    }
  }

  for (I = 1; I <= NM1; I++) {
    if ((E[I].abs() < EPS) || (I == NM1)) {
      NSUB++;
      IWORK[NSUB] = ST;

      // Subproblem found. First determine its size and then
      // apply divide and conquer on it.

      if (I < NM1) {
        // A subproblem with E[I] small for I < NM1.

        NSIZE = I - ST + 1;
        IWORK[SIZEI + NSUB - 1] = NSIZE;
      } else if (E[I].abs() >= EPS) {
        // A subproblem with E[NM1] not too small but I = NM1.

        NSIZE = N - ST + 1;
        IWORK[SIZEI + NSUB - 1] = NSIZE;
      } else {
        // A subproblem with E[NM1] small. This implies an
        // 1-by-1 subproblem at D[N], which is not solved
        // explicitly.

        NSIZE = I - ST + 1;
        IWORK[SIZEI + NSUB - 1] = NSIZE;
        NSUB++;
        IWORK[NSUB] = N;
        IWORK[SIZEI + NSUB - 1] = 1;
        dcopy(NRHS, B(N, 1).asArray(), LDB, WORK(BX + NM1), N);
      }
      ST1 = ST - 1;
      if (NSIZE == 1) {
        // This is a 1-by-1 subproblem and is not solved
        // explicitly.

        dcopy(NRHS, B(ST, 1).asArray(), LDB, WORK(BX + ST1), N);
      } else if (NSIZE <= SMLSIZ) {
        // This is a small subproblem and is solved by DLASDQ.

        dlaset('A', NSIZE, NSIZE, ZERO, ONE, WORK(VT + ST1).asMatrix(N), N);
        dlasdq(
            'U',
            0,
            NSIZE,
            NSIZE,
            0,
            NRHS,
            D(ST),
            E(ST),
            WORK(VT + ST1).asMatrix(N),
            N,
            WORK(NWORK).asMatrix(N),
            N,
            B(ST, 1),
            LDB,
            WORK(NWORK),
            INFO);
        if (INFO.value != 0) {
          return;
        }
        dlacpy('A', NSIZE, NRHS, B(ST, 1), LDB, WORK(BX + ST1).asMatrix(N), N);
      } else {
        // A large problem. Solve it using divide and conquer.

        dlasda(
            ICMPQ1,
            SMLSIZ,
            NSIZE,
            SQRE,
            D(ST),
            E(ST),
            WORK(U + ST1).asMatrix(N),
            N,
            WORK(VT + ST1).asMatrix(N),
            IWORK(K + ST1),
            WORK(DIFL + ST1).asMatrix(N),
            WORK(DIFR + ST1).asMatrix(N),
            WORK(Z + ST1).asMatrix(N),
            WORK(POLES + ST1).asMatrix(N),
            IWORK(GIVPTR + ST1),
            IWORK(GIVCOL + ST1).asMatrix(N),
            N,
            IWORK(PERM + ST1).asMatrix(N),
            WORK(GIVNUM + ST1).asMatrix(N),
            WORK(C + ST1),
            WORK(S + ST1),
            WORK(NWORK),
            IWORK(IWK),
            INFO);
        if (INFO.value != 0) {
          return;
        }
        BXST = BX + ST1;
        dlalsa(
            ICMPQ2,
            SMLSIZ,
            NSIZE,
            NRHS,
            B(ST, 1),
            LDB,
            WORK(BXST).asMatrix(N),
            N,
            WORK(U + ST1).asMatrix(N),
            N,
            WORK(VT + ST1).asMatrix(N),
            IWORK(K + ST1),
            WORK(DIFL + ST1).asMatrix(N),
            WORK(DIFR + ST1).asMatrix(N),
            WORK(Z + ST1).asMatrix(N),
            WORK(POLES + ST1).asMatrix(N),
            IWORK(GIVPTR + ST1),
            IWORK(GIVCOL + ST1).asMatrix(N),
            N,
            IWORK(PERM + ST1).asMatrix(N),
            WORK(GIVNUM + ST1).asMatrix(N),
            WORK(C + ST1),
            WORK(S + ST1),
            WORK(NWORK),
            IWORK(IWK),
            INFO);
        if (INFO.value != 0) {
          return;
        }
      }
      ST = I + 1;
    }
  }

  // Apply the singular values and treat the tiny ones as zero.

  TOL = RCND * D[idamax(N, D, 1)].abs();

  for (I = 1; I <= N; I++) {
    // Some of the elements in D can be negative because 1-by-1
    // subproblems were not solved explicitly.

    if (D[I].abs() <= TOL) {
      dlaset('A', 1, NRHS, ZERO, ZERO, WORK(BX + I - 1).asMatrix(N), N);
    } else {
      RANK.value++;
      dlascl(
          'G', 0, 0, D[I], ONE, 1, NRHS, WORK(BX + I - 1).asMatrix(N), N, INFO);
    }
    D[I] = D[I].abs();
  }

  // Now apply back the right singular vectors.

  ICMPQ2 = 1;
  for (I = 1; I <= NSUB; I++) {
    ST = IWORK[I];
    ST1 = ST - 1;
    NSIZE = IWORK[SIZEI + I - 1];
    BXST = BX + ST1;
    if (NSIZE == 1) {
      dcopy(NRHS, WORK(BXST), N, B(ST, 1).asArray(), LDB);
    } else if (NSIZE <= SMLSIZ) {
      dgemm('T', 'N', NSIZE, NRHS, NSIZE, ONE, WORK(VT + ST1).asMatrix(N), N,
          WORK(BXST).asMatrix(N), N, ZERO, B(ST, 1), LDB);
    } else {
      dlalsa(
          ICMPQ2,
          SMLSIZ,
          NSIZE,
          NRHS,
          WORK(BXST).asMatrix(N),
          N,
          B(ST, 1),
          LDB,
          WORK(U + ST1).asMatrix(N),
          N,
          WORK(VT + ST1).asMatrix(N),
          IWORK(K + ST1),
          WORK(DIFL + ST1).asMatrix(N),
          WORK(DIFR + ST1).asMatrix(N),
          WORK(Z + ST1).asMatrix(N),
          WORK(POLES + ST1).asMatrix(N),
          IWORK(GIVPTR + ST1),
          IWORK(GIVCOL + ST1).asMatrix(N),
          N,
          IWORK(PERM + ST1).asMatrix(N),
          WORK(GIVNUM + ST1).asMatrix(N),
          WORK(C + ST1),
          WORK(S + ST1),
          WORK(NWORK),
          IWORK(IWK),
          INFO);
      if (INFO.value != 0) {
        return;
      }
    }
  }

  // Unscale and sort the singular values.

  dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D.asMatrix(N), N, INFO);
  dlasrt('D', N, D, INFO);
  dlascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO);
}
