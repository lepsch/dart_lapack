// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dtrsyl.dart';
import 'package:lapack/src/dtrsyl3.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/range.dart';
import 'package:test/test.dart';

import '../matgen/dlatmr.dart';
import '../test_driver.dart';

void dsyl01(
  final double THRESH,
  final Array<int> NFAIL_,
  final Array<double> RMAX_,
  final Array<int> NINFO_,
  final Box<int> KNT,
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NFAIL = NFAIL_.having();
  final RMAX = RMAX_.having();
  final NINFO = NINFO_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXM = 245, MAXN = 192, LDSWORK = 36;
  final DUML = Array<double>(MAXM),
      DUMR = Array<double>(MAXN),
      D = Array<double>(max(MAXM, MAXN)),
      DUM = Array<double>(MAXN);
  final ISEED = Array<int>(4), IWORK = Array<int>(MAXM + MAXN + 2);
  final A = Matrix<double>(MAXM, MAXM),
      B = Matrix<double>(MAXN, MAXN),
      C = Matrix<double>(MAXM, MAXN),
      CC = Matrix<double>(MAXM, MAXN),
      X = Matrix<double>(MAXM, MAXN),
      SWORK = Matrix<double>(LDSWORK, 126);
  final INFO = Box(0), IINFO = Box(0);

  // Get machine parameters

  final EPS = dlamch('P');
  final SMLNUM = dlamch('S') / EPS;
  final BIGNUM = ONE / SMLNUM;

  const VM = [ONE, 0.000001];

  // Begin test loop

  NINFO[1] = 0;
  NINFO[2] = 0;
  NFAIL[1] = 0;
  NFAIL[2] = 0;
  NFAIL[3] = 0;
  RMAX[1] = ZERO;
  RMAX[2] = ZERO;
  KNT.value = 0;
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = 1;
  }
  const LIWORK = MAXM + MAXN + 2;
  for (final J in [1, 2]) {
    for (final ISGN in [-1, 1]) {
      // Reset seed (overwritten by LATMR)
      for (var I = 1; I <= 4; I++) {
        ISEED[I] = 1;
      }
      for (final M in 32.through(MAXM, step: 71)) {
        final KLA = 0;
        final KUA = M - 1;
        dlatmr(
            M,
            M,
            'S',
            ISEED,
            'N',
            D,
            6,
            ONE,
            ONE,
            'T',
            'N',
            DUML,
            1,
            ONE,
            DUMR,
            1,
            ONE,
            'N',
            IWORK,
            KLA,
            KUA,
            ZERO,
            ONE,
            'NO',
            A,
            MAXM,
            IWORK,
            IINFO);
        for (var I = 1; I <= M; I++) {
          A[I][I] *= VM[J - 1];
        }
        final ANRM = dlange('M', M, M, A, MAXM, DUM);
        for (final N in 51.through(MAXN, step: 47)) {
          final KLB = 0;
          final KUB = N - 1;
          dlatmr(
              N,
              N,
              'S',
              ISEED,
              'N',
              D,
              6,
              ONE,
              ONE,
              'T',
              'N',
              DUML,
              1,
              ONE,
              DUMR,
              1,
              ONE,
              'N',
              IWORK,
              KLB,
              KUB,
              ZERO,
              ONE,
              'NO',
              B,
              MAXN,
              IWORK,
              IINFO);
          final BNRM = dlange('M', N, N, B, MAXN, DUM);
          final TNRM = max(ANRM, BNRM);
          dlatmr(
              M,
              N,
              'S',
              ISEED,
              'N',
              D,
              6,
              ONE,
              ONE,
              'T',
              'N',
              DUML,
              1,
              ONE,
              DUMR,
              1,
              ONE,
              'N',
              IWORK,
              M,
              N,
              ZERO,
              ONE,
              'NO',
              C,
              MAXM,
              IWORK,
              IINFO);
          for (final TRANA in ['N', 'T']) {
            for (final TRANB in ['N', 'T']) {
              final ctx = (
                A: A.copy(),
                B: B.copy(),
                C: C.copy(),
              );
              test('(PARAM=$J ISGN=$ISGN M=$M N=$N TRANA=$TRANA TRANB=$TRANB)',
                  () {
                final (:A, :B, :C) = ctx;
                KNT.value++;

                dlacpy('All', M, N, C, MAXM, X, MAXM);
                dlacpy('All', M, N, C, MAXM, CC, MAXM);
                final SCALE = Box(ONE);
                dtrsyl(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM,
                    SCALE, IINFO);
                if (IINFO.value != 0) NINFO[1]++;
                var XNRM = dlange('M', M, N, X, MAXM, DUM);
                var RMUL = ONE;
                if (XNRM > ONE && TNRM > ONE) {
                  if (XNRM > BIGNUM / TNRM) {
                    RMUL = ONE / max(XNRM, TNRM);
                  }
                }
                dgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM,
                    -SCALE.value * RMUL, CC, MAXM);
                dgemm('N', TRANB, M, N, N, ISGN * RMUL, X, MAXM, B, MAXN, ONE,
                    CC, MAXM);
                var RES1 = dlange('M', M, N, CC, MAXM, DUM);
                var RES = RES1 /
                    max(SMLNUM,
                        max(SMLNUM * XNRM, ((RMUL * TNRM) * EPS) * XNRM));
                if (RES > THRESH) NFAIL[1]++;
                if (RES > RMAX[1]) RMAX[1] = RES;
                test.expect(RES, lessThan(THRESH));

                dlacpy('All', M, N, C, MAXM, X, MAXM);
                dlacpy('All', M, N, C, MAXM, CC, MAXM);
                final SCALE3 = Box(ONE);
                dtrsyl3(TRANA, TRANB, ISGN, M, N, A, MAXM, B, MAXN, X, MAXM,
                    SCALE3, IWORK, LIWORK, SWORK, Box(LDSWORK), INFO);
                if (INFO.value != 0) NINFO[2]++;
                XNRM = dlange('M', M, N, X, MAXM, DUM);
                RMUL = ONE;
                if (XNRM > ONE && TNRM > ONE) {
                  if (XNRM > BIGNUM / TNRM) {
                    RMUL = ONE / max(XNRM, TNRM);
                  }
                }
                dgemm(TRANA, 'N', M, N, M, RMUL, A, MAXM, X, MAXM,
                    -SCALE3.value * RMUL, CC, MAXM);
                dgemm('N', TRANB, M, N, N, ISGN * RMUL, X, MAXM, B, MAXN, ONE,
                    CC, MAXM);
                RES1 = dlange('M', M, N, CC, MAXM, DUM);
                RES = RES1 /
                    max(SMLNUM,
                        max(SMLNUM * XNRM, ((RMUL * TNRM) * EPS) * XNRM));
                // Verify that TRSYL3 only flushes if TRSYL flushes (but
                // there may be cases where TRSYL3 avoid flushing).
                if (SCALE3.value == ZERO && SCALE.value > ZERO ||
                    IINFO.value != INFO.value) {
                  NFAIL[3]++;
                }
                if (RES > THRESH || disnan(RES)) NFAIL[2]++;
                if (RES > RMAX[2]) RMAX[2] = RES;
                test.expect(RES, lessThan(THRESH));
              });
            }
          }
        }
      }
    }
  }
}
