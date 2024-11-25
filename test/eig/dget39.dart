// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dasum.dart';
import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/ddot.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlaqtr.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:test/test.dart';

import '../test_driver.dart';

void dget39(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
  final TestDriver test,
  final double THRESH,
) {
  const LDT = 10, LDT2 = 2 * LDT;
  const ZERO = 0.0, ONE = 1.0;
  final T = Matrix<double>(LDT, LDT);
  final B = Array<double>(LDT),
      D = Array<double>(LDT2),
      WORK = Array<double>(LDT),
      X = Array<double>(LDT2),
      Y = Array<double>(LDT2);
  const IDIM = [4, 5, 5, 5, 5, 5];
  final IVAL = Matrix3d.fromData([
    3, 0, 0, 0, 0, 1, 1, -1, 0, 0, 3, 2, 1, 0, 0, 4, 3, 2, 2, //
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 3, 3, 4, //
    0, 0, 4, 2, 2, 3, 0, 1, 1, 1, 1, 5, 1, 0, 0, 0, 0, 2, 4, //
    -2, 0, 0, 3, 3, 4, 0, 0, 4, 2, 2, 3, 0, 1, 1, 1, 1, 1, 1, //
    0, 0, 0, 0, 2, 1, -1, 0, 0, 9, 8, 1, 0, 0, 4, 9, 1, 2, -1, //
    2, 2, 2, 2, 2, 9, 0, 0, 0, 0, 6, 4, 0, 0, 0, 3, 2, 1, 1, //
    0, 5, 1, -1, 1, 0, 2, 2, 2, 2, 2, 4, 0, 0, 0, 0, 2, 2, 0, //
    0, 0, 1, 4, 4, 0, 0, 2, 4, 2, 2, -1, 2, 2, 2, 2, 2,
  ], (
    5, 5, 6 //
  ));

  // Get machine parameters

  final EPS = dlamch('P');
  final SMLNUM = dlamch('S') / EPS;
  final BIGNUM = ONE / dlamch('S');

  // Set up test case parameters

  final sqrtSmallNum = sqrt(SMLNUM);
  final sqrtBigNum = sqrt(BIGNUM);
  final VM1 = [ONE, sqrtSmallNum, sqrtSmallNum, sqrtBigNum, sqrtBigNum];
  final VM2 = [ONE, sqrtSmallNum, sqrtSmallNum, sqrtBigNum, sqrtBigNum];
  final VM3 = [ONE, sqrtSmallNum, sqrtSmallNum, sqrtBigNum, sqrtBigNum];
  final VM4 = [ONE, sqrtSmallNum, sqrtSmallNum, sqrtBigNum, sqrtBigNum];
  final VM5 = [ONE, EPS, sqrtSmallNum];

  // Initialization

  KNT.value = 0;
  RMAX.value = ZERO;
  NINFO.value = 0;

  // Begin test loop

  for (final IVM5 in 1.through(3)) {
    for (final IVM4 in 1.through(5)) {
      for (final IVM3 in 1.through(5)) {
        for (final IVM2 in 1.through(5)) {
          for (final IVM1 in 1.through(5)) {
            for (final NDIM in 1.through(6)) {
              test(
                  'DGET39 IVM5=$IVM5 IVM4=$IVM4 IVM3=$IVM3 IVM2=$IVM2 IVM1=$IVM1 NDIM=$NDIM',
                  () {
                final INFO = Box(0);

                final N = IDIM[NDIM - 1];
                for (var I = 1; I <= N; I++) {
                  for (var J = 1; J <= N; J++) {
                    T[I][J] = IVAL[I][J][NDIM] * VM1[IVM1 - 1];
                    if (I >= J) T[I][J] *= VM5[IVM5 - 1];
                  }
                }

                final W = ONE * VM2[IVM2 - 1];

                for (var I = 1; I <= N; I++) {
                  B[I] = cos(I) * VM3[IVM3 - 1];
                }

                for (var I = 1; I <= 2 * N; I++) {
                  D[I] = sin(I) * VM4[IVM4 - 1];
                }

                final NORM = dlange('1', N, N, T, LDT, WORK);
                final K = idamax(N, B, 1);
                final NORMTB = NORM + B[K].abs() + W.abs();

                dcopy(N, D, 1, X, 1);
                KNT.value++;
                final DUM = Array<double>(1);
                final SCALE = Box(0.0);
                dlaqtr(false, true, N, T, LDT, DUM, 0, SCALE, X, WORK, INFO);
                if (INFO.value != 0) NINFO.value++;

                // || T*x - scale*d || /
                //    max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)

                dcopy(N, D, 1, Y, 1);
                dgemv('No transpose', N, N, ONE, T, LDT, X, 1, -SCALE.value, Y,
                    1);
                var XNORM = dasum(N, X, 1);
                var RESID = dasum(N, Y, 1);
                var DOMIN = max(
                    SMLNUM, max((SMLNUM / EPS) * NORM, (NORM * EPS) * XNORM));
                RESID /= DOMIN;
                test.expect(RESID, lessThan(THRESH));
                if (RESID > RMAX.value) {
                  RMAX.value = RESID;
                  LMAX.value = KNT.value;
                }

                dcopy(N, D, 1, X, 1);
                KNT.value++;
                dlaqtr(true, true, N, T, LDT, DUM, 0, SCALE, X, WORK, INFO);
                if (INFO.value != 0) NINFO.value++;

                // || T*x - scale*d || /
                //    max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)

                dcopy(N, D, 1, Y, 1);
                dgemv('Transpose', N, N, ONE, T, LDT, X, 1, -SCALE.value, Y, 1);
                XNORM = dasum(N, X, 1);
                RESID = dasum(N, Y, 1);
                DOMIN = max(
                    SMLNUM, max((SMLNUM / EPS) * NORM, (NORM * EPS) * XNORM));
                RESID /= DOMIN;
                test.expect(RESID, lessThan(THRESH));
                if (RESID > RMAX.value) {
                  RMAX.value = RESID;
                  LMAX.value = KNT.value;
                }

                dcopy(2 * N, D, 1, X, 1);
                KNT.value++;
                dlaqtr(false, false, N, T, LDT, B, W, SCALE, X, WORK, INFO);
                if (INFO.value != 0) NINFO.value++;

                // ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
                //    max(ulp*(||T||+||B||)*(||x1||+||x2||),
                //             smlnum/ulp * (||T||+||B||), smlnum )

                dcopy(2 * N, D, 1, Y, 1);
                Y[1] = ddot(N, B, 1, X(1 + N), 1) + SCALE.value * Y[1];
                for (var I = 2; I <= N; I++) {
                  Y[I] = W * X[I + N] + SCALE.value * Y[I];
                }
                dgemv('No transpose', N, N, ONE, T, LDT, X, 1, -ONE, Y, 1);

                Y[1 + N] = ddot(N, B, 1, X, 1) - SCALE.value * Y[1 + N];
                for (var I = 2; I <= N; I++) {
                  Y[I + N] = W * X[I] - SCALE.value * Y[I + N];
                }
                dgemv('No transpose', N, N, ONE, T, LDT, X(1 + N), 1, ONE,
                    Y(1 + N), 1);

                RESID = dasum(2 * N, Y, 1);
                DOMIN = max(
                    SMLNUM,
                    max(
                      (SMLNUM / EPS) * NORMTB,
                      EPS * (NORMTB * dasum(2 * N, X, 1)),
                    ));
                RESID /= DOMIN;
                if (RESID > RMAX.value) {
                  RMAX.value = RESID;
                  LMAX.value = KNT.value;
                }

                dcopy(2 * N, D, 1, X, 1);
                KNT.value++;
                dlaqtr(true, false, N, T, LDT, B, W, SCALE, X, WORK, INFO);
                if (INFO.value != 0) NINFO.value++;

                // ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
                //    max(ulp*(||T||+||B||)*(||x1||+||x2||),
                //             smlnum/ulp * (||T||+||B||), smlnum )

                dcopy(2 * N, D, 1, Y, 1);
                Y[1] = B[1] * X[1 + N] - SCALE.value * Y[1];
                for (var I = 2; I <= N; I++) {
                  Y[I] = B[I] * X[1 + N] + W * X[I + N] - SCALE.value * Y[I];
                }
                dgemv('Transpose', N, N, ONE, T, LDT, X, 1, ONE, Y, 1);

                Y[1 + N] = B[1] * X[1] + SCALE.value * Y[1 + N];
                for (var I = 2; I <= N; I++) {
                  Y[I + N] = B[I] * X[1] + W * X[I] + SCALE.value * Y[I + N];
                }
                dgemv('Transpose', N, N, ONE, T, LDT, X(1 + N), 1, -ONE,
                    Y(1 + N), 1);

                RESID = dasum(2 * N, Y, 1);
                DOMIN = max(
                    SMLNUM,
                    max(
                      (SMLNUM / EPS) * NORMTB,
                      EPS * (NORMTB * dasum(2 * N, X, 1)),
                    ));
                RESID /= DOMIN;
                test.expect(RESID, lessThan(THRESH));
                if (RESID > RMAX.value) {
                  RMAX.value = RESID;
                  LMAX.value = KNT.value;
                }
              });
            }
          }
        }
      }
    }
  }
}
