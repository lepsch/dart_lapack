// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dtrsyl.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/range.dart';
import 'package:test/test.dart';

import '../test_driver.dart';

void dget35(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
  final TestDriver test,
  final double THRESH,
) {
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0, FOUR = 4.0;
  final INFO = Box(0);
  final SCALE = Box(0.0);
  final A = Matrix<double>(6, 6),
      B = Matrix<double>(6, 6),
      C = Matrix<double>(6, 6),
      CC = Matrix<double>(6, 6);
  const IDIM = [1, 2, 3, 4, 3, 3, 6, 4];
  final IVAL = Matrix3d<double>.fromData(<double>[
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, //
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    1, 0, 0, 0, 0, 0, 5, 1, 2, 0, 0, 0, -8, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 4, 0, 0, 0, 0, -5, 3, 0, 0, 0, 0, //
    1, 2, 1, 4, 0, 0, -3, -9, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, //
    1, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0, 0, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 3, -4, 0, 0, 0, //
    2, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    1, 2, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 5, 6, 3, 4, 0, 0, -1, -9, -5, 2, 0,
    0, //
    8, 8, 8, 8, 5, 6, 9, 9, 9, 9, -7, 5, 1, 0, 0, 0, 0, 0, 1, 5, 2, 0, 0, 0, //
    2, -21, 5, 0, 0, 0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
  ], (
    6, 6, 8 //
  ));

  // Get machine parameters

  final EPS = dlamch('P');
  final SMLNUM = dlamch('S') * FOUR / EPS;
  final BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  final VM1 = [sqrt(SMLNUM), ONE, sqrt(BIGNUM)];
  final VM2 = [ONE, ONE + TWO * EPS, TWO];

  KNT.value = 0;
  NINFO.value = 0;
  LMAX.value = 0;
  RMAX.value = ZERO;

  // Begin test loop

  for (final TRANA in ['N', 'T']) {
    for (final TRANB in ['N', 'T']) {
      for (final ISGN in [-1, 1]) {
        for (final IMA in 1.through(8)) {
          for (final IMLDA1 in 1.through(3)) {
            for (final IMLDA2 in 1.through(3)) {
              for (final IMLOFF in 1.through(2)) {
                for (final IMB in 1.through(8)) {
                  for (final IMLDB1 in 1.through(3)) {
                    test(
                        'DGET35 TRANA=$TRANA TRANB=$TRANB ISGN=$ISGN IMA=$IMA IMLDA1=$IMLDA1 IMLDA2=$IMLDA2 IMLOFF=$IMLOFF IMB=$IMB IMLDB1=$IMLDB1',
                        () {
                      final M = IDIM[IMA - 1];
                      final N = IDIM[IMB - 1];
                      var TNRM = ZERO;
                      for (var I = 1; I <= M; I++) {
                        for (var J = 1; J <= M; J++) {
                          A[I][J] = IVAL[I][J][IMA];
                          if ((I - J).abs() <= 1) {
                            A[I][J] *= VM1[IMLDA1 - 1];
                            A[I][J] *= VM2[IMLDA2 - 1];
                          } else {
                            A[I][J] *= VM1[IMLOFF - 1];
                          }
                          TNRM = max(TNRM, A[I][J].abs());
                        }
                      }
                      for (var I = 1; I <= N; I++) {
                        for (var J = 1; J <= N; J++) {
                          B[I][J] = IVAL[I][J][IMB];
                          if ((I - J).abs() <= 1) {
                            B[I][J] *= VM1[IMLDB1 - 1];
                          } else {
                            B[I][J] *= VM1[IMLOFF - 1];
                          }
                          TNRM = max(TNRM, B[I][J].abs());
                        }
                      }
                      var CNRM = ZERO;
                      for (var I = 1; I <= M; I++) {
                        for (var J = 1; J <= N; J++) {
                          C[I][J] = sin(I * J);
                          CNRM = max(CNRM, C[I][J]);
                          CC[I][J] = C[I][J];
                        }
                      }
                      KNT.value++;
                      dtrsyl(TRANA, TRANB, ISGN, M, N, A, 6, B, 6, C, 6, SCALE,
                          INFO);
                      if (INFO.value != 0) NINFO.value++;
                      final DUM = Array<double>(1);
                      final XNRM = dlange('M', M, N, C, 6, DUM);
                      var RMUL = ONE;
                      if (XNRM > ONE && TNRM > ONE) {
                        if (XNRM > BIGNUM / TNRM) {
                          RMUL = ONE / max(XNRM, TNRM);
                        }
                      }
                      dgemm(TRANA, 'N', M, N, M, RMUL, A, 6, C, 6,
                          -SCALE.value * RMUL, CC, 6);
                      dgemm('N', TRANB, M, N, N, ISGN * RMUL, C, 6, B, 6, ONE,
                          CC, 6);

                      final RES1 = dlange('M', M, N, CC, 6, DUM);
                      final RES = RES1 /
                          max(SMLNUM,
                              max(SMLNUM * XNRM, ((RMUL * TNRM) * EPS) * XNRM));
                      if (RES > RMAX.value) {
                        LMAX.value = KNT.value;
                        RMAX.value = RES;
                      }

                      test.expect(RES, lessThan(THRESH));
                    });
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
