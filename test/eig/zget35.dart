// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/ztrsyl.dart';

Future<void> zget35(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
  final Nin NIN,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const LDT = 10;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  const LARGE = 1.0e6;
  String TRANA = '', TRANB = '';
  int I, IMLA, IMLAD, IMLB, IMLC, ISGN, ITRANA, ITRANB, J, M, N;
  double BIGNUM, EPS, RES, RES1, SMLNUM, TNRM, XNRM;
  Complex RMUL;
  final DUM = Array<double>(1), VM1 = Array<double>(3), VM2 = Array<double>(3);
  final A = Matrix<Complex>(LDT, LDT),
      ATMP = Matrix<Complex>(LDT, LDT),
      B = Matrix<Complex>(LDT, LDT),
      BTMP = Matrix<Complex>(LDT, LDT),
      C = Matrix<Complex>(LDT, LDT),
      CSAV = Matrix<Complex>(LDT, LDT),
      CTMP = Matrix<Complex>(LDT, LDT);
  final INFO = Box(0);
  final SCALE = Box(0.0);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  VM1[1] = sqrt(SMLNUM);
  VM1[2] = ONE;
  VM1[3] = LARGE;
  VM2[1] = ONE;
  VM2[2] = ONE + TWO * EPS;
  VM2[3] = TWO;

  KNT.value = 0;
  NINFO.value = 0;
  LMAX.value = 0;
  RMAX.value = ZERO;

  // Begin test loop
  while (true) {
    (M, N) = await NIN.readInt2();
    if (N == 0) return;
    await NIN.readMatrix(ATMP, M, M);
    await NIN.readMatrix(BTMP, N, N);
    await NIN.readMatrix(CTMP, M, N);
    for (IMLA = 1; IMLA <= 3; IMLA++) {
      for (IMLAD = 1; IMLAD <= 3; IMLAD++) {
        for (IMLB = 1; IMLB <= 3; IMLB++) {
          for (IMLC = 1; IMLC <= 3; IMLC++) {
            for (ITRANA = 1; ITRANA <= 2; ITRANA++) {
              for (ITRANB = 1; ITRANB <= 2; ITRANB++) {
                for (ISGN = -1; ISGN <= 1; ISGN += 2) {
                  if (ITRANA == 1) TRANA = 'N';
                  if (ITRANA == 2) TRANA = 'C';
                  if (ITRANB == 1) TRANB = 'N';
                  if (ITRANB == 2) TRANB = 'C';
                  TNRM = ZERO;
                  for (I = 1; I <= M; I++) {
                    for (J = 1; J <= M; J++) {
                      A[I][J] = ATMP[I][J] * VM1[IMLA].toComplex();
                      TNRM = max(TNRM, A[I][J].abs());
                    }
                    A[I][I] *= VM2[IMLAD].toComplex();
                    TNRM = max(TNRM, A[I][I].abs());
                  }
                  for (I = 1; I <= N; I++) {
                    for (J = 1; J <= N; J++) {
                      B[I][J] = BTMP[I][J] * VM1[IMLB].toComplex();
                      TNRM = max(TNRM, B[I][J].abs());
                    }
                  }
                  if (TNRM == ZERO) TNRM = ONE;
                  for (I = 1; I <= M; I++) {
                    for (J = 1; J <= N; J++) {
                      C[I][J] = CTMP[I][J] * VM1[IMLC].toComplex();
                      CSAV[I][J] = C[I][J];
                    }
                  }
                  KNT.value++;
                  ztrsyl(TRANA, TRANB, ISGN, M, N, A, LDT, B, LDT, C, LDT,
                      SCALE, INFO);
                  if (INFO.value != 0) NINFO.value++;
                  XNRM = zlange('M', M, N, C, LDT, DUM);
                  RMUL = Complex.one;
                  if (XNRM > ONE && TNRM > ONE) {
                    if (XNRM > BIGNUM / TNRM) {
                      RMUL = max(XNRM, TNRM).toComplex();
                      RMUL = Complex.one / RMUL;
                    }
                  }
                  zgemm(TRANA, 'N', M, N, M, RMUL, A, LDT, C, LDT,
                      -SCALE.value.toComplex() * RMUL, CSAV, LDT);
                  zgemm('N', TRANB, M, N, N, ISGN.toComplex() * RMUL, C, LDT, B,
                      LDT, Complex.one, CSAV, LDT);
                  RES1 = zlange('M', M, N, CSAV, LDT, DUM);
                  RES = RES1 /
                      max(
                          SMLNUM,
                          max(SMLNUM * XNRM,
                              ((RMUL.abs() * TNRM) * EPS) * XNRM));
                  if (RES > RMAX.value) {
                    LMAX.value = KNT.value;
                    RMAX.value = RES;
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
