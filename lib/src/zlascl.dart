// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zlascl(
  final String TYPE,
  final int KL,
  final int KU,
  final double CFROM,
  final double CTO,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  const ZERO = 0.0, ONE = 1.0;
  bool DONE;
  int I, ITYPE, J, K1, K2, K3, K4;
  double BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM;

  // Test the input arguments

  INFO.value = 0;

  if (lsame(TYPE, 'G')) {
    ITYPE = 0;
  } else if (lsame(TYPE, 'L')) {
    ITYPE = 1;
  } else if (lsame(TYPE, 'U')) {
    ITYPE = 2;
  } else if (lsame(TYPE, 'H')) {
    ITYPE = 3;
  } else if (lsame(TYPE, 'B')) {
    ITYPE = 4;
  } else if (lsame(TYPE, 'Q')) {
    ITYPE = 5;
  } else if (lsame(TYPE, 'Z')) {
    ITYPE = 6;
  } else {
    ITYPE = -1;
  }

  if (ITYPE == -1) {
    INFO.value = -1;
  } else if (CFROM == ZERO || disnan(CFROM)) {
    INFO.value = -4;
  } else if (disnan(CTO)) {
    INFO.value = -5;
  } else if (M < 0) {
    INFO.value = -6;
  } else if (N < 0 || (ITYPE == 4 && N != M) || (ITYPE == 5 && N != M)) {
    INFO.value = -7;
  } else if (ITYPE <= 3 && LDA < max(1, M)) {
    INFO.value = -9;
  } else if (ITYPE >= 4) {
    if (KL < 0 || KL > max(M - 1, 0)) {
      INFO.value = -2;
    } else if (KU < 0 ||
        KU > max(N - 1, 0) ||
        ((ITYPE == 4 || ITYPE == 5) && KL != KU)) {
      INFO.value = -3;
    } else if ((ITYPE == 4 && LDA < KL + 1) ||
        (ITYPE == 5 && LDA < KU + 1) ||
        (ITYPE == 6 && LDA < 2 * KL + KU + 1)) {
      INFO.value = -9;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZLASCL', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || M == 0) return;

  // Get machine parameters

  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;

  CFROMC = CFROM;
  CTOC = CTO;

  do {
    CFROM1 = CFROMC * SMLNUM;
    if (CFROM1 == CFROMC) {
      // CFROMC is an inf.  Multiply by a correctly signed zero for
      // finite CTOC, or a NaN if CTOC is infinite.
      MUL = CTOC / CFROMC;
      DONE = true;
      CTO1 = CTOC;
    } else {
      CTO1 = CTOC / BIGNUM;
      if (CTO1 == CTOC) {
        // CTOC is either 0 or an inf.  In both cases, CTOC itself
        // serves as the correct multiplication factor.
        MUL = CTOC;
        DONE = true;
        CFROMC = ONE;
      } else if (CFROM1.abs() > CTOC.abs() && CTOC != ZERO) {
        MUL = SMLNUM;
        DONE = false;
        CFROMC = CFROM1;
      } else if (CTO1.abs() > CFROMC.abs()) {
        MUL = BIGNUM;
        DONE = false;
        CTOC = CTO1;
      } else {
        MUL = CTOC / CFROMC;
        DONE = true;
        if (MUL == ONE) return;
      }
    }

    if (ITYPE == 0) {
      // Full matrix

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          A[I][J] *= MUL.toComplex();
        }
      }
    } else if (ITYPE == 1) {
      // Lower triangular matrix

      for (J = 1; J <= N; J++) {
        for (I = J; I <= M; I++) {
          A[I][J] *= MUL.toComplex();
        }
      }
    } else if (ITYPE == 2) {
      // Upper triangular matrix

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= min(J, M); I++) {
          A[I][J] *= MUL.toComplex();
        }
      }
    } else if (ITYPE == 3) {
      // Upper Hessenberg matrix

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= min(J + 1, M); I++) {
          A[I][J] *= MUL.toComplex();
        }
      }
    } else if (ITYPE == 4) {
      // Lower half of a symmetric band matrix

      K3 = KL + 1;
      K4 = N + 1;
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= min(K3, K4 - J); I++) {
          A[I][J] *= MUL.toComplex();
        }
      }
    } else if (ITYPE == 5) {
      // Upper half of a symmetric band matrix

      K1 = KU + 2;
      K3 = KU + 1;
      for (J = 1; J <= N; J++) {
        for (I = max(K1 - J, 1); I <= K3; I++) {
          A[I][J] *= MUL.toComplex();
        }
      }
    } else if (ITYPE == 6) {
      // Band matrix

      K1 = KL + KU + 2;
      K2 = KL + 1;
      K3 = 2 * KL + KU + 1;
      K4 = KL + KU + 1 + M;
      for (J = 1; J <= N; J++) {
        for (I = max(K1 - J, K2); I <= min(K3, K4 - J); I++) {
          A[I][J] *= MUL.toComplex();
        }
      }
    }
  } while (!DONE);
}
