// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dasum.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlange.dart';
import 'package:dart_lapack/src/dlarf.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlaror.dart';
import 'dlaord.dart';

void dqrt15(
  final int SCALE,
  final int RKSEL,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> S_,
  final Box<int> RANK,
  final Box<double> NORMA,
  final Box<double> NORMB,
  final Array<int> ISEED_,
  final Array<double> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final S = S_.having();
  final ISEED = ISEED_.having(length: 4);
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, SVMIN = 0.1;
  final DUMMY = Array<double>(1);
  final INFO = Box(0);

  final MN = min(M, N);
  if (LWORK < max(M + MN, max(MN * NRHS, 2 * N + M))) {
    xerbla('DQRT15', 16);
    return;
  }

  var SMLNUM = dlamch('Safe minimum');
  var BIGNUM = ONE / SMLNUM;
  final EPS = dlamch('Epsilon');
  SMLNUM = (SMLNUM / EPS) / EPS;
  BIGNUM = ONE / SMLNUM;

  // Determine rank and (unscaled) singular values
  if (RKSEL == 1) {
    RANK.value = MN;
  } else if (RKSEL == 2) {
    RANK.value = (3 * MN) ~/ 4;
    for (var J = RANK.value + 1; J <= MN; J++) {
      S[J] = ZERO;
    }
  } else {
    xerbla('DQRT15', 2);
  }

  if (RANK.value > 0) {
    // Nontrivial case

    S[1] = ONE;
    for (var J = 2; J <= RANK.value; J++) {
      while (true) {
        final TEMP = dlarnd(1, ISEED);
        if (TEMP > SVMIN) {
          S[J] = TEMP.abs();
          break;
        }
      }
    }
    dlaord('Decreasing', RANK.value, S, 1);

    // Generate 'rank' columns of a random orthogonal matrix in A

    dlarnv(2, ISEED, M, WORK);
    dscal(M, ONE / dnrm2(M, WORK, 1), WORK, 1);
    dlaset('Full', M, RANK.value, ZERO, ONE, A, LDA);
    dlarf('Left', M, RANK.value, WORK, 1, TWO, A, LDA, WORK(M + 1));

    // workspace used: m+mn

    // Generate consistent rhs in the range space of A

    dlarnv(2, ISEED, RANK.value * NRHS, WORK);
    dgemm('No transpose', 'No transpose', M, NRHS, RANK.value, ONE, A, LDA,
        WORK.asMatrix(), RANK.value, ZERO, B, LDB);

    // work space used: <= mn *nrhs

    // generate (unscaled) matrix A

    for (var J = 1; J <= RANK.value; J++) {
      dscal(M, S[J], A(1, J).asArray(), 1);
    }
    if (RANK.value < N) {
      dlaset('Full', M, N - RANK.value, ZERO, ZERO, A(1, RANK.value + 1), LDA);
    }
    dlaror('Right', 'No initialization', M, N, A, LDA, ISEED, WORK, INFO);
  } else {
    // work space used 2*n+m

    // Generate null matrix and rhs

    for (var J = 1; J <= MN; J++) {
      S[J] = ZERO;
    }
    dlaset('Full', M, N, ZERO, ZERO, A, LDA);
    dlaset('Full', M, NRHS, ZERO, ZERO, B, LDB);
  }

  // Scale the matrix

  if (SCALE != 1) {
    NORMA.value = dlange('Max', M, N, A, LDA, DUMMY);
    if (NORMA.value != ZERO) {
      if (SCALE == 2) {
        // matrix scaled up

        dlascl('General', 0, 0, NORMA.value, BIGNUM, M, N, A, LDA, INFO);
        dlascl('General', 0, 0, NORMA.value, BIGNUM, MN, 1, S.asMatrix(), MN,
            INFO);
        dlascl('General', 0, 0, NORMA.value, BIGNUM, M, NRHS, B, LDB, INFO);
      } else if (SCALE == 3) {
        // matrix scaled down

        dlascl('General', 0, 0, NORMA.value, SMLNUM, M, N, A, LDA, INFO);
        dlascl('General', 0, 0, NORMA.value, SMLNUM, MN, 1, S.asMatrix(), MN,
            INFO);
        dlascl('General', 0, 0, NORMA.value, SMLNUM, M, NRHS, B, LDB, INFO);
      } else {
        xerbla('DQRT15', 1);
        return;
      }
    }
  }

  NORMA.value = dasum(MN, S, 1);
  NORMB.value = dlange('One-norm', M, NRHS, B, LDB, DUMMY);
}
