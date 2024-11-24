// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/lapack.dart';

import '../matgen/dlarnd.dart';
import '../matgen/zlaror.dart';
import 'dlaord.dart';

void zqrt15(
  final int SCALE,
  final int RKSEL,
  final int M,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> S_,
  final Box<int> RANK,
  final Box<double> NORMA,
  final Box<double> NORMB,
  final Array<int> ISEED_,
  final Array<Complex> WORK_,
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
    xerbla('ZQRT15', 16);
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
    xerbla('ZQRT15', 2);
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

    zlarnv(2, ISEED, M, WORK);
    zdscal(M, ONE / dznrm2(M, WORK, 1), WORK, 1);
    zlaset('Full', M, RANK.value, Complex.zero, Complex.one, A, LDA);
    zlarf('Left', M, RANK.value, WORK, 1, TWO.toComplex(), A, LDA, WORK(M + 1));

    // workspace used: m+mn

    // Generate consistent rhs in the range space of A

    zlarnv(2, ISEED, RANK.value * NRHS, WORK);
    zgemm('No transpose', 'No transpose', M, NRHS, RANK.value, Complex.one, A,
        LDA, WORK.asMatrix(), RANK.value, Complex.zero, B, LDB);

    // work space used: <= mn *nrhs

    // generate (unscaled) matrix A

    for (var J = 1; J <= RANK.value; J++) {
      zdscal(M, S[J], A(1, J).asArray(), 1);
    }
    if (RANK.value < N) {
      zlaset('Full', M, N - RANK.value, Complex.zero, Complex.zero,
          A(1, RANK.value + 1), LDA);
    }
    zlaror('Right', 'No initialization', M, N, A, LDA, ISEED, WORK, INFO);
  } else {
    // work space used 2*n+m

    // Generate null matrix and rhs

    for (var J = 1; J <= MN; J++) {
      S[J] = ZERO;
    }
    zlaset('Full', M, N, Complex.zero, Complex.zero, A, LDA);
    zlaset('Full', M, NRHS, Complex.zero, Complex.zero, B, LDB);
  }

  // Scale the matrix

  if (SCALE != 1) {
    NORMA.value = zlange('Max', M, N, A, LDA, DUMMY);
    if (NORMA.value != ZERO) {
      if (SCALE == 2) {
        // matrix scaled up

        zlascl('General', 0, 0, NORMA.value, BIGNUM, M, N, A, LDA, INFO);
        dlascl('General', 0, 0, NORMA.value, BIGNUM, MN, 1, S.asMatrix(), MN,
            INFO);
        zlascl('General', 0, 0, NORMA.value, BIGNUM, M, NRHS, B, LDB, INFO);
      } else if (SCALE == 3) {
        // matrix scaled down

        zlascl('General', 0, 0, NORMA.value, SMLNUM, M, N, A, LDA, INFO);
        dlascl('General', 0, 0, NORMA.value, SMLNUM, MN, 1, S.asMatrix(), MN,
            INFO);
        zlascl('General', 0, 0, NORMA.value, SMLNUM, M, NRHS, B, LDB, INFO);
      } else {
        xerbla('ZQRT15', 1);
        return;
      }
    }
  }

  NORMA.value = dasum(MN, S, 1);
  NORMB.value = zlange('One-norm', M, NRHS, B, LDB, DUMMY);
}
