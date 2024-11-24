// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dasum.dart';
import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dbdsqr.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgebd2.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlascl.dart';
import 'package:dart_lapack/src/zlaset.dart';

double zqrt12(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> S_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final S = S_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();

  const ZERO = 0.0, ONE = 1.0;
  final DUMMY = Array<double>(1);
  final INFO = Box(0);

  // Test that enough workspace is supplied

  if (LWORK < M * N + 2 * min(M, N) + max(M, N)) {
    xerbla('ZQRT12', 7);
    return 0;
  }

  // Quick return if possible

  final MN = min(M, N);
  if (MN <= ZERO) return 0;

  final NRMSVL = dnrm2(MN, S, 1);

  // Copy upper triangle of A into work

  zlaset('Full', M, N, Complex.zero, Complex.zero, WORK.asMatrix(), M);
  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= min(J, M); I++) {
      WORK[(J - 1) * M + I] = A[I][J];
    }
  }

  // Get machine parameters

  final SMLNUM = dlamch('S') / dlamch('P');
  final BIGNUM = ONE / SMLNUM;

  // Scale work if max entry outside range [SMLNUM,BIGNUM]

  final ANRM = zlange('M', M, N, WORK.asMatrix(), M, DUMMY);
  final int ISCL;
  if (ANRM > ZERO && ANRM < SMLNUM) {
    // Scale matrix norm up to SMLNUM

    zlascl('G', 0, 0, ANRM, SMLNUM, M, N, WORK.asMatrix(), M, INFO);
    ISCL = 1;
  } else if (ANRM > BIGNUM) {
    // Scale matrix norm down to BIGNUM

    zlascl('G', 0, 0, ANRM, BIGNUM, M, N, WORK.asMatrix(), M, INFO);
    ISCL = 1;
  } else {
    ISCL = 0;
  }

  if (ANRM != ZERO) {
    // Compute SVD of work

    zgebd2(M, N, WORK.asMatrix(), M, RWORK(1), RWORK(MN + 1), WORK(M * N + 1),
        WORK(M * N + MN + 1), WORK(M * N + 2 * MN + 1), INFO);
    dbdsqr('Upper', MN, 0, 0, 0, RWORK(1), RWORK(MN + 1), DUMMY.asMatrix(), MN,
        DUMMY.asMatrix(), 1, DUMMY.asMatrix(), MN, RWORK(2 * MN + 1), INFO);

    if (ISCL == 1) {
      if (ANRM > BIGNUM) {
        dlascl('G', 0, 0, BIGNUM, ANRM, MN, 1, RWORK(1).asMatrix(), MN, INFO);
      }
      if (ANRM < SMLNUM) {
        dlascl('G', 0, 0, SMLNUM, ANRM, MN, 1, RWORK(1).asMatrix(), MN, INFO);
      }
    }
  } else {
    for (var I = 1; I <= MN; I++) {
      RWORK[I] = ZERO;
    }
  }

  // Compare s and singular values of work

  daxpy(MN, -ONE, S, 1, RWORK(1), 1);
  final result = dasum(MN, RWORK(1), 1) / (dlamch('Epsilon') * max(M, N));

  return NRMSVL != ZERO ? result / NRMSVL : result;
}
