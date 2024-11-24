// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zhpr.dart';
import 'package:dart_lapack/src/blas/zhpr2.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlanhb.dart';
import 'package:dart_lapack/src/zlanhp.dart';

void zhbt21(
  final String UPLO,
  final int N,
  final int KA,
  final int KS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> U_,
  final int LDU,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final D = D_.having();
  final E = E_.having();
  final RESULT = RESULT_.having(length: 2);
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER;
  String CUPLO;
  int IKA, J, JC, JR;
  double ANORM, ULP, UNFL, WNORM;

  // Constants

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  IKA = max(0, min(N - 1, KA));

  if (lsame(UPLO, 'U')) {
    LOWER = false;
    CUPLO = 'U';
  } else {
    LOWER = true;
    CUPLO = 'L';
  }

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');

  // Some Error Checks

  // Do Test 1

  // Norm of A:

  ANORM = max(zlanhb('1', CUPLO, N, IKA, A, LDA, RWORK), UNFL);

  // Compute error matrix:    Error = A - U S U**H

  // Copy A from SB to SP storage format.

  J = 0;
  for (JC = 1; JC <= N; JC++) {
    if (LOWER) {
      for (JR = 1; JR <= min(IKA + 1, N + 1 - JC); JR++) {
        J++;
        WORK[J] = A[JR][JC];
      }
      for (JR = IKA + 2; JR <= N + 1 - JC; JR++) {
        J++;
        WORK[J] = Complex.zero;
      }
    } else {
      for (JR = IKA + 2; JR <= JC; JR++) {
        J++;
        WORK[J] = Complex.zero;
      }
      for (JR = min(IKA, JC - 1); JR >= 0; JR--) {
        J++;
        WORK[J] = A[IKA + 1 - JR][JC];
      }
    }
  }

  for (J = 1; J <= N; J++) {
    zhpr(CUPLO, N, -D[J], U(1, J).asArray(), 1, WORK);
  }

  if (N > 1 && KS == 1) {
    for (J = 1; J <= N - 1; J++) {
      zhpr2(CUPLO, N, -E[J].toComplex(), U(1, J).asArray(), 1,
          U(1, J + 1).asArray(), 1, WORK);
    }
  }
  WNORM = zlanhp('1', CUPLO, N, WORK, RWORK);

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (N * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, N) / (N * ULP);
    }
  }

  // Do Test 2

  // Compute  U U**H - I

  zgemm('N', 'C', N, N, N, Complex.one, U, LDU, U, LDU, Complex.zero,
      WORK.asMatrix(), N);

  for (J = 1; J <= N; J++) {
    WORK[(N + 1) * (J - 1) + 1] -= Complex.one;
  }

  RESULT[2] = min(zlange('1', N, N, WORK.asMatrix(), N, RWORK), N) / (N * ULP);
}
