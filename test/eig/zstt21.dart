// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zher.dart';
import 'package:lapack/src/blas/zher2.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlaset.dart';

void zstt21(
  final int N,
  final int KBAND,
  final Array<double> AD_,
  final Array<double> AE_,
  final Array<double> SD_,
  final Array<double> SE_,
  final Matrix<Complex> U_,
  final int LDU,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final AD = AD_.having();
  final AE = AE_.having();
  final SD = SD_.having();
  final SE = SE_.having();
  final RESULT = RESULT_.having();
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int J;
  double ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM;

  // 1)      Constants

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');

  // Do Test 1

  // Copy A & Compute its 1-Norm:

  zlaset('Full', N, N, Complex.zero, Complex.zero, WORK.asMatrix(), N);

  ANORM = ZERO;
  TEMP1 = ZERO;

  for (J = 1; J <= N - 1; J++) {
    WORK[(N + 1) * (J - 1) + 1] = AD[J].toComplex();
    WORK[(N + 1) * (J - 1) + 2] = AE[J].toComplex();
    TEMP2 = AE[J].abs();
    ANORM = max(ANORM, AD[J].abs() + TEMP1 + TEMP2);
    TEMP1 = TEMP2;
  }

  WORK[pow(N, 2).toInt()] = AD[N].toComplex();
  ANORM = max(ANORM, max(AD[N].abs() + TEMP1, UNFL));

  // Norm of A - USU*

  for (J = 1; J <= N; J++) {
    zher('L', N, -SD[J], U(1, J).asArray(), 1, WORK.asMatrix(), N);
  }

  if (N > 1 && KBAND == 1) {
    for (J = 1; J <= N - 1; J++) {
      zher2('L', N, -SE[J].toComplex(), U(1, J).asArray(), 1,
          U(1, J + 1).asArray(), 1, WORK.asMatrix(), N);
    }
  }

  WNORM = zlanhe('1', 'L', N, WORK.asMatrix(), N, RWORK);

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

  RESULT[2] = min(N, zlange('1', N, N, WORK.asMatrix(), N, RWORK)) / (N * ULP);
}
