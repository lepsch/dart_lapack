// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlassq.dart';

double zlanhe(
  final String NORM,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> WORK_,
) {
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  double ABSA, VALUE = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N == 0) {
    VALUE = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J - 1; I++) {
          SUM.value = A[I][J].abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        }
        SUM.value = A[J][J].real.abs();
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    } else {
      for (J = 1; J <= N; J++) {
        SUM.value = A[J][J].real.abs();
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        for (I = J + 1; I <= N; I++) {
          SUM.value = A[I][J].abs();
          if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
        }
      }
    }
  } else if ((lsame(NORM, 'I')) || (lsame(NORM, 'O')) || (NORM == '1')) {
    // Find normI(A) ( = norm1(A), since A is hermitian).

    VALUE = ZERO;
    if (lsame(UPLO, 'U')) {
      for (J = 1; J <= N; J++) {
        SUM.value = ZERO;
        for (I = 1; I <= J - 1; I++) {
          ABSA = A[I][J].abs();
          SUM.value += ABSA;
          WORK[I] += ABSA;
        }
        WORK[J] = SUM.value + A[J][J].real.abs();
      }
      for (I = 1; I <= N; I++) {
        SUM.value = WORK[I];
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    } else {
      for (I = 1; I <= N; I++) {
        WORK[I] = ZERO;
      }
      for (J = 1; J <= N; J++) {
        SUM.value = WORK[J] + A[J][J].real.abs();
        for (I = J + 1; I <= N; I++) {
          ABSA = A[I][J].abs();
          SUM.value += ABSA;
          WORK[I] += ABSA;
        }
        if (VALUE < SUM.value || disnan(SUM.value)) VALUE = SUM.value;
      }
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    if (lsame(UPLO, 'U')) {
      for (J = 2; J <= N; J++) {
        zlassq(J - 1, A(1, J).asArray(), 1, SCALE, SUM);
      }
    } else {
      for (J = 1; J <= N - 1; J++) {
        zlassq(N - J, A(J + 1, J).asArray(), 1, SCALE, SUM);
      }
    }
    SUM.value = 2 * SUM.value;
    for (I = 1; I <= N; I++) {
      if (A[I][I].real != ZERO) {
        ABSA = A[I][I].real.abs();
        if (SCALE.value < ABSA) {
          SUM.value = ONE + SUM.value * pow(SCALE.value / ABSA, 2);
          SCALE.value = ABSA;
        } else {
          SUM.value += pow(ABSA / SCALE.value, 2);
        }
      }
    }
    VALUE = SCALE.value * sqrt(SUM.value);
  }

  return VALUE;
}
