// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlacp2(
  final String UPLO,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  int I, J;

  if (lsame(UPLO, 'U')) {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= min(J, M); I++) {
        B[I][J] = A[I][J].toComplex();
      }
    }
  } else if (lsame(UPLO, 'L')) {
    for (J = 1; J <= N; J++) {
      for (I = J; I <= M; I++) {
        B[I][J] = A[I][J].toComplex();
      }
    }
  } else {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= M; I++) {
        B[I][J] = A[I][J].toComplex();
      }
    }
  }
}
