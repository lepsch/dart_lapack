// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlaipd(
  final int N,
  final Array<Complex> A_,
  final int INDA,
  final int VINDA,
) {
  final A = A_.having();

  final BIGNUM = dlamch('Epsilon') / dlamch('Safe minimum');
  var IA = 1;
  var IXA = INDA;
  for (var I = 1; I <= N; I++) {
    A[IA] = Complex(A[IA].real, BIGNUM);
    IA += IXA;
    IXA += VINDA;
  }
}
