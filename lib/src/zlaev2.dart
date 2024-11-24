// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlaev2.dart';

void zlaev2(
  final Complex A,
  final Complex B,
  final Complex C,
  final Box<double> RT1,
  final Box<double> RT2,
  final Box<double> CS1,
  final Box<Complex> SN1,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  final T = Box(0.0);
  Complex W;

  if (B.abs() == ZERO) {
    W = Complex.one;
  } else {
    W = B.conjugate() / B.abs().toComplex();
  }
  dlaev2(A.real, B.abs(), C.real, RT1, RT2, CS1, T);
  SN1.value = W * T.value.toComplex();
}
