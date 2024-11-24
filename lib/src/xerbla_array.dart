// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/xerbla.dart';

void xerbla_array(
  final String SRNAME_ARRAY,
  final int SRNAME_LEN,
  final int INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I;
  String SRNAME;

  SRNAME = ' ';
  for (I = 1; I <= min(SRNAME_LEN, SRNAME.length); I++) {
    SRNAME += SRNAME_ARRAY[I - 1];
  }

  xerbla(SRNAME, INFO);
}
