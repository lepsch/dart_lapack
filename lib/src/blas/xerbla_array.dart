// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/xerbla.dart';

void xerbla_array(
  final List<int> SRNAME_ARRAY,
  final int SRNAME_LEN,
  final int INFO,
) {
  var SRNAME = ' ';
  for (var I = 1; I <= min(SRNAME_LEN, SRNAME.length); I++) {
    SRNAME += String.fromCharCode(SRNAME_ARRAY[I]);
  }

  xerbla(SRNAME, INFO);
}
