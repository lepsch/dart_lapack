// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

void Function(String SRNAME, int INFO) xerbla = _xerbla;

void _xerbla(final String SRNAME, final int INFO) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  throw ArgumentError(
      ' ** On entry to $SRNAME parameter number $INFO had an illegal value');
}
