// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';

void ilaver(final Box<int> VERS_MAJOR, final Box<int> VERS_MINOR,
    final Box<int> VERS_PATCH) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  VERS_MAJOR.value = 3;
  VERS_MINOR.value = 12;
  VERS_PATCH.value = 0;
}
