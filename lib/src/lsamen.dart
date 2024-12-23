// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';

bool lsamen(final int N, final String CA, final String CB) {
  if (CA.length < N || CB.length < N) return false;

  // Do for each character in the two strings.

  for (var I = 1; I <= N; I++) {
    // Test if the characters are equal using lsame.
    if (!lsame(CA[I - 1], CB[I - 1])) return false;
  }
  return true;
}
