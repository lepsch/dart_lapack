// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlaord(
  final String JOB,
  final int N,
  final Array<double> X_,
  final int INCX,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();

  final INC = INCX.abs();
  if (lsame(JOB, 'I')) {
    // Sort in increasing order

    sortLoop:
    for (var I = 2; I <= N; I++) {
      var IX = 1 + (I - 1) * INC;
      while (true) {
        if (IX == 1) continue sortLoop;
        final IXNEXT = IX - INC;
        if (X[IX] > X[IXNEXT]) {
          continue sortLoop;
        } else {
          final TEMP = X[IX];
          X[IX] = X[IXNEXT];
          X[IXNEXT] = TEMP;
        }
        IX = IXNEXT;
      }
    }
  } else if (lsame(JOB, 'D')) {
    // Sort in decreasing order

    sortLoop:
    for (var I = 2; I <= N; I++) {
      var IX = 1 + (I - 1) * INC;
      while (true) {
        if (IX == 1) continue sortLoop;
        final IXNEXT = IX - INC;
        if (X[IX] < X[IXNEXT]) {
          continue sortLoop;
        } else {
          final TEMP = X[IX];
          X[IX] = X[IXNEXT];
          X[IXNEXT] = TEMP;
        }
        IX = IXNEXT;
      }
    }
  }
}
