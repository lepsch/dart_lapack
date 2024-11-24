// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

final _timer = Stopwatch()..start();

double dsecnd() {
  // -- LAPACK auxiliary routine --
  // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  // =====================================================================

  return _timer.elapsed.inMilliseconds / 1000;
}
