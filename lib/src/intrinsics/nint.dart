// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

int nint(final double x) {
  return (x >= 0 ? (x + .5).floor() : -(.5 - x).floor()).toInt();
}
