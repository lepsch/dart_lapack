// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

T sign<T extends num>(final T a, final num b) {
  return a.abs() * (b.isNegative ? -1 : 1) as T;
}
