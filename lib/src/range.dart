// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

extension RangeExtension on int {
  /// Returns a generator from [this] up to, but not including, the [end].
  Iterable<int> to(int end, {int step = 1}) sync* {
    for (var i = this; step < 0 ? i > end : i < end; i += step) {
      yield i;
    }
  }

  /// Returns a generator from [this] through the [end].
  Iterable<int> through(int end, {int step = 1}) sync* {
    for (var i = this; step < 0 ? i >= end : i <= end; i += step) {
      yield i;
    }
  }
}
