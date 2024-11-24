// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

N radix<N extends num>(final N n) {
  return switch (n) {
    double() => 2.0,
    int() => 2,
  } as N;
}

const radix32 = radix;
const radix64 = radix;
