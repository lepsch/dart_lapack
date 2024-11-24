// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:typed_data';

const N Function<N extends num>(N n) huge = huge64;

N huge64<N extends num>(final N n) {
  return switch (n) {
    double() => double.maxFinite,
    int() => 0x7fffffffffffffff,
  } as N;
}

N huge32<N extends num>(final N n) {
  return switch (n) {
    double() => Uint32List.fromList([0x7F7FFFFF]).buffer.asFloat32List().first,
    int() => 0x7fffffff,
  } as N;
}
