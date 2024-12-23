// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:typed_data';

const double Function(double _) tiny = tiny64;

double tiny64(final double _) {
  // minimum **normal** positive 64 bits float (exp=1, m=0)
  return Uint32List.fromList([0x00000000, 0x00100000])
      .buffer
      .asFloat64List(0)
      .first;
}

double tiny32(final double _) {
  // minimum **normal** positive 32 bits float (exp=1, m=0)
  return Uint32List.fromList([0x00800000]).buffer.asFloat32List(0).first;
}
