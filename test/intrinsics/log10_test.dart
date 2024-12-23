// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/intrinsics/log10.dart';
import 'package:test/test.dart';

void main() {
  group('Log10', () {
    test('Known', () {
      expect(log10(0), double.negativeInfinity);
      expect(log10(-1), isNaN);
      expect(log10(1000000), 6);
      expect(log(1000000) / ln10, isNot(6));
      expect(log(1000000) / ln10, closeTo(6, 0.0000000001));
    });
  });
}
