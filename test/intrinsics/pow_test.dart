// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/intrinsics/pow.dart';
import 'package:test/test.dart';

void main() {
  group('Power **', () {
    for (final (exponent, expected, delta) in [
      (1, 0.1500128399472628881, 0.1e-16),
      (2, 0.2250385214904311245e-1, 0.1e-17),
      (3, 0.3375866770631272273e-2, 0.1e-18),
      (4, 0.5064233615459923492e-3, 0.1e-19),
      (5, 0.7597000668115379176e-4, 0.1e-20),
      (6, 0.1139647645305241703e-4, 0.1e-20),
      (7, 0.1709617798114502428e-5, 0.1e-21),
      (8, 0.2564646211195428605e-6, 0.1e-22),
      (9, 0.3847298616014139761e-7, 0.1e-23),
      (10, 0.5771441915134552535e-8, 0.1e-24),
      (11, 0.8657903922800039815e-9, 0.1e-25),
      (12, 0.1298796755449781962e-9, 0.1e-25),
      (13, 0.1948361897993124646e-10, 0.1e-26),
      (14, 0.2922793015629879813e-11, 0.1e-27),
      (15, 0.4384564808526629749e-12, 0.1e-28),
      (16, 0.6577410188599066658e-13, 0.1e-29),
    ]) {
      test('int exponent', () async {
        const alpha = 0.1500128399472628881;
        expect(
          pow(alpha, exponent),
          closeTo(expected, delta),
          reason: 'alpha^$exponent',
        );
      });
    }
  });
}
