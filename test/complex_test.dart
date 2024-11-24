// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:test/test.dart';

void main() {
  group('Complex', () {
    test('Addition', () {
      expect(Complex(0, 0) + Complex.zero, Complex.zero);
      expect(Complex(1.1, 2.2) + Complex(3.3, 4.4),
          CloseTo(Complex(4.4, 6.6), 0.0001));
      expect(Complex(3.3, 4.4) + Complex(1.1, 2.2),
          CloseTo(Complex(4.4, 6.6), 0.0001));
      expect(Complex(1.1, 2.2) + 123.toComplex(),
          CloseTo(Complex(124.1, 2.2), 0.0001));
      expect(
          Complex(1.1, 2.2) + Complex.zero, CloseTo(Complex(1.1, 2.2), 0.0001));
    });

    test('Subtraction', () {
      expect(Complex(0, 0) - Complex.zero, Complex.zero);
      expect(Complex(1.1, 2.2) - Complex(4.4, 3.3),
          CloseTo(Complex(-3.3, -1.1), 0.0001));
      expect(
          Complex(1.1, 2.2) - Complex.zero, CloseTo(Complex(1.1, 2.2), 0.0001));
    });

    test('Negative', () {
      expect(-Complex.zero, Complex.zero);
      expect(-Complex(4.4, 3.3), CloseTo(Complex(-4.4, -3.3), 0.0001));
      expect(-Complex(-4.4, -3.3), CloseTo(Complex(4.4, 3.3), 0.0001));
    });

    test('Multiplication', () {
      expect(Complex(0, 0) * Complex.zero, Complex.zero);
      expect(Complex(1.1, 2.2) * Complex.zero, Complex.zero);
      expect(Complex(1.1, 2.2) * Complex(3.3, 4.4),
          CloseTo(Complex(-6.05, 12.1), 0.0001));
      expect(Complex(3.3, 4.4) * Complex(1.1, 2.2),
          CloseTo(Complex(-6.05, 12.1), 0.0001));
      expect(
          Complex(3.3, 4.4) * Complex.one, CloseTo(Complex(3.3, 4.4), 0.0001));
      expect(Complex(3.3, 4.4) * 2.toComplex(),
          CloseTo(Complex(6.6, 8.8), 0.0001));
    });

    test('Division', () {
      final c1 = Complex(0, 0) / Complex.zero;
      expect(c1, predicate((Complex c) => c.isNaN));
      final c2 = Complex(1.1, 2.2) / Complex.zero;
      expect(c2, predicate((Complex c) => c.isNaN));
      expect(
          Complex(3.3, 4.4) / Complex.one, CloseTo(Complex(3.3, 4.4), 0.0001));
      expect(Complex(1.1, 2.2) / Complex(3.3, 4.4),
          CloseTo(Complex(0.44, 0.08), 0.0001));
      expect(Complex(3.3, 4.4) / Complex(1.1, 2.2),
          CloseTo(Complex(2.2, -0.4), 0.0001));
    });

    test('Equality', () {
      expect(Complex(0, 0), Complex.zero);
      expect(Complex(1.1, 2.2), Complex(1.1, 2.2));
      expect((1.1, 2.2), Complex(1.1, 2.2));
      expect((1, 2), Complex(1, 2));
      expect(123, Complex(123.0));
      expect(1.23, Complex(1.23));
    });

    test('Conjugate', () {
      expect(Complex.zero.conjugate(), Complex.zero);
      expect(Complex.one.conjugate(), Complex.one);
      expect(Complex(3.3, 4.4).conjugate(), Complex(3.3, -4.4));
      expect(Complex(0, 4.4).conjugate(), Complex(0, -4.4));
      expect(Complex(3.3).conjugate(), Complex(3.3));
    });

    test('Sqrt', () {
      expect(Complex.zero.sqrt(), Complex.zero);
      expect(Complex.one.sqrt(), Complex.one);
      expect(Complex(-1).sqrt(), Complex(0, 1));
      expect(
          Complex(3.3, 4.4).sqrt(), CloseTo(Complex(2.0976, 1.0488), 0.0001));
    });

    // test('Pow', () {
    //   // final c20 = Complex(2.2, 3.3).pow(2.0);
    //   // expect(c20, CloseTo(Complex(-6.05, 14.52), 0.0001));

    //   expect(Complex.zero.pow(0), Complex.one);
    //   expect(Complex.zero.pow(0.0), predicate((Complex c) => c.isNaN));
    //   expect(Complex.zero.pow(1), Complex.zero);
    //   expect(
    //     Complex.zero.pow(1.0),
    //     predicate(
    //       (Complex c) =>
    //           c == Complex.zero && !c.real.isNegative && c.imaginary.isNegative,
    //     ),
    //   );
    //   expect(Complex.zero.pow(2), Complex.zero);
    //   expect(
    //     Complex.zero.pow(2.0),
    //     predicate(
    //       (Complex c) =>
    //           c == Complex.zero && !c.real.isNegative && c.imaginary.isNegative,
    //     ),
    //   );

    //   expect(Complex(2.2, 3.3).pow(0), Complex.one);
    //   expect(Complex(2.2, 3.3).pow(0.0), Complex.one);
    //   expect(Complex(2.2, 3.3).pow(1), Complex(2.2, 3.3));
    //   expect(Complex(2.2, 3.3).pow(1.0), Complex(2.2, 3.3));
    //   expect(Complex(2.2, 3.3).pow(2), CloseTo(Complex(-6.05, 14.52), 0.0001));
    //   expect(
    //       Complex(2.2, 3.3).pow(2.0), CloseTo(Complex(-6.05, 14.52), 0.0001));

    //   expect(Complex.one.pow(4), Complex.one);
    //   expect(Complex.one.pow(-4), Complex.one);
    //   expect(Complex.one.pow(4.0), Complex.one);
    //   expect(Complex.one.pow(-4.0), Complex.one);

    //   // expect(Complex(0, 1).pow(0), Complex.one);
    //   // expect(Complex(0, 1).pow(0.0), Complex.one);
    //   // expect(Complex(0, 1).pow(1), Complex(0, 1));
    //   // expect(Complex(0, 1).pow(1.0), Complex(0, 1));
    //   // expect(Complex(0, 1).pow(2), Complex(0, 1));
    //   // expect(Complex(0, 1).pow(2.0), Complex(0, 1));

    //   expect(Complex(0, 2).pow(0), Complex.one);
    //   expect(Complex(0, 2).pow(0.0), Complex.one);
    //   expect(Complex(0, 2).pow(1), Complex(0, 2));
    //   expect(Complex(0, 2).pow(1.0), CloseTo(Complex(-8.7422E-08, 2), 0.0001));
    //   expect(Complex(0, 2).pow(2), Complex(-4, 0));
    //   expect(Complex(0, 2).pow(2.0), CloseTo(Complex(-4, -3.4969E-07), 0.0001));
    //   expect(Complex(0, 2).pow(3), Complex(0, -8));
    //   expect(Complex(0, 2).pow(3.0), CloseTo(Complex(9.5399E-08, -8), 0.0001));
    //   expect(Complex(0, 2).pow(4), Complex(16, 0));
    //   expect(Complex(0, 2).pow(4.0), CloseTo(Complex(16, 2.7975E-06), 0.0001));
    // });
  });
}

class CloseTo extends CustomMatcher {
  CloseTo(Complex expected, num delta)
      : super('Complex with', 'complex', [
          closeTo(expected.real, delta),
          closeTo(expected.imaginary, delta)
        ]);

  @override
  Object? featureValueOf(dynamic actual) =>
      [(actual as Complex).real, actual.imaginary];
}
