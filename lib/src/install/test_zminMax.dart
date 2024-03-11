import 'dart:math';

import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/intrinsics/huge.dart';
import 'package:lapack/src/matrix.dart';

void main() {
// -- LAPACK test routine --
// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  const zero = 0.0;
  int i, nFailingTests, nTests;
  double aInf, aNaN, OV, R;

  // .. Initialize error counts ..
  nFailingTests = 0;
  nTests = 0;

  // .. Inf and NaN entries ..
  OV = huge(0.0);
  aInf = OV * 2;
  aNaN = aInf / aInf;
  final X = Array.fromList([-aInf, zero, -aInf, zero, aInf, aInf, zero, aNaN]);
  final Y = Array.fromList([zero, aInf, aInf, -aInf, zero, -aInf, aNaN, zero]);

  // .. Tests ..

  for (i = 1; i <= 3; i++) {
    nTests = nTests + 2;
    R = min(X[i], Y[i]);
    if (R != X[i]) {
      nFailingTests++;
      _print9998('i', i, 'MIN', X[i], Y[i], R);
    }
    R = max(X[i], Y[i]);
    if (R != Y[i]) {
      nFailingTests++;
      _print9998('i', i, 'MAX', X[i], Y[i], R);
    }
  }
  for (i = 4; i <= 6; i++) {
    nTests = nTests + 2;
    R = min(X[i], Y[i]);
    if (R != Y[i]) {
      nFailingTests++;
      _print9998('i', i, 'MIN', X[i], Y[i], R);
    }
    R = max(X[i], Y[i]);
    if (R != X[i]) {
      nFailingTests++;
      _print9998('i', i, 'MAX', X[i], Y[i], R);
    }
  }
  for (i = 7; i <= 8; i++) {
    nTests = nTests + 2;
    R = min(X[i], Y[i]);
    if (R == R) {
      nFailingTests++;
      _print9998('i', i, 'MIN', X[i], Y[i], R);
    }
    R = max(X[i], Y[i]);
    if (R == R) {
      nFailingTests++;
      _print9998('i', i, 'MAX', X[i], Y[i], R);
    }
  }

  if (nFailingTests > 0) {
    print(
        '# ${nTests - nFailingTests} tests out of $nTests pass for intrinsic MIN and MAX,$nFailingTests fail.');
  } else {
    print('# All tests pass for intrinsic MIN and MAX.');
  }
}

void _print9998(String s, int i, String f, double d1, double d2, double d3) {
  print('[${s.a1}${i.i1}] ${f.a3}(${d1.f5_0},${d2.f5_0}) = ${d3.f5_0}');
}
