import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/intrinsics/huge.dart';
import 'package:lapack/src/matrix.dart';

void main() {
// -- LAPACK test routine --
// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  const nNaN = 3, nInf = 5;
  const czero = Complex.zero, cone = Complex.one;
  int i, nFailingTests, nTests;
  double aInf, aNaN, OV;
  Complex Y, R;
  final cInf = Array<Complex>(nInf), cNaN = Array<Complex>(nNaN);

  // .. Intrinsic Functions ..
  // intrinsic HUGE, DCMPLX

  // .. Initialize error counts ..
  nFailingTests = 0;
  nTests = 0;

  // .. Inf entries ..
  OV = huge(0.0);
  aInf = OV * 2;
  cInf[1] = Complex(aInf, 0.0);
  cInf[2] = Complex(-aInf, 0.0);
  cInf[3] = Complex(0.0, aInf);
  cInf[4] = Complex(0.0, -aInf);
  cInf[5] = Complex(aInf, aInf);

  // .. NaN entries ..
  aNaN = aInf / aInf;
  cNaN[1] = Complex(aNaN, 0.0);
  cNaN[2] = Complex(0.0, aNaN);
  cNaN[3] = Complex(aNaN, aNaN);

  // .. Tests ..

  // Test (a) Infs
  for (i = 1; i <= nInf; i++) {
    nTests = nTests + 3;
    Y = cInf[i];
    R = czero * Y;
    if (R == R) {
      nFailingTests++;
      _print9998('ia', i, czero, Y, R, 'NaN');
    }
    R = cone * Y;
    if ((R != Y) && (R == R)) {
      nFailingTests++;
      _print9998('ib', i, cone, Y, R, 'the input and NaN');
    }
    R = Y * Y;
    if ((i == 1) || (i == 2)) {
      if ((R != cInf[1]) && (R == R)) {
        nFailingTests++;
        _print9998('ic', i, Y, Y, R, 'Inf and NaN');
      }
    } else if ((i == 3) || (i == 4)) {
      if ((R != cInf[2]) && (R == R)) {
        nFailingTests++;
        _print9998('ic', i, Y, Y, R, '-Inf and NaN');
      }
    } else {
      if (R == R) {
        nFailingTests++;
        _print9998('ic', i, Y, Y, R, 'NaN');
      }
    }
  }

  // Test (b) NaNs
  for (i = 1; i <= nNaN; i++) {
    nTests = nTests + 3;
    Y = cNaN[i];
    R = czero * Y;
    if (R == R) {
      nFailingTests++;
      _print9998('na', i, czero, Y, R, 'NaN');
    }
    R = cone * Y;
    if (R == R) {
      nFailingTests++;
      _print9998('nb', i, cone, Y, R, 'NaN');
    }
    R = Y * Y;
    if (R == R) {
      nFailingTests++;
      _print9998('nc', i, Y, Y, R, 'NaN');
    }
  }

  if (nFailingTests > 0) {
    print(
        '# ${nTests - nFailingTests} tests out of $nTests pass for Complex multiplication,$nFailingTests fail.');
  } else {
    print('# All tests pass for Complex multiplication.');
  }
}

void _print9998(
    String s1, int i, Complex c1, Complex c2, Complex c3, String s2) {
  print(
      '[${s1.a2}${i.i1}] (${c1.real.e24_16e3}${c1.imaginary.sp}${c1.imaginary.e24_16e3}*I) * (${c2.real.e24_16e3}${c2.imaginary.sp}${c2.imaginary.e24_16e3}*I) = (${c3.real.e24_16e3}${c3.imaginary.sp}${c3.imaginary.e24_16e3}*I) differs from ${s2.a17}');
}
