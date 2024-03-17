import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/huge.dart';
import 'package:lapack/src/intrinsics/maxexponent.dart';
import 'package:lapack/src/intrinsics/minexponent.dart';
import 'package:lapack/src/intrinsics/radix.dart';
import 'package:lapack/src/intrinsics/tiny.dart';
import 'package:lapack/src/matrix.dart';

void main() {
// -- LAPACK test routine --
// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  const debug = false;
  const N = 4, nNaN = 3, nInf = 5;
  // const threeFourth = 3.0 / 4, fiveFourth = 5.0 / 4;
  int i,
      min,
      Max,
      m,
      subnormalTreatedAs0,
      caseAFails,
      caseBFails,
      caseCFails,
      caseDFails,
      caseEFails,
      caseFFails,
      caseInfFails,
      caseNaNFails,
      nFailingTests,
      nTests;
  double aInf,
      aNaN,
      b,
      // eps,
      blueMin,
      blueMax,
      OV,
      Xj;
  Complex Y, Y2, R;
  final cInf = Array<Complex>(nInf), cNaN = Array<Complex>(nNaN);
  final X = Array<double>(N), stepX = Array<double>(N), limX = Array<double>(N);

  // .. Initialize error counts ..
  subnormalTreatedAs0 = 0;
  caseAFails = 0;
  caseBFails = 0;
  caseCFails = 0;
  caseDFails = 0;
  caseEFails = 0;
  caseFFails = 0;
  caseInfFails = 0;
  caseNaNFails = 0;
  nFailingTests = 0;
  nTests = 0;

  // .. Initialize machine constants ..
  min = minexponent(0.0);
  Max = maxexponent(0.0);
  m = digits(0.0);
  b = (radix(0.0)).toDouble();
  // eps = epsilon(0.0);
  blueMin = pow(b, ((min - 1) * 0.5).ceil()).toDouble();
  blueMax = pow(b, ((Max - m + 1) * 0.5).floor()).toDouble();
  OV = huge(0.0);

  // .. Vector X ..
  X[1] = tiny(0.0) * pow(b, (1 - m).toDouble());
  X[2] = tiny(0.0);
  X[3] = OV;
  X[4] = pow(b, (Max - 1)).toDouble();

  // .. Then modify X using the step ..
  stepX[1] = 2.0;
  stepX[2] = 2.0;
  stepX[3] = 0.0;
  stepX[4] = 0.5;

  // .. Up to the value ..
  limX[1] = X[2];
  limX[2] = 1.0;
  limX[3] = 0.0;
  limX[4] = 2.0;

  // .. Inf entries ..
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

  // ignore: dead_code
  if (debug) {
    print('# X := $X');
    print('# Blue min constant := $blueMin');
    print('# Blue max constant := $blueMax');
  }

  Xj = X[1];
  if (Xj == 0.0) {
    subnormalTreatedAs0++;
    if (debug || subnormalTreatedAs0 == 1) {
      print('!! fl( subnormal ) may be 0');
    }
  } else {
    for (i = 1; i <= N; i++) {
      // 100
      Xj = X[i];
      if (Xj == 0.0) {
        subnormalTreatedAs0++;
        if (debug || subnormalTreatedAs0 == 1) {
          print('!! fl( subnormal ) may be 0');
        }
      }
    } // 100
  }

  // Test (a) y = x + 0 * I, y/y = 1
  for (i = 1; i <= N; i++) {
    // 10
    Xj = X[i];
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [a] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        nTests++;
        Y = Complex(Xj, 0.0);
        R = Y / Y;
        if (R != Complex.one) {
          caseAFails++;
          if (caseAFails == 1) {
            print('!! Some (x+0*I)/(x+0*I) differ from 1');
          }
          _print9999('a', i, Xj, '(x+0*I)/(x+0*I)', R, Complex.one);
        }
        Xj *= stepX[i];
      }
    }
  } // 10

  // Test (b) y = 0 + x * I, y/y = 1
  for (i = 1; i <= N; i++) {
    // 20
    Xj = X[i];
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [b] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        nTests++;
        Y = Complex(0.0, Xj);
        R = Y / Y;
        if (R != Complex.one) {
          caseBFails++;
          if (caseBFails == 1) {
            print('!! Some (0+x*I)/(0+x*I) differ from 1');
          }
          _print9999('b', i, Xj, '(0+x*I)/(0+x*I)', R, Complex.one);
        }
        Xj *= stepX[i];
      }
    }
  } // 20

  // Test (c) y = x + x * I, y/y = 1
  for (i = 1; i <= N; i++) {
    // 30
    Xj = X[i];
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [c] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        nTests++;
        Y = Complex(Xj, Xj);
        R = Y / Y;
        if (R != Complex.one) {
          caseCFails++;
          if (caseCFails == 1) {
            print('!! Some (x+x*I)/(x+x*I) differ from 1');
          }
          _print9999('c', i, Xj, '(x+x*I)/(x+x*I)', R, Complex.one);
        }
        Xj *= stepX[i];
      }
    }
  } // 30

  // Test (d) y1 = 0 + x * I, y2 = x + 0 * I, y1/y2 = I
  for (i = 1; i <= N; i++) {
    // 40
    Xj = X[i];
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [d] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        nTests++;
        Y = Complex(0.0, Xj);
        Y2 = Complex(Xj, 0.0);
        R = Y / Y2;
        if (R != Complex(0.0, 1.0)) {
          caseDFails++;
          if (caseDFails == 1) {
            print('!! Some (0+x*I)/(x+0*I) differ from I');
          }
          _print9999('d', i, Xj, '(0+x*I)/(x+0*I)', R, Complex(0.0, 1.0));
        }
        Xj *= stepX[i];
      }
    }
  } // 40

  // Test (e) y1 = 0 + x * I, y2 = x + 0 * I, y2/y1 = -I
  for (i = 1; i <= N; i++) {
    // 50
    Xj = X[i];
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [e] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        nTests++;
        Y = Complex(0.0, Xj);
        Y2 = Complex(Xj, 0.0);
        R = Y2 / Y;
        if (R != Complex(0.0, -1.0)) {
          caseEFails++;
          if (caseEFails == 1) {
            print('!! Some (x+0*I)/(0+x*I) differ from -I');
          }
          _print9999('e', i, Xj, '(x+0*I)/(0+x*I)', R, Complex(0.0, -1.0));
        }
        Xj *= stepX[i];
      }
    }
  } // 50

  // Test (f) y = x + x * I, y/conj(y) = I
  for (i = 1; i <= N; i++) {
    // 60
    Xj = X[i];
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [f] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        nTests++;
        Y = Complex(Xj, Xj);
        R = Y / Y.conjugate();
        if (R != Complex(0.0, 1.0)) {
          caseFFails++;
          if (caseFFails == 1) {
            print('!! Some (x+x*I)/(x-x*I) differ from I');
          }
          _print9999('f', i, Xj, '(x+x*I)/(x-x*I)', R, Complex(0.0, 1.0));
        }
        Xj *= stepX[i];
      }
    }
  } // 60

  // Test (g) Infs
  for (i = 1; i <= nInf; i++) {
    // 70
    nTests += 3;
    Y = cInf[i];
    R = Complex.zero / Y;
    if ((R != Complex.zero) && (R == R)) {
      caseInfFails++;
      _print9998('ia', i, Complex.zero, Y, R, 'NaN and 0');
    }
    R = Complex.one / Y;
    if ((R != Complex.zero) && (R == R)) {
      caseInfFails++;
      _print9998('ib', i, Complex.one, Y, R, 'NaN and 0');
    }
    R = Y / Y;
    if (R == R) {
      caseInfFails++;
      _print9998('ic', i, Y, Y, R, 'NaN');
    }
  } // 70

  // Test (h) NaNs
  for (i = 1; i <= nNaN; i++) {
    // 80
    nTests += 3;
    Y = cNaN[i];
    R = Complex.zero / Y;
    if (R == R) {
      caseNaNFails++;
      _print9998('na', i, Complex.zero, Y, R, 'NaN');
    }
    R = Complex.one / Y;
    if (R == R) {
      caseNaNFails++;
      _print9998('nb', i, Complex.one, Y, R, 'NaN');
    }
    R = Y / Y;
    if (R == R) {
      caseNaNFails++;
      _print9998('nc', i, Y, Y, R, 'NaN');
    }
  } // 80

  // If any test fails, displays a message
  nFailingTests = caseAFails +
      caseBFails +
      caseCFails +
      caseDFails +
      caseEFails +
      caseFFails +
      caseInfFails +
      caseNaNFails;
  if (nFailingTests > 0) {
    print(
        '# ${nTests - nFailingTests} tests out of $nTests pass for Complex division, $nFailingTests fail.');
  } else {
    print('# All tests pass for Complex division.');
  }

  // If anything was written to stderr, print the message
  if ((caseAFails > 0) ||
      (caseBFails > 0) ||
      (caseCFails > 0) ||
      (caseDFails > 0) ||
      (caseEFails > 0) ||
      (caseFFails > 0)) {
    print('# Please check the failed divisions in [stderr]');
  }
}

void _print9999(String s1, int i, double v, String s2, Complex c1, Complex c2) {
  print(
      '[${s1.a2}${i.i1}] X = ${v.e24_16e3} : ${s2.a15} = ${c1.real.e24_16e3}${c1.imaginary.sp}${c1.imaginary.e24_16e3}*I differs from ${c2.real.e24_16e3}${c2.imaginary.sp}${c2.imaginary.e24_16e3}*I');
}

void _print9998(
    String s1, int i, Complex c1, Complex c2, Complex c3, String s2) {
  print(
      '[${s1.a2}${i.i1}] ${c1.real.e24_16e3}${c1.imaginary.sp}${c1.imaginary.e24_16e3}*I * ${c2.real.e24_16e3}${c2.imaginary.sp}${c2.imaginary.e24_16e3}*I = ${c3.real.e24_16e3}${c3.imaginary.sp}${c3.imaginary.e24_16e3}*I differs from ${s2.a10}');
}
