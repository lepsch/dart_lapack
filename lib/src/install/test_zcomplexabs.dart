import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/intrinsics/digits.dart';
import 'package:lapack/src/intrinsics/epsilon.dart';
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
  const threeFourth = 3.0 / 4, fiveFourth = 5.0 / 4, oneHalf = 1.0 / 2;
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
      nFailingTests,
      nTests;
  double R, answerC, answerD, aInf, aNaN, relDiff, b, eps, blueMin, blueMax, Xj;
  final X = Array<double>(N), stepX = Array<double>(N), limX = Array<double>(N);
  Complex Y;
  final cInf = Array<Complex>(nInf), cNaN = Array<Complex>(nNaN);

  // .. Intrinsic Functions ..
  // intrinsic ABS, DBLE, RADIX, CEILING, TINY, DIGITS, SQRT, MAXEXPONENT, MINEXPONENT, FLOOR, HUGE, DCMPLX, EPSILON

  // .. Initialize error counts ..
  subnormalTreatedAs0 = 0;
  caseAFails = 0;
  caseBFails = 0;
  caseCFails = 0;
  caseDFails = 0;
  caseEFails = 0;
  caseFFails = 0;
  nFailingTests = 0;
  nTests = 0;

  // .. Initialize machine constants ..
  min = minexponent(0.0);
  Max = maxexponent(0.0);
  m = digits(0.0);
  b = (radix(0.0)).toDouble();
  eps = epsilon(0.0);
  blueMin = pow(b, ((min - 1) * 0.5).ceil()).toDouble();
  blueMax = pow(b, ((Max - m + 1) * 0.5).floor()).toDouble();

  // .. Vector X ..
  X[1] = tiny(0.0) * pow(b, (1 - m).toDouble());
  X[2] = tiny(0.0);
  X[3] = huge(0.0);
  X[4] = pow(b, (Max - 1).toDouble()).toDouble();

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
  aInf = X[3] * 2;
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
    print('# X :=$X');
    print('# Blue min constant :=$blueMin');
    print('# Blue max constant :=$blueMax');
  }

  Xj = X[1];
  if (Xj == 0.0) {
    subnormalTreatedAs0++;
    if (debug || subnormalTreatedAs0 == 1) {
      print('!! fl( subnormal ) may be 0');
    }
  } else {
    for (i = 1; i <= N; i++) {
      Xj = X[i];
      if (Xj == 0.0) {
        subnormalTreatedAs0++;
        if (debug || subnormalTreatedAs0 == 1) {
          print('!! fl( subnormal ) may be 0');
        }
      }
    }
  }

  // Test (a) y = x + 0 * I, |y| = x
  for (i = 1; i <= N; i++) {
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
        R = Y.abs();
        if (R != Xj) {
          caseAFails++;
          if (caseAFails == 1) {
            print('!! Some (x+0*I).abs() differ from x.abs()');
          }
          _print9999('a', i, Xj, '(1+0*I)', R, Xj);
        }
        Xj *= stepX[i];
      }
    }
  }

  // Test (b) y = 0 + x * I, |y| = x
  for (i = 1; i <= N; i++) {
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
        R = Y.abs();
        if (R != Xj) {
          caseBFails++;
          if (caseBFails == 1) {
            print('!! Some (0+x*I).abs() differ from x.abs()');
          }
          _print9999('b', i, Xj, '(0+1*I)', R, Xj);
        }
        Xj *= stepX[i];
      }
    }
  }

  // Test (c) y = (3/4)*x + x * I, |y| = (5/4)*x
  for (i = 1; i <= N; i++) {
    if (i == 3) continue;
    if (i == 1) {
      Xj = 4 * X[i];
    } else {
      Xj = X[i];
    }
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [c] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        nTests++;
        answerC = fiveFourth * Xj;
        Y = Complex(threeFourth * Xj, Xj);
        R = Y.abs();
        if (R != answerC) {
          caseCFails++;
          if (caseCFails == 1) {
            print('!! Some (x*(3/4+I).abs()) differ from (5/4)*x.abs()');
          }
          _print9999('c', i, Xj, '(3/4+I)', R, answerC);
        }
        Xj *= stepX[i];
      }
    }
  }

  // Test (d) y = (1/2)*x + (1/2)*x * I, |y| = (1/2)*x*sqrt(2)
  for (i = 1; i <= N; i++) {
    if (i == 1) {
      Xj = 2 * X[i];
    } else {
      Xj = X[i];
    }
    if (Xj == 0.0) {
      subnormalTreatedAs0++;
      if (debug || subnormalTreatedAs0 == 1) {
        print('!! [d] fl( subnormal ) may be 0');
      }
    } else {
      while (Xj != limX[i]) {
        answerD = (oneHalf * Xj) * sqrt(2.0);
        if (answerD == 0.0) {
          subnormalTreatedAs0++;
          if (debug || subnormalTreatedAs0 == 1) {
            print('!! [d] fl( subnormal ) may be 0');
          }
        } else {
          nTests++;
          Y = Complex(oneHalf * Xj, oneHalf * Xj);
          R = Y.abs();
          relDiff = (R - answerD).abs() / answerD;
          if (relDiff >= (0.5 * eps)) {
            caseDFails++;
            if (caseDFails == 1) {
              print('!! Some (x*(1+I).abs()) differ from sqrt(2)*x.abs()');
            }
            _print9999('d', i, (oneHalf * Xj), '(1+1*I)', R, answerD);
          }
        }
        Xj *= stepX[i];
      }
    }
  }

  // Test (e) Infs
  for (i = 1; i <= nInf; i++) {
    nTests++;
    Y = cInf[i];
    R = Y.abs();
    if (!(R > huge(0.0))) {
      caseEFails++;
      print(
          '[i$i] ABS(${Y.real.e8_1}${Y.imaginary.sp}${Y.imaginary.e8_1}*I ) = ${R.e8_1} differs from Inf');
    }
  }

  // Test (f) NaNs
  for (i = 1; i <= nNaN; i++) {
    nTests++;
    Y = cNaN[i];
    R = Y.abs();
    if (R == R) {
      caseFFails++;
      print(
          '[n$i] ABS(${Y.real.e8_1}${Y.imaginary.sp}${Y.imaginary.e8_1}*I ) = ${R.e8_1} differs from NaN');
    }
  }

  // If any test fails, displays a message
  nFailingTests = caseAFails +
      caseBFails +
      caseCFails +
      caseDFails +
      caseEFails +
      caseFFails;
  if (nFailingTests > 0) {
    print(
        '# ${nTests - nFailingTests} tests out of $nTests pass for (a+b*I).abs(),$nFailingTests tests fail.');
  } else {
    print('# All tests pass for (a+b*I).abs()');
  }

  // If anything was written to stderr, print the message
  if ((caseAFails > 0) ||
      (caseBFails > 0) ||
      (caseCFails > 0) ||
      (caseDFails > 0)) {
    print('# Please check the failed (a+b*I).abs() in [stderr]');
  }

  // .. Formats ..
//  9997 FORMAT(  );

//  9998 FORMAT(  );

//  9999 FORMAT(  );
}

void _print9999(
    String s1, int i, double d1, String s2, double d2, double answer) {
  print(
      '[${s1.a1}${i.i1}] (${d1.e24_16e3} * ${s2.a7} ).abs() = ${d2.e24_16e3} differs from ${answer.e24_16e3}');
}
