import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarnv.dart';

import 'dlaran.dart';
import 'zlarnd.dart';

void zlatm1(
  final int MODE,
  final double COND,
  final int IRSIGN,
  final int IDIST,
  final Array<int> ISEED_,
  final Array<Complex> D_,
  final int N,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.dim(4);
  final D = D_.dim();

  const ONE = 1.0;
  int I;
  double ALPHA, TEMP;
  Complex CTEMP;

  // Decode and Test the input parameters. Initialize flags & seed.

  INFO.value = 0;

  // Quick return if possible

  if (N == 0) return;

  // Set INFO.value if an error

  if (MODE < -6 || MODE > 6) {
    INFO.value = -1;
  } else if ((MODE != -6 && MODE != 0 && MODE != 6) &&
      (IRSIGN != 0 && IRSIGN != 1)) {
    INFO.value = -2;
  } else if ((MODE != -6 && MODE != 0 && MODE != 6) && COND < ONE) {
    INFO.value = -3;
  } else if ((MODE == 6 || MODE == -6) && (IDIST < 1 || IDIST > 4)) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -7;
  }

  if (INFO.value != 0) {
    xerbla('ZLATM1', -INFO.value);
    return;
  }

  // Compute D according to COND and MODE

  if (MODE != 0) {
    //  GO TO ( 10, 30, 50, 70, 90, 110 )( MODE ).abs();
    switch (MODE.abs()) {
      case 1:
        // One large D value:

        for (I = 1; I <= N; I++) {
          // 20
          D[I] = (ONE / COND).toComplex();
        } // 20
        D[1] = Complex.one;
        break;

      case 2:

        // One small D value:

        for (I = 1; I <= N; I++) {
          // 40
          D[I] = Complex.one;
        } // 40
        D[N] = (ONE / COND).toComplex();
        break;

      case 3:

        // Exponentially distributed D values:

        D[1] = Complex.one;
        if (N > 1) {
          ALPHA = pow(COND, (-ONE / (N - 1).toDouble())).toDouble();
          for (I = 2; I <= N; I++) {
            // 60
            D[I] = pow(ALPHA, I - 1).toComplex();
          } // 60
        }
        break;

      case 4:

        // Arithmetically distributed D values:

        D[1] = Complex.one;
        if (N > 1) {
          TEMP = ONE / COND;
          ALPHA = (ONE - TEMP) / (N - 1).toDouble();
          for (I = 2; I <= N; I++) {
            // 80
            D[I] = ((N - I).toDouble() * ALPHA + TEMP).toComplex();
          } // 80
        }
        break;

      case 5:

        // Randomly distributed D values on ( 1/COND , 1):

        ALPHA = log(ONE / COND);
        for (I = 1; I <= N; I++) {
          // 100
          D[I] = exp(ALPHA * dlaran(ISEED)).toComplex();
        } // 100
        break;

      case 6:

        // Randomly distributed D values from IDIST

        zlarnv(IDIST, ISEED, N, D);
        break;
    } // 120

    // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
    // random signs to D

    if ((MODE != -6 && MODE != 0 && MODE != 6) && IRSIGN == 1) {
      for (I = 1; I <= N; I++) {
        // 130
        CTEMP = zlarnd(3, ISEED);
        D[I] = D[I] * (CTEMP / CTEMP.abs().toComplex());
      } // 130
    }

    // Reverse if MODE < 0

    if (MODE < 0) {
      for (I = 1; I <= N / 2; I++) {
        // 140
        CTEMP = D[I];
        D[I] = D[N + 1 - I];
        D[N + 1 - I] = CTEMP;
      } // 140
    }
  }
}
