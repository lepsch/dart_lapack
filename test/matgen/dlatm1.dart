import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

import 'dlaran.dart';

void dlatm1(
  final int MODE,
  final double COND,
  final int IRSIGN,
  final int IDIST,
  final Array<int> ISEED_,
  final Array<double> D_,
  final int N,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having();
  final D = D_.having();
  const ONE = 1.0;
  const HALF = 0.5;
  int I;
  double ALPHA, TEMP;

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
  } else if ((MODE == 6 || MODE == -6) && (IDIST < 1 || IDIST > 3)) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -7;
  }

  if (INFO.value != 0) {
    xerbla('DLATM1', -INFO.value);
    return;
  }

  if (MODE == 0) return;

  // Compute D according to COND and MODE
  switch (MODE.abs()) {
    case 1:
      // One large D value:
      for (I = 1; I <= N; I++) {
        D[I] = ONE / COND;
      }
      D[1] = ONE;
      break;

    case 2:
      // One small D value:
      for (I = 1; I <= N; I++) {
        D[I] = ONE;
      }
      D[N] = ONE / COND;
      break;

    case 3:
      // Exponentially distributed D values:
      D[1] = ONE;
      if (N > 1) {
        ALPHA = pow(COND, -ONE / (N - 1)).toDouble();
        for (I = 2; I <= N; I++) {
          D[I] = pow(ALPHA, I - 1).toDouble();
        }
      }
      break;

    case 4:
      // Arithmetically distributed D values:
      D[1] = ONE;
      if (N > 1) {
        TEMP = ONE / COND;
        ALPHA = (ONE - TEMP) / (N - 1).toDouble();
        for (I = 2; I <= N; I++) {
          D[I] = (N - I).toDouble() * ALPHA + TEMP;
        }
      }
      break;

    case 5:
      // Randomly distributed D values on ( 1/COND , 1):
      ALPHA = log(ONE / COND);
      for (I = 1; I <= N; I++) {
        D[I] = exp(ALPHA * dlaran(ISEED));
      }
      break;

    case 6:
      // Randomly distributed D values from IDIST
      dlarnv(IDIST, ISEED, N, D);
      break;
  }

  // If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
  // random signs to D

  if ((MODE != -6 && MODE != 0 && MODE != 6) && IRSIGN == 1) {
    for (I = 1; I <= N; I++) {
      TEMP = dlaran(ISEED);
      if (TEMP > HALF) D[I] = -D[I];
    }
  }

  // Reverse if MODE < 0

  if (MODE < 0) {
    for (I = 1; I <= N / 2; I++) {
      TEMP = D[I];
      D[I] = D[N + 1 - I];
      D[N + 1 - I] = TEMP;
    }
  }
}
