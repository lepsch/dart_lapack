import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zpoequb(
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> S_,
  final Box<double> SCOND,
  final Box<double> AMAX,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final S = S_.dim();
  const ZERO = 0.0, ONE = 1.0;
  int I;
  double SMIN, BASE, TMP;

  // Test the input parameters.

  // Positive definite only performs 1 pass of equilibration.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (LDA < max(1, N)) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('ZPOEQUB', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) {
    SCOND.value = ONE;
    AMAX.value = ZERO;
    return;
  }

  BASE = dlamch('B');
  TMP = -0.5 / log(BASE);

  // Find the minimum and maximum diagonal elements.

  S[1] = A[1][1].toDouble();
  SMIN = S[1];
  AMAX.value = S[1];
  for (I = 2; I <= N; I++) {
    // 10
    S[I] = A[I][I].toDouble();
    SMIN = min(SMIN, S[I]);
    AMAX.value = max(AMAX.value, S[I]);
  } // 10

  if (SMIN <= ZERO) {
    // Find the first non-positive diagonal element and return.

    for (I = 1; I <= N; I++) {
      // 20
      if (S[I] <= ZERO) {
        INFO.value = I;
        return;
      }
    } // 20
  } else {
    // Set the scale factors to the reciprocals
    // of the diagonal elements.

    for (I = 1; I <= N; I++) {
      // 30
      S[I] = pow(BASE, (TMP * log(S[I])).toInt()).toDouble();
    } // 30

    // Compute SCOND.value = min(S(I)) / max(S(I)).

    SCOND.value = sqrt(SMIN) / sqrt(AMAX.value);
  }
}
