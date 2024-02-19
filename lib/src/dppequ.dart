import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dppequ(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> S_,
  final Box<double> SCOND,
  final Box<double> AMAX,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.dim();
  final S = S_.dim();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int I, JJ;
  double SMIN;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('DPPEQU', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) {
    SCOND.value = ONE;
    AMAX.value = ZERO;
    return;
  }

  // Initialize SMIN and AMAX.value.

  S[1] = AP[1];
  SMIN = S[1];
  AMAX.value = S[1];

  if (UPPER) {
    // UPLO = 'U':  Upper triangle of A is stored.
    // Find the minimum and maximum diagonal elements.

    JJ = 1;
    for (I = 2; I <= N; I++) {
      // 10
      JJ = JJ + I;
      S[I] = AP[JJ];
      SMIN = min(SMIN, S[I]);
      AMAX.value = max(AMAX.value, S[I]);
    } // 10
  } else {
    // UPLO = 'L':  Lower triangle of A is stored.
    // Find the minimum and maximum diagonal elements.

    JJ = 1;
    for (I = 2; I <= N; I++) {
      // 20
      JJ = JJ + N - I + 2;
      S[I] = AP[JJ];
      SMIN = min(SMIN, S[I]);
      AMAX.value = max(AMAX.value, S[I]);
    } // 20
  }

  if (SMIN <= ZERO) {
    // Find the first non-positive diagonal element and return.

    for (I = 1; I <= N; I++) {
      // 30
      if (S[I] <= ZERO) {
        INFO.value = I;
        return;
      }
    } // 30
  } else {
    // Set the scale factors to the reciprocals
    // of the diagonal elements.

    for (I = 1; I <= N; I++) {
      // 40
      S[I] = ONE / sqrt(S[I]);
    } // 40

    // Compute SCOND.value = min(S(I)) / max(S(I))

    SCOND.value = sqrt(SMIN) / sqrt(AMAX.value);
  }
}
