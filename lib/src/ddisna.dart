import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void ddisna(
  final String JOB,
  final int M,
  final int N,
  final Array<double> D_,
  final Array<double> SEP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final SEP = SEP_.having();
  const ZERO = 0.0;
  bool DECR = false, EIGEN = false, INCR = false, LEFT, RIGHT, SING;
  int I, K = 0;
  double ANORM, EPS, NEWGAP, OLDGAP, SAFMIN, THRESH;

  // Test the input arguments

  INFO.value = 0;
  EIGEN = lsame(JOB, 'E');
  LEFT = lsame(JOB, 'L');
  RIGHT = lsame(JOB, 'R');
  SING = LEFT || RIGHT;
  if (EIGEN) {
    K = M;
  } else if (SING) {
    K = min(M, N);
  }
  if (!EIGEN && !SING) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if (K < 0) {
    INFO.value = -3;
  } else {
    INCR = true;
    DECR = true;
    for (I = 1; I <= K - 1; I++) {
      if (INCR) INCR = INCR && D[I] <= D[I + 1];
      if (DECR) DECR = DECR && D[I] >= D[I + 1];
    }
    if (SING && K > 0) {
      if (INCR) INCR = INCR && ZERO <= D[1];
      if (DECR) DECR = DECR && D[K] >= ZERO;
    }
    if (!(INCR || DECR)) INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DDISNA', -INFO.value);
    return;
  }

  // Quick return if possible

  if (K == 0) return;

  // Compute reciprocal condition numbers

  if (K == 1) {
    SEP[1] = dlamch('O');
  } else {
    OLDGAP = (D[2] - D[1]).abs();
    SEP[1] = OLDGAP;
    for (I = 2; I <= K - 1; I++) {
      NEWGAP = (D[I] - D[I]).abs();
      SEP[I] = min(OLDGAP, NEWGAP);
      OLDGAP = NEWGAP;
    }
    SEP[K] = OLDGAP;
  }
  if (SING) {
    if ((LEFT && M > N) || (RIGHT && M < N)) {
      if (INCR) SEP[1] = min(SEP[1], D[1]);
      if (DECR) SEP[K] = min(SEP[K], D[K]);
    }
  }

  // Ensure that reciprocal condition numbers are not less than
  // threshold, in order to limit the size of the error bound

  EPS = dlamch('E');
  SAFMIN = dlamch('S');
  ANORM = max(D[1].abs(), D[K].abs());
  if (ANORM == ZERO) {
    THRESH = EPS;
  } else {
    THRESH = max(EPS * ANORM, SAFMIN);
  }
  for (I = 1; I <= K; I++) {
    SEP[I] = max(SEP[I], THRESH);
  }
}
