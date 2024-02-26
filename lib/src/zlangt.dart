import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlangt(
  final String NORM,
  final int N,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.dim();
  final D = D_.dim();
  final DU = DU_.dim();
  const ONE = 1.0, ZERO = 0.0;
  int I;
  double ANORM = 0, TEMP;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N <= 0) {
    ANORM = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    ANORM = (D[N]).abs();
    for (I = 1; I <= N - 1; I++) {
      // 10
      if (ANORM < (DL[I]).abs() || disnan((DL[I]).abs())) ANORM = (DL[I]).abs();
      if (ANORM < (D[I]).abs() || disnan((D[I]).abs())) ANORM = (D[I]).abs();
      if (ANORM < (DU[I]).abs() || disnan((DU[I]).abs())) ANORM = (DU[I]).abs();
    } // 10
  } else if (lsame(NORM, 'O') || NORM == '1') {
    // Find norm1(A).

    if (N == 1) {
      ANORM = (D[1]).abs();
    } else {
      ANORM = (D[1]).abs() + (DL[1]).abs();
      TEMP = (D[N]).abs() + (DU[N - 1]).abs();
      if (ANORM < TEMP || disnan(TEMP)) ANORM = TEMP;
      for (I = 2; I <= N - 1; I++) {
        // 20
        TEMP = (D[I]).abs() + (DL[I]).abs() + (DU[I - 1]).abs();
        if (ANORM < TEMP || disnan(TEMP)) ANORM = TEMP;
      } // 20
    }
  } else if (lsame(NORM, 'I')) {
    // Find normI(A).

    if (N == 1) {
      ANORM = (D[1]).abs();
    } else {
      ANORM = (D[1]).abs() + (DU[1]).abs();
      TEMP = (D[N]).abs() + (DL[N - 1]).abs();
      if (ANORM < TEMP || disnan(TEMP)) ANORM = TEMP;
      for (I = 2; I <= N - 1; I++) {
        // 30
        TEMP = (D[I]).abs() + (DU[I]).abs() + (DL[I - 1]).abs();
        if (ANORM < TEMP || disnan(TEMP)) ANORM = TEMP;
      } // 30
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    zlassq(N, D, 1, SCALE, SUM);
    if (N > 1) {
      zlassq(N - 1, DL, 1, SCALE, SUM);
      zlassq(N - 1, DU, 1, SCALE, SUM);
    }
    ANORM = SCALE.value * sqrt(SUM.value);
  }

  return ANORM;
}
