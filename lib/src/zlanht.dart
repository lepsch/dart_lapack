import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlassq.dart';

double zlanht(
  final String NORM,
  final int N,
  final Array<double> D_,
  final Array<Complex> E_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.dim();
  final E = E_.dim();
  const ONE = 1.0, ZERO = 0.0;
  int I;
  double ANORM = 0;
  final SCALE = Box(0.0), SUM = Box(0.0);

  if (N <= 0) {
    ANORM = ZERO;
  } else if (lsame(NORM, 'M')) {
    // Find max(abs(A(i,j))).

    ANORM = D[N].abs();
    for (I = 1; I <= N - 1; I++) {
      // 10
      SUM.value = D[I].abs();
      if (ANORM < SUM.value || disnan(SUM.value)) ANORM = SUM.value;
      SUM.value = E[I].abs();
      if (ANORM < SUM.value || disnan(SUM.value)) ANORM = SUM.value;
    } // 10
  } else if (lsame(NORM, 'O') || NORM == '1' || lsame(NORM, 'I')) {
    // Find norm1(A).

    if (N == 1) {
      ANORM = D[1].abs();
    } else {
      ANORM = D[1].abs() + E[1].abs();
      SUM.value = E[N - 1].abs() + D[N].abs();
      if (ANORM < SUM.value || disnan(SUM.value)) ANORM = SUM.value;
      for (I = 2; I <= N - 1; I++) {
        // 20
        SUM.value = D[I].abs() + E[I].abs() + E[I - 1].abs();
        if (ANORM < SUM.value || disnan(SUM.value)) ANORM = SUM.value;
      } // 20
    }
  } else if ((lsame(NORM, 'F')) || (lsame(NORM, 'E'))) {
    // Find normF(A).

    SCALE.value = ZERO;
    SUM.value = ONE;
    if (N > 1) {
      zlassq(N - 1, E, 1, SCALE, SUM);
      SUM.value = 2 * SUM.value;
    }
    dlassq(N, D, 1, SCALE, SUM);
    ANORM = SCALE.value * sqrt(SUM.value);
  }

  return ANORM;
}
