import 'package:lapack/src/matrix.dart';

double dsdot(
  final int N,
  final Array<double> SX,
  final int INCX,
  final Array<double> SY,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// Authors:
// ========
// Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
// Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)

  int I, KX, KY, NS;
  var result = 0.0;
  if (N <= 0) return result;
  if (INCX == INCY && INCX > 0) {
    // Code for equal, positive, non-unit increments.

    NS = N * INCX;
    for (I = 1; INCX < 0 ? I >= NS : I <= NS; I += INCX) {
      result += SX[I].toDouble() * SY[I].toDouble();
    }
  } else {
    // Code for unequal or nonpositive increments.

    KX = 1;
    KY = 1;
    if (INCX < 0) KX = 1 + (1 - N) * INCX;
    if (INCY < 0) KY = 1 + (1 - N) * INCY;
    for (I = 1; I <= N; I++) {
      result += SX[KX].toDouble() * SY[KY].toDouble();
      KX = KX + INCX;
      KY = KY + INCY;
    }
  }
  return result;
}
