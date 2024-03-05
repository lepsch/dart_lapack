import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

double dzsum1(final int N, final Array<Complex> CX_, final int INCX) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final CX = CX_.having();
  int I, NINCX;
  double STEMP;

  STEMP = 0.0;
  if (N <= 0) return 0;
  if (INCX != 1) {
    // CODE FOR INCREMENT NOT EQUAL TO 1

    NINCX = N * INCX;
    for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
      // NEXT LINE MODIFIED.

      STEMP = STEMP + CX[I].abs();
    }
    return STEMP;

    // CODE FOR INCREMENT EQUAL TO 1
  }
  for (I = 1; I <= N; I++) {
    // NEXT LINE MODIFIED.

    STEMP = STEMP + CX[I].abs();
  }
  return STEMP;
}
