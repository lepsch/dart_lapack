import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zcopy(
  final int N,
  final Array<Complex> ZX,
  final int INCX,
  final Array<Complex> ZY,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, IX, IY;

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    for (I = 1; I <= N; I++) {
      ZY[I] = ZX[I];
    }
  } else {
    // code for unequal increments or equal increments
    // not equal to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      ZY[IY] = ZX[IX];
      IX = IX + INCX;
      IY = IY + INCY;
    }
  }
}
