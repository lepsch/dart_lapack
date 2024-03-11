import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zdrot(
  final int N,
  final Array<Complex> ZX_,
  final int INCX,
  final Array<Complex> ZY_,
  final int INCY,
  final double C,
  final double S,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();
  final ZY = ZY_.having();
  int I, IX, IY;
  Complex CTEMP;

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    for (I = 1; I <= N; I++) {
      CTEMP = C.toComplex() * ZX[I] + S.toComplex() * ZY[I];
      ZY[I] = C.toComplex() * ZY[I] - S.toComplex() * ZX[I];
      ZX[I] = CTEMP;
    }
  } else {
    // code for unequal increments or equal increments not equal
    // to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      CTEMP = C.toComplex() * ZX[IX] + S.toComplex() * ZY[IY];
      ZY[IY] = C.toComplex() * ZY[IY] - S.toComplex() * ZX[IX];
      ZX[IX] = CTEMP;
      IX += INCX;
      IY += INCY;
    }
  }
}
