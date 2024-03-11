import 'package:lapack/src/matrix.dart';

void drot(
  final int N,
  final Array<double> DX_,
  final int INCX,
  final Array<double> DY_,
  final int INCY,
  final double C,
  final double S,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DX = DX_.having();
  final DY = DY_.having();
  double DTEMP;
  int I, IX, IY;

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    for (I = 1; I <= N; I++) {
      DTEMP = C * DX[I] + S * DY[I];
      DY[I] = C * DY[I] - S * DX[I];
      DX[I] = DTEMP;
    }
  } else {
    // code for unequal increments or equal increments not equal
    // to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      DTEMP = C * DX[IX] + S * DY[IY];
      DY[IY] = C * DY[IY] - S * DX[IX];
      DX[IX] = DTEMP;
      IX += INCX;
      IY += INCY;
    }
  }
}
