import 'package:lapack/src/matrix.dart';

void dswap(
  final int N,
  final Array<double> DX,
  final int INCX,
  final Array<double> DY,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  double DTEMP;
  int I, IX, IY, M, MP1;

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    // clean-up loop

    M = (N % 3);
    if (M != 0) {
      for (I = 1; I <= M; I++) {
        DTEMP = DX[I];
        DX[I] = DY[I];
        DY[I] = DTEMP;
      }
      if (N < 3) return;
    }
    MP1 = M + 1;
    for (I = MP1; I <= N; I += 3) {
      DTEMP = DX[I];
      DX[I] = DY[I];
      DY[I] = DTEMP;
      DTEMP = DX[I + 1];
      DX[I + 1] = DY[I + 1];
      DY[I + 1] = DTEMP;
      DTEMP = DX[I + 2];
      DX[I + 2] = DY[I + 2];
      DY[I + 2] = DTEMP;
    }
  } else {
    // code for unequal increments or equal increments not equal
    // to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      DTEMP = DX[IX];
      DX[IX] = DY[IY];
      DY[IY] = DTEMP;
      IX = IX + INCX;
      IY = IY + INCY;
    }
  }
}