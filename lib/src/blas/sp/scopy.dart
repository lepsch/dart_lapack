import 'package:lapack/src/matrix.dart';

void scopy(
  final int N,
  final Array<double> SX,
  final int INCX,
  final Array<double> SY,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, IX, IY, M, MP1;

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    // clean-up loop

    M = (N % 7);
    if (M != 0) {
      for (I = 1; I <= M; I++) {
        SY[I] = SX[I];
      }
      if (N < 7) return;
    }
    MP1 = M + 1;
    for (I = MP1; I <= N; I += 7) {
      SY[I] = SX[I];
      SY[I + 1] = SX[I + 1];
      SY[I + 2] = SX[I + 2];
      SY[I + 3] = SX[I + 3];
      SY[I + 4] = SX[I + 4];
      SY[I + 5] = SX[I + 5];
      SY[I + 6] = SX[I + 6];
    }
  } else {
    // code for unequal increments or equal increments
    // not equal to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      SY[IY] = SX[IX];
      IX = IX + INCX;
      IY = IY + INCY;
    }
  }
}
