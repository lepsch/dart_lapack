import 'package:lapack/src/matrix.dart';

double ddot(
  final int N,
  final Array<double> DX_,
  final int INCX,
  final Array<double> DY_,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DX = DX_.having();
  final DY = DY_.having();

  double DTEMP = 0.0;
  if (N <= 0) return 0.0;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    // clean-up loop

    final M = (N % 5);
    if (M != 0) {
      for (var I = 1; I <= M; I++) {
        DTEMP += DX[I] * DY[I];
      }
      if (N < 5) {
        return DTEMP;
      }
    }
    final MP1 = M + 1;
    for (var I = MP1; I <= N; I += 5) {
      DTEMP += DX[I] * DY[I] +
          DX[I + 1] * DY[I + 1] +
          DX[I + 2] * DY[I + 2] +
          DX[I + 3] * DY[I + 3] +
          DX[I + 4] * DY[I + 4];
    }
  } else {
    // code for unequal increments or equal increments
    // not equal to 1

    var IX = 1;
    var IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (var I = 1; I <= N; I++) {
      DTEMP += DX[IX] * DY[IY];
      IX += INCX;
      IY += INCY;
    }
  }
  return DTEMP;
}
