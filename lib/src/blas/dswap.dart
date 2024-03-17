import 'package:lapack/src/matrix.dart';

void dswap(
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

  if (N <= 0) return;
  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    // clean-up loop

    final M = N % 3;
    if (M != 0) {
      for (var I = 1; I <= M; I++) {
        final DTEMP = DX[I];
        DX[I] = DY[I];
        DY[I] = DTEMP;
      }
      if (N < 3) return;
    }
    final MP1 = M + 1;
    for (var I = MP1; I <= N; I += 3) {
      var DTEMP = DX[I];
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

    var IX = INCX < 0 ? (-N + 1) * INCX + 1 : 1;
    var IY = INCY < 0 ? (-N + 1) * INCY + 1 : 1;
    for (var I = 1; I <= N; I++) {
      final DTEMP = DX[IX];
      DX[IX] = DY[IY];
      DY[IY] = DTEMP;
      IX += INCX;
      IY += INCY;
    }
  }
}
