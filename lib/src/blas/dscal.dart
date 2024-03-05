import 'package:lapack/src/matrix.dart';

void dscal(
  final int N,
  final double DA,
  final Array<double> DX_,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DX = DX_.having();
  int I, M, MP1, NINCX;
  const ONE = 1.0;

  if (N <= 0 || INCX <= 0 || DA == ONE) return;

  if (INCX == 1) {
    // code for increment equal to 1

    // clean-up loop

    M = (N % 5);
    if (M != 0) {
      for (I = 1; I <= M; I++) {
        DX[I] = DA * DX[I];
      }
      if (N < 5) return;
    }
    MP1 = M + 1;
    for (I = MP1; I <= N; I += 5) {
      DX[I] = DA * DX[I];
      DX[I + 1] = DA * DX[I + 1];
      DX[I + 2] = DA * DX[I + 2];
      DX[I + 3] = DA * DX[I + 3];
      DX[I + 4] = DA * DX[I + 4];
    }
  } else {
    // code for increment not equal to 1

    NINCX = N * INCX;
    for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
      DX[I] = DA * DX[I];
    }
  }
}
