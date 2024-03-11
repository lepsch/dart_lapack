import 'package:lapack/src/matrix.dart';

double dasum(final int N, Array<double> DX, final int INCX) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  double DTEMP;
  int I, M, MP1, NINCX;
  DTEMP = 0.0;

  if (N <= 0 || INCX <= 0) return 0.0;
  if (INCX == 1) {
    // code for increment equal to 1

    // clean-up loop

    M = (N % 6);
    if (M != 0) {
      for (I = 1; I <= M; I++) {
        DTEMP += DX[I].abs();
      }
      if (N < 6) {
        return DTEMP;
      }
    }
    MP1 = M + 1;
    for (I = MP1; I <= N; I += 6) {
      DTEMP = DTEMP +
          DX[I].abs() +
          DX[I + 1].abs() +
          DX[I + 2].abs() +
          DX[I + 3].abs() +
          DX[I + 4].abs() +
          DX[I + 5].abs();
    }
  } else {
    // code for increment not equal to 1

    NINCX = N * INCX;
    for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
      DTEMP += DX[I].abs();
    }
  }
  return DTEMP;
}
