import 'package:lapack/src/matrix.dart';

int idamax(final int N, final Array<double> DX_, final int INCX) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DX = DX_.having();
  double DMAX;
  int IX;

  var result = 0;
  if (N < 1 || INCX <= 0) return result;
  result = 1;
  if (N == 1) return result;
  if (INCX == 1) {
    // code for increment equal to 1
    DMAX = DX[1].abs();
    for (var I = 2; I <= N; I++) {
      if (DX[I].abs() > DMAX) {
        result = I;
        DMAX = DX[I].abs();
      }
    }
  } else {
    // code for increment not equal to 1
    IX = 1;
    DMAX = DX[1].abs();
    IX += INCX;
    for (var I = 2; I <= N; I++) {
      if (DX[IX].abs() > DMAX) {
        result = I;
        DMAX = DX[IX].abs();
      }
      IX += INCX;
    }
  }
  return result;
}
