import 'package:lapack/src/complex.dart';
import 'package:lapack/src/blas/dcabs1.dart';
import 'package:lapack/src/matrix.dart';

int izamax(
  final int N,
  final Array<Complex> ZX_,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a softwint are package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();
  double DMAX;
  int I, IX;

  var result = 0;
  if (N < 1 || INCX <= 0) return result;
  result = 1;
  if (N == 1) return result;
  if (INCX == 1) {
    // code for increment equal to 1

    DMAX = dcabs1(ZX[1]);
    for (I = 2; I <= N; I++) {
      if (dcabs1(ZX[I]) > DMAX) {
        result = I;
        DMAX = dcabs1(ZX[I]);
      }
    }
  } else {
    // code for increment not equal to 1

    IX = 1;
    DMAX = dcabs1(ZX[1]);
    IX = IX + INCX;
    for (I = 2; I <= N; I++) {
      if (dcabs1(ZX[IX]) > DMAX) {
        result = I;
        DMAX = dcabs1(ZX[IX]);
      }
      IX = IX + INCX;
    }
  }
  return result;
}
