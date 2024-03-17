import 'package:lapack/src/complex.dart';
import 'package:lapack/src/blas/dcabs1.dart';
import 'package:lapack/src/matrix.dart';

int izamax(final int N, final Array<Complex> ZX_, final int INCX) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a softwint are package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();

  if (N < 1 || INCX <= 0) return 0;

  if (N == 1) return 1;

  var index = 1;
  if (INCX == 1) {
    // code for increment equal to 1

    var DMAX = dcabs1(ZX[1]);
    for (var I = 2; I <= N; I++) {
      if (dcabs1(ZX[I]) > DMAX) {
        index = I;
        DMAX = dcabs1(ZX[I]);
      }
    }
  } else {
    // code for increment not equal to 1

    var IX = 1;
    var DMAX = dcabs1(ZX[1]);
    IX += INCX;
    for (var I = 2; I <= N; I++) {
      if (dcabs1(ZX[IX]) > DMAX) {
        index = I;
        DMAX = dcabs1(ZX[IX]);
      }
      IX += INCX;
    }
  }
  return index;
}
