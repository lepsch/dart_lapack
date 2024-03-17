import 'package:lapack/src/complex.dart';
import 'package:lapack/src/blas/dcabs1.dart';
import 'package:lapack/src/matrix.dart';

void zaxpy(
  final int N,
  final Complex ZA,
  final Array<Complex> ZX_,
  final int INCX,
  final Array<Complex> ZY_,
  final int INCY,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.having();
  final ZY = ZY_.having();

  if (N <= 0) return;

  if (dcabs1(ZA) == 0.0) return;

  if (INCX == 1 && INCY == 1) {
    // code for both increments equal to 1

    for (var I = 1; I <= N; I++) {
      ZY[I] += ZA * ZX[I];
    }
  } else {
    // code for unequal increments or equal increments
    // not equal to 1

    var IX = INCX < 0 ? (-N + 1) * INCX + 1 : 1;
    var IY = INCY < 0 ? (-N + 1) * INCY + 1 : 1;
    for (var I = 1; I <= N; I++) {
      ZY[IY] += ZA * ZX[IX];
      IX += INCX;
      IY += INCY;
    }
  }
}
