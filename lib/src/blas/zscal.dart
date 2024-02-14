import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zscal(
  final int N,
  final Complex ZA,
  final Array<Complex> ZX_,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.dim();
  int I, NINCX;

  if (N <= 0 || INCX <= 0 || ZA == Complex.one) return;
  if (INCX == 1) {
    // code for increment equal to 1

    for (I = 1; I <= N; I++) {
      ZX[I] = ZA * ZX[I];
    }
  } else {
    // code for increment not equal to 1

    NINCX = N * INCX;
    for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
      ZX[I] = ZA * ZX[I];
    }
  }
}
