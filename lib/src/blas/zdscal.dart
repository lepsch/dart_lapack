import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zdscal(
  final int N,
  final double DA,
  final Array<Complex> ZX,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, NINCX;
  const ONE = 1.0;

  if (N <= 0 || INCX <= 0 || DA == ONE) return;
  if (INCX == 1) {
    // code for increment equal to 1

    for (I = 1; I <= N; I++) {
      ZX[I] = Complex(DA * ZX[I].real, DA * ZX[I].imaginary);
    }
  } else {
    // code for increment not equal to 1

    NINCX = N * INCX;
    for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
      ZX[I] = Complex(DA * ZX[I].real, DA * ZX[I].imaginary);
    }
  }
}
