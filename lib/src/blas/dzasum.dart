import 'package:lapack/src/complex.dart';
import 'package:lapack/src/blas/dcabs1.dart';
import 'package:lapack/src/matrix.dart';

double dzasum(
  final int N,
  final Array<Complex> ZX_,
  final int INCX,
) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ZX = ZX_.dim();
  double STEMP;
  int I, NINCX;
  STEMP = 0.0;
  if (N <= 0 || INCX <= 0) return 0.0;
  if (INCX == 1) {
    // code for increment equal to 1

    for (I = 1; I <= N; I++) {
      STEMP = STEMP + dcabs1(ZX[I]);
    }
  } else {
    // code for increment not equal to 1

    NINCX = N * INCX;
    for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) {
      STEMP = STEMP + dcabs1(ZX[I]);
    }
  }
  return STEMP;
}
