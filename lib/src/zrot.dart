import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zrot(
  final int N,
  final Array<Complex> CX,
  final int INCX,
  final Array<Complex> CY,
  final int INCY,
  final double C,
  final Complex S,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, IX, IY;
  Complex STEMP;

  if (N <= 0) return;
  if (INCX != 1 || INCY != 1) {
    // Code for unequal increments or equal increments not equal to 1

    IX = 1;
    IY = 1;
    if (INCX < 0) IX = (-N + 1) * INCX + 1;
    if (INCY < 0) IY = (-N + 1) * INCY + 1;
    for (I = 1; I <= N; I++) {
      // 10
      STEMP = C.toComplex() * CX[IX] + S * CY[IY];
      CY[IY] = C.toComplex() * CY[IY] - S.conjugate() * CX[IX];
      CX[IX] = STEMP;
      IX = IX + INCX;
      IY = IY + INCY;
    } // 10
    return;
  } // 20

  // Code for both increments equal to 1

  for (I = 1; I <= N; I++) {
    // 30
    STEMP = C.toComplex() * CX[I] + S * CY[I];
    CY[I] = C.toComplex() * CY[I] - S.conjugate() * CX[I];
    CX[I] = STEMP;
  } // 30
}
