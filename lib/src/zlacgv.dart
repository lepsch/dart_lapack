import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlacgv(
  final int N,
  final Array<Complex> X_,
  final int INCX,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  int I, IOFF;
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC DCONJG

  if (INCX == 1) {
    for (I = 1; I <= N; I++) {
      // 10
      X[I] = X[I].conjugate();
    } // 10
  } else {
    IOFF = 1;
    if (INCX < 0) IOFF = 1 - (N - 1) * INCX;
    for (I = 1; I <= N; I++) {
      // 20
      X[IOFF] = X[IOFF].conjugate();
      IOFF = IOFF + INCX;
    } // 20
  }
}
