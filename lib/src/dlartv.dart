import 'package:lapack/src/matrix.dart';

void dlartv(
  final int N,
  final Array<double> X,
  final int INCX,
  final Array<double> Y,
  final int INCY,
  final Array<double> C,
  final Array<double> S,
  final int INCC,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, IC, IX, IY;
  double XI, YI;

  IX = 1;
  IY = 1;
  IC = 1;
  for (I = 1; I <= N; I++) {
    XI = X[IX];
    YI = Y[IY];
    X[IX] = C[IC] * XI + S[IC] * YI;
    Y[IY] = C[IC] * YI - S[IC] * XI;
    IX = IX + INCX;
    IY = IY + INCY;
    IC = IC + INCC;
  }
}
