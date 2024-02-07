import 'package:lapack/src/matrix.dart';

void dlar2v(
  final int N,
  final Array<double> X,
  final Array<double> Y,
  final Array<double> Z,
  final int INCX,
  final Array<double> C,
  final Array<double> S,
  final int INCC,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, IC, IX;
  double CI, SI, T1, T2, T3, T4, T5, T6, XI, YI, ZI;

  IX = 1;
  IC = 1;
  for (I = 1; I <= N; I++) {
    XI = X[IX];
    YI = Y[IX];
    ZI = Z[IX];
    CI = C[IC];
    SI = S[IC];
    T1 = SI * ZI;
    T2 = CI * ZI;
    T3 = T2 - SI * XI;
    T4 = T2 + SI * YI;
    T5 = CI * XI + T1;
    T6 = CI * YI - T1;
    X[IX] = CI * T5 + SI * T4;
    Y[IX] = CI * T6 - SI * T3;
    Z[IX] = CI * T4 - SI * T5;
    IX = IX + INCX;
    IC = IC + INCC;
  }
}
