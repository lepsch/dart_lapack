import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zlaipd(
  final int N,
  final Array<Complex> A_,
  final int INDA,
  final int VINDA,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();

  final BIGNUM = dlamch('Epsilon') / dlamch('Safe minimum');
  var IA = 1;
  var IXA = INDA;
  for (var I = 1; I <= N; I++) {
    A[IA] = Complex(A[IA].real, BIGNUM);
    IA = IA + IXA;
    IXA = IXA + VINDA;
  }
}
