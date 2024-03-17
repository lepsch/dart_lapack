import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlaev2.dart';

void zlaev2(
  final Complex A,
  final Complex B,
  final Complex C,
  final Box<double> RT1,
  final Box<double> RT2,
  final Box<double> CS1,
  final Box<Complex> SN1,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  final T = Box(0.0);
  Complex W;

  if (B.abs() == ZERO) {
    W = Complex.one;
  } else {
    W = B.conjugate() / B.abs().toComplex();
  }
  dlaev2(A.toDouble(), B.abs(), C.toDouble(), RT1, RT2, CS1, T);
  SN1.value = W * T.value.toComplex();
}
