import 'package:lapack/blas.dart';

double dcabs1(final Complex Z) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  return Z.real.abs() + Z.imaginary.abs();
}
