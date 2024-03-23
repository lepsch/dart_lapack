import 'package:lapack/src/dpotrf.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void spotrf(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  dpotrf(UPLO, N, A_, LDA, INFO);
}
