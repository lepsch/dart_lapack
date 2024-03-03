import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgetrs.dart';
import 'package:lapack/src/matrix.dart';

void sgetrs(
  final String TRANS,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  dgetrs(TRANS, N, NRHS, A_, LDA, IPIV_, B_, LDB, INFO);
}
