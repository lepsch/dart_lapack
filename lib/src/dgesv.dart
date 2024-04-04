import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgetrf.dart';
import 'package:lapack/src/dgetrs.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgesv(
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DGESV', -INFO.value);
    return;
  }

  // Compute the LU factorization of A.

  dgetrf(N, N, A, LDA, IPIV, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    dgetrs('No transpose', N, NRHS, A, LDA, IPIV, B, LDB, INFO);
  }
}
