import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/variants/cholesky/top/zpotrf.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zpotrs.dart';

void zposv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZPOSV ', -INFO.value);
    return;
  }

  // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

  zpotrf(UPLO, N, A, LDA, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    zpotrs(UPLO, N, NRHS, A, LDA, B, LDB, INFO);
  }
}
