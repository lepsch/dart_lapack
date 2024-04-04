import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zsptrf.dart';
import 'package:lapack/src/zsptrs.dart';

void zspsv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZSPSV', -INFO.value);
    return;
  }

  // Compute the factorization A = U*D*U**T or A = L*D*L**T.

  zsptrf(UPLO, N, AP, IPIV, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    zsptrs(UPLO, N, NRHS, AP, IPIV, B, LDB, INFO);
  }
}
