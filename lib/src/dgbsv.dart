import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbtrf.dart';
import 'package:lapack/src/dgbtrs.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgbsv(
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (KL < 0) {
    INFO.value = -2;
  } else if (KU < 0) {
    INFO.value = -3;
  } else if (NRHS < 0) {
    INFO.value = -4;
  } else if (LDAB < 2 * KL + KU + 1) {
    INFO.value = -6;
  } else if (LDB < max(N, 1)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DGBSV ', -INFO.value);
    return;
  }

  // Compute the LU factorization of the band matrix A.

  dgbtrf(N, N, KL, KU, AB, LDAB, IPIV, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    dgbtrs('No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO);
  }
}
