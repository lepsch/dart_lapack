import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dsytrf.dart';
import 'package:lapack/src/dsytrs.dart';
import 'package:lapack/src/dsytrs2.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsysv(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  bool LQUERY;
  int LWKOPT = 0;

  // Test the input parameters.

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -10;
  }

  if (INFO.value == 0) {
    if (N == 0) {
      LWKOPT = 1;
    } else {
      dsytrf(UPLO, N, A, LDA, IPIV, WORK, -1, INFO);
      LWKOPT = WORK[1].toInt();
    }
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DSYSV ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Compute the factorization A = U*D*U**T or A = L*D*L**T.

  dsytrf(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    if (LWORK < N) {
      // Solve with TRS ( Use Level BLAS 2)

      dsytrs(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
    } else {
      // Solve with TRS2 ( Use Level BLAS 3)

      dsytrs2(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, INFO);
    }
  }

  WORK[1] = LWKOPT.toDouble();
}
