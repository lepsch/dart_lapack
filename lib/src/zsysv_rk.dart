import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zsytrf_rk.dart';
import 'package:lapack/src/zsytrs_3.dart';

void zsysv_rk(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> E_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final E = E_.having();
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
    INFO.value = -9;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -11;
  }

  if (INFO.value == 0) {
    if (N == 0) {
      LWKOPT = 1;
    } else {
      zsytrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, -1, INFO);
      LWKOPT = WORK[1].toInt();
    }
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZSYSV_RK', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Compute the factorization A = P*U*D*(U**T)*(P**T) or
  // A = P*U*D*(U**T)*(P**T).

  zsytrf_rk(UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO);

  if (INFO.value == 0) {
    // Solve the system A*X = B with BLAS3 solver, overwriting B with X.

    zsytrs_3(UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO);
  }

  WORK[1] = LWKOPT.toComplex();
}
