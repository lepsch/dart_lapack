import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrf_rook.dart';
import 'package:lapack/src/zhetrs_rook.dart';

void zhesv_rook(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
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
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  bool LQUERY;
  int LWKOPT = 0, NB;

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
      NB = ilaenv(1, 'ZHETRF_ROOK', UPLO, N, -1, -1, -1);
      LWKOPT = N * NB;
    }
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZHESV_ROOK ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Compute the factorization A = U*D*U**H or A = L*D*L**H.

  zhetrf_rook(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    // Solve with TRS ( Use Level BLAS 2)

    zhetrs_rook(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
  }

  WORK[1] = LWKOPT.toComplex();
}
