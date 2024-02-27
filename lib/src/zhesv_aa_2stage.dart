import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrf_aa_2stage.dart';
import 'package:lapack/src/zhetrs_aa_2stage.dart';

void zhesv_aa_2stage(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TB_,
  final int LTB,
  final Array<int> IPIV_,
  final Array<int> IPIV2_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final IPIV = IPIV_.dim();
  final IPIV2 = IPIV2_.dim();
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
  final TB = TB_.dim();
  bool UPPER, TQUERY, WQUERY;
  int LWKOPT = 0, LWKMIN;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  WQUERY = (LWORK == -1);
  TQUERY = (LTB == -1);
  LWKMIN = max(1, N);
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LTB < max(1, 4 * N) && !TQUERY) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -11;
  } else if (LWORK < LWKMIN && !WQUERY) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    zhetrf_aa_2stage(UPLO, N, A, LDA, TB, -1, IPIV, IPIV2, WORK, -1, INFO);
    LWKOPT = max(LWKMIN, WORK[1].toInt());
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZHESV_AA_2STAGE', -INFO.value);
    return;
  } else if (WQUERY || TQUERY) {
    return;
  }

  // Compute the factorization A = U**H*T*U or A = L*T*L**H.

  zhetrf_aa_2stage(UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO);
  if (INFO.value == 0) {
    // Solve the system A*X = B, overwriting B with X.

    zhetrs_aa_2stage(UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO);
  }

  WORK[1] = LWKOPT.toComplex();
}
