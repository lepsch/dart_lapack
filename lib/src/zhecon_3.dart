import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrs_3.dart';
import 'package:lapack/src/zlacn2.dart';

void zhecon_3(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> E_,
  final Array<int> IPIV_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final E = E_.having();
  final WORK = WORK_.having();

  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int I;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (ANORM < ZERO) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZHECON_3', -INFO.value);
    return;
  }

  // Quick return if possible

  RCOND.value = ZERO;
  if (N == 0) {
    RCOND.value = ONE;
    return;
  } else if (ANORM <= ZERO) {
    return;
  }

  // Check that the diagonal matrix D is nonsingular.

  if (UPPER) {
    // Upper triangular storage: examine D from bottom to top

    for (I = N; I >= 1; I--) {
      if (IPIV[I] > 0 && A[I][I] == Complex.zero) return;
    }
  } else {
    // Lower triangular storage: examine D from top to bottom.

    for (I = 1; I <= N; I++) {
      if (IPIV[I] > 0 && A[I][I] == Complex.zero) return;
    }
  }

  // Estimate the 1-norm of the inverse.

  KASE.value = 0;
  while (true) {
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    // Multiply by inv(L*D*L**H) or inv(U*D*U**H).

    zhetrs_3(UPLO, N, 1, A, LDA, E, IPIV, WORK.asMatrix(N), N, INFO);
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
