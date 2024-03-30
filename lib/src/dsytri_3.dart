import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dsytri_3x.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytri_3(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> E_,
  final Array<int> IPIV_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  bool UPPER, LQUERY;
  int LWKOPT, NB = 0;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);

  // Determine the block size

  if (N == 0) {
    LWKOPT = 1;
  } else {
    NB = max(1, ilaenv(1, 'DSYTRI_3', UPLO, N, -1, -1, -1));
    LWKOPT = (N + NB + 1) * (NB + 3);
  }
  WORK[1] = LWKOPT.toDouble();

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < LWKOPT && !LQUERY) {
    INFO.value = -8;
  }

  if (INFO.value != 0) {
    xerbla('DSYTRI_3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  dsytri_3x(UPLO, N, A, LDA, E, IPIV, WORK.asMatrix(NB), NB, INFO);

  WORK[1] = LWKOPT.toDouble();
}
