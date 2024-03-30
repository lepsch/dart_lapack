import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsyr.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpbtf2(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<double> AB_,
  final int LDAB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int J, KLD, KN;
  double AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (LDAB < KD + 1) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DPBTF2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  KLD = max(1, LDAB - 1);

  if (UPPER) {
    // Compute the Cholesky factorization A = U**T*U.

    for (J = 1; J <= N; J++) {
      // Compute U(J,J) and test for non-positive-definiteness.

      AJJ = AB[KD + 1][J];
      if (AJJ <= ZERO) {
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[KD + 1][J] = AJJ;

      // Compute elements J+1:J+KN of row J and update the
      // trailing submatrix within the band.

      KN = min(KD, N - J);
      if (KN > 0) {
        dscal(KN, ONE / AJJ, AB(KD, J + 1).asArray(), KLD);
        dsyr('Upper', KN, -ONE, AB(KD, J + 1).asArray(), KLD, AB(KD + 1, J + 1),
            KLD);
      }
    }
  } else {
    // Compute the Cholesky factorization A = L*L**T.

    for (J = 1; J <= N; J++) {
      // Compute L(J,J) and test for non-positive-definiteness.

      AJJ = AB[1][J];
      if (AJJ <= ZERO) {
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[1][J] = AJJ;

      // Compute elements J+1:J+KN of column J and update the
      // trailing submatrix within the band.

      KN = min(KD, N - J);
      if (KN > 0) {
        dscal(KN, ONE / AJJ, AB(2, J).asArray(), 1);
        dsyr('Lower', KN, -ONE, AB(2, J).asArray(), 1, AB(1, J + 1), KLD);
      }
    }
  }
}
