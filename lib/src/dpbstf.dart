import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsyr.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpbstf(
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
  int J, KLD, KM, M;
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
    xerbla('DPBSTF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  KLD = max(1, LDAB - 1);

  // Set the splitting point m.

  M = (N + KD) ~/ 2;

  if (UPPER) {
    // Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).

    for (J = N; J >= M + 1; J--) {
      // Compute s(j,j) and test for non-positive-definiteness.

      AJJ = AB[KD + 1][J];
      if (AJJ <= ZERO) {
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[KD + 1][J] = AJJ;
      KM = min(J - 1, KD);

      // Compute elements j-km:j-1 of the j-th column and update the
      // the leading submatrix within the band.

      dscal(KM, ONE / AJJ, AB(KD + 1 - KM, J).asArray(), 1);
      dsyr('Upper', KM, -ONE, AB(KD + 1 - KM, J).asArray(), 1,
          AB(KD + 1, J - KM), KLD);
    }

    // Factorize the updated submatrix A(1:m,1:m) as U**T*U.

    for (J = 1; J <= M; J++) {
      // Compute s(j,j) and test for non-positive-definiteness.

      AJJ = AB[KD + 1][J];
      if (AJJ <= ZERO) {
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[KD + 1][J] = AJJ;
      KM = min(KD, M - J);

      // Compute elements j+1:j+km of the j-th row and update the
      // trailing submatrix within the band.

      if (KM > 0) {
        dscal(KM, ONE / AJJ, AB(KD, J + 1).asArray(), KLD);
        dsyr('Upper', KM, -ONE, AB(KD, J + 1).asArray(), KLD, AB(KD + 1, J + 1),
            KLD);
      }
    }
  } else {
    // Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).

    for (J = N; J >= M + 1; J--) {
      // Compute s(j,j) and test for non-positive-definiteness.

      AJJ = AB[1][J];
      if (AJJ <= ZERO) {
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[1][J] = AJJ;
      KM = min(J - 1, KD);

      // Compute elements j-km:j-1 of the j-th row and update the
      // trailing submatrix within the band.

      dscal(KM, ONE / AJJ, AB(KM + 1, J - KM).asArray(), KLD);
      dsyr('Lower', KM, -ONE, AB(KM + 1, J - KM).asArray(), KLD, AB(1, J - KM),
          KLD);
    }

    // Factorize the updated submatrix A(1:m,1:m) as U**T*U.

    for (J = 1; J <= M; J++) {
      // Compute s(j,j) and test for non-positive-definiteness.

      AJJ = AB[1][J];
      if (AJJ <= ZERO) {
        INFO.value = J;
        return;
      }
      AJJ = sqrt(AJJ);
      AB[1][J] = AJJ;
      KM = min(KD, M - J);

      // Compute elements j+1:j+km of the j-th column and update the
      // trailing submatrix within the band.

      if (KM > 0) {
        dscal(KM, ONE / AJJ, AB(2, J).asArray(), 1);
        dsyr('Lower', KM, -ONE, AB(2, J).asArray(), 1, AB(1, J + 1), KLD);
      }
    }
  }
}
