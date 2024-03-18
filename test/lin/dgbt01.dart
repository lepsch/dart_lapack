import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dgbt01(
  final int M,
  final int N,
  final int KL,
  final int KU,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AFAC_,
  final int LDAFAC,
  final Array<int> IPIV_,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;

  // Quick exit if M = 0 or N = 0.

  RESID.value = ZERO;
  if (M <= 0 || N <= 0) return;

  // Determine EPS and the norm of A.

  final EPS = dlamch('Epsilon');
  var KD = KU + 1;
  var ANORM = ZERO;
  for (var J = 1; J <= N; J++) {
    final I1 = max(KD + 1 - J, 1);
    final I2 = min(KD + M - J, KL + KD);
    if (I2 >= I1) ANORM = max(ANORM, dasum(I2 - I1 + 1, A(I1, J).asArray(), 1));
  }

  // Compute one column at a time of L*U - A.

  KD = KL + KU + 1;
  for (var J = 1; J <= N; J++) {
    // Copy the J-th column of U to WORK.

    final JU = min(KL + KU, J - 1);
    final JL = min(KL, M - J);
    final LENJ = min(M, J) - J + JU + 1;
    if (LENJ > 0) {
      dcopy(LENJ, AFAC(KD - JU, J).asArray(), 1, WORK, 1);
      for (var I = LENJ + 1; I <= JU + JL + 1; I++) {
        WORK[I] = ZERO;
      }

      // Multiply by the unit lower triangular matrix L.  Note that L
      // is stored as a product of transformations and permutations.

      for (var I = min(M - 1, J); I >= J - JU; I--) {
        final IL = min(KL, M - I);
        if (IL > 0) {
          final IW = I - J + JU + 1;
          final T = WORK[IW];
          daxpy(IL, T, AFAC(KD + 1, I).asArray(), 1, WORK(IW + 1), 1);
          var IP = IPIV[I];
          if (I != IP) {
            IP -= J - JU - 1;
            WORK[IW] = WORK[IP];
            WORK[IP] = T;
          }
        }
      }

      // Subtract the corresponding column of A.

      final JUA = min(JU, KU);
      if (JUA + JL + 1 > 0) {
        daxpy(JUA + JL + 1, -ONE, A(KU + 1 - JUA, J).asArray(), 1,
            WORK(JU + 1 - JUA), 1);
      }

      // Compute the 1-norm of the column.

      RESID.value = max(RESID.value, dasum(JU + JL + 1, WORK, 1));
    }
  }

  // Compute norm(L*U - A) / ( N * norm(A) * EPS )

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }
}
