import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dbdt05(
  final int M,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Array<double> S,
  final int NS,
  final Matrix<double> U,
  final int LDU,
  final Matrix<double> VT,
  final int LDVT,
  final Array<double> WORK,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, EPS;

  // Quick return if possible.

  RESID.value = ZERO;
  if (min(M, N) <= 0 || NS <= 0) return;

  EPS = dlamch('Precision');
  ANORM = dlange('M', M, N, A, LDA, WORK);

  // Compute U' * A * V.

  dgemm('N', 'T', M, NS, N, ONE, A, LDA, VT, LDVT, ZERO,
      WORK(1 + NS * NS).asMatrix(M), M);
  dgemm('T', 'N', NS, NS, M, -ONE, U, LDU, WORK(1 + NS * NS).asMatrix(M), M,
      ZERO, WORK.asMatrix(NS), NS);

  // norm(S - U' * B * V)

  J = 0;
  for (I = 1; I <= NS; I++) {
    WORK[J + I] = WORK[J + I] + S[I];
    RESID.value = max(RESID.value, dasum(NS, WORK(J + 1), 1));
    J = J + NS;
  }

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (ANORM >= RESID.value) {
      RESID.value = (RESID.value / ANORM) / (N.toDouble() * EPS);
    } else {
      if (ANORM < ONE) {
        RESID.value = (min(RESID.value, (N).toDouble() * ANORM) / ANORM) /
            (N.toDouble() * EPS);
      } else {
        RESID.value =
            min(RESID.value / ANORM, (N).toDouble()) / (N.toDouble() * EPS);
      }
    }
  }
}
