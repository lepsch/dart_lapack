import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dbdt04(
  final String UPLO,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> S_,
  final int NS,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final Array<double> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having(length: N);
  final E = E_.having(length: N - 1);
  final S = S_.having(length: NS);
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final WORK = WORK_.having(length: 2 * N);
  const ZERO = 0.0, ONE = 1.0;
  int I, J, K;
  double BNORM, EPS;

  // Quick return if possible.

  RESID.value = ZERO;
  if (N <= 0 || NS <= 0) return;

  EPS = dlamch('Precision');

  // Compute S - U' * B * V.

  BNORM = ZERO;

  if (lsame(UPLO, 'U')) {
    // B is upper bidiagonal.

    K = 0;
    for (I = 1; I <= NS; I++) {
      for (J = 1; J <= N - 1; J++) {
        K++;
        WORK[K] = D[J] * VT[I][J] + E[J] * VT[I][J + 1];
      }
      K++;
      WORK[K] = D[N] * VT[I][N];
    }
    BNORM = (D[1]).abs();
    for (I = 2; I <= N; I++) {
      BNORM = max(BNORM, (D[I]).abs() + (E[I - 1]).abs());
    }
  } else {
    // B is lower bidiagonal.

    K = 0;
    for (I = 1; I <= NS; I++) {
      K++;
      WORK[K] = D[1] * VT[I][1];
      for (J = 1; J <= N - 1; J++) {
        K++;
        WORK[K] = E[J] * VT[I][J] + D[J + 1] * VT[I][J + 1];
      }
    }
    BNORM = (D[N]).abs();
    for (I = 1; I <= N - 1; I++) {
      BNORM = max(BNORM, (D[I]).abs() + (E[I]).abs());
    }
  }

  dgemm('T', 'N', NS, NS, N, -ONE, U, LDU, WORK(1).asMatrix(N), N, ZERO,
      WORK(1 + N * NS).asMatrix(NS), NS);

  // norm(S - U' * B * V)

  K = N * NS;
  for (I = 1; I <= NS; I++) {
    WORK[K + I] += S[I];
    RESID.value = max(RESID.value, dasum(NS, WORK(K + 1), 1));
    K += NS;
  }

  if (BNORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (BNORM >= RESID.value) {
      RESID.value = (RESID.value / BNORM) / (N.toDouble() * EPS);
    } else {
      if (BNORM < ONE) {
        RESID.value = (min(RESID.value, (N).toDouble() * BNORM) / BNORM) /
            (N.toDouble() * EPS);
      } else {
        RESID.value =
            min(RESID.value / BNORM, (N).toDouble()) / (N.toDouble() * EPS);
      }
    }
  }
}
