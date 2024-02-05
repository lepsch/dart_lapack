import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dbdt01(
  final int M,
  final int N,
  final int KD,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> Q,
  final int LDQ,
  final Array<double> D,
  final Array<double> E,
  final Matrix<double> PT,
  final int LDPT,
  final Array<double> WORK,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, EPS;

  // Quick return if possible

  if (M <= 0 || N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Compute A - Q * B * P**T one column at a time.

  RESID.value = ZERO;
  if (KD != 0) {
    // B is bidiagonal.

    if (KD != 0 && M >= N) {
      // B is upper bidiagonal and M >= N.

      for (J = 1; J <= N; J++) {
        dcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= N - 1; I++) {
          WORK[M + I] = D[I] * PT[I][J] + E[I] * PT[I + 1][J];
        }
        WORK[M + N] = D[N] * PT[N][J];
        dgemv('No transpose', M, N, -ONE, Q, LDQ, WORK(M + 1), 1, ONE, WORK, 1);
        RESID.value = max(RESID.value, dasum(M, WORK, 1));
      }
    } else if (KD < 0) {
      // B is upper bidiagonal and M < N.

      for (J = 1; J <= N; J++) {
        dcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= M - 1; I++) {
          WORK[M + I] = D[I] * PT[I][J] + E[I] * PT[I + 1][J];
        }
        WORK[M + M] = D[M] * PT[M][J];
        dgemv('No transpose', M, M, -ONE, Q, LDQ, WORK(M + 1), 1, ONE, WORK, 1);
        RESID.value = max(RESID.value, dasum(M, WORK, 1));
      }
    } else {
      // B is lower bidiagonal.

      for (J = 1; J <= N; J++) {
        dcopy(M, A(1, J).asArray(), 1, WORK, 1);
        WORK[M + 1] = D[1] * PT[1][J];
        for (I = 2; I <= M; I++) {
          WORK[M + I] = E[I - 1] * PT[I - 1][J] + D[I] * PT[I][J];
        }
        dgemv('No transpose', M, M, -ONE, Q, LDQ, WORK(M + 1), 1, ONE, WORK, 1);
        RESID.value = max(RESID.value, dasum(M, WORK, 1));
      }
    }
  } else {
    // B is diagonal.

    if (M >= N) {
      for (J = 1; J <= N; J++) {
        dcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= N; I++) {
          WORK[M + I] = D[I] * PT[I][J];
        }
        dgemv('No transpose', M, N, -ONE, Q, LDQ, WORK(M + 1), 1, ONE, WORK, 1);
        RESID.value = max(RESID.value, dasum(M, WORK, 1));
      }
    } else {
      for (J = 1; J <= N; J++) {
        dcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= M; I++) {
          WORK[M + I] = D[I] * PT[I][J];
        }
        dgemv('No transpose', M, M, -ONE, Q, LDQ, WORK(M + 1), 1, ONE, WORK, 1);
        RESID.value = max(RESID.value, dasum(M, WORK, 1));
      }
    }
  }

  // Compute norm(A - Q * B * P**T) / ( n * norm(A) * EPS )

  ANORM = dlange('1', M, N, A, LDA, WORK);
  EPS = dlamch('Precision');

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
