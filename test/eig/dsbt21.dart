import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dspr.dart';
import 'package:lapack/src/blas/dspr2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansb.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dsbt21(
  final String UPLO,
  final int N,
  final int KA,
  final int KS,
  final Matrix<double> A,
  final int LDA,
  final Array<double> D,
  final Array<double> E,
  final Matrix<double> U,
  final int LDU,
  final Array<double> WORK,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  bool LOWER;
  String CUPLO;
  int IKA, J, JC, JR, LW;
  double ANORM, ULP, UNFL, WNORM;

  // Constants

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  IKA = max(0, min(N - 1, KA));
  LW = (N * (N + 1)) ~/ 2;

  if (lsame(UPLO, 'U')) {
    LOWER = false;
    CUPLO = 'U';
  } else {
    LOWER = true;
    CUPLO = 'L';
  }

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');

  // Some Error Checks

  // Do Test 1

  // Norm of A:

  ANORM = max(dlansb('1', CUPLO, N, IKA, A, LDA, WORK), UNFL);

  // Compute error matrix:    Error = A - U S U**T

  // Copy A from SB to SP storage format.

  J = 0;
  for (JC = 1; JC <= N; JC++) {
    if (LOWER) {
      for (JR = 1; JR <= min(IKA + 1, N + 1 - JC); JR++) {
        J = J + 1;
        WORK[J] = A[JR][JC];
      }
      for (JR = IKA + 2; JR <= N + 1 - JC; JR++) {
        J = J + 1;
        WORK[J] = ZERO;
      }
    } else {
      for (JR = IKA + 2; JR <= JC; JR++) {
        J = J + 1;
        WORK[J] = ZERO;
      }
      for (JR = min(IKA, JC - 1); JR >= 0; JR--) {
        J = J + 1;
        WORK[J] = A[IKA + 1 - JR][JC];
      }
    }
  }

  for (J = 1; J <= N; J++) {
    dspr(CUPLO, N, -D[J], U(1, J).asArray(), 1, WORK);
  }

  if (N > 1 && KS == 1) {
    for (J = 1; J <= N - 1; J++) {
      dspr2(
        CUPLO,
        N,
        -E[J],
        U(1, J).asArray(),
        1,
        U(1, J + 1).asArray(),
        1,
        WORK,
      );
    }
  }
  WNORM = dlansp('1', CUPLO, N, WORK, WORK[LW + 1]);

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (N * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, N.toDouble()) / (N * ULP);
    }
  }

  // Do Test 2

  // Compute  U U**T - I

  dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK.asMatrix(N), N);

  for (J = 1; J <= N; J++) {
    WORK[(N + 1) * (J - 1) + 1] = WORK[(N + 1) * (J - 1) + 1] - ONE;
  }

  RESULT[2] = min(
        dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1)),
        N.toDouble(),
      ) /
      (N * ULP);
}
