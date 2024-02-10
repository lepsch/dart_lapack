import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggsvd3.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dgsvts3(
  final int M,
  final int P,
  final int N,
  final Matrix<double> A,
  final Matrix<double> AF,
  final int LDA,
  final Matrix<double> B,
  final Matrix<double> BF,
  final int LDB,
  final Matrix<double> U,
  final int LDU,
  final Matrix<double> V,
  final int LDV,
  final Matrix<double> Q,
  final int LDQ,
  final Array<double> ALPHA,
  final Array<double> BETA,
  final Matrix<double> R,
  final int LDR,
  final Array<int> IWORK,
  final Array<double> WORK,
  final int LWORK,
  final Array<double> RWORK,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL;
  final INFO = Box(0), K = Box(0), L = Box(0);

  ULP = dlamch('Precision');
  ULPINV = ONE / ULP;
  UNFL = dlamch('Safe minimum');

  // Copy the matrix A to the array AF.

  dlacpy('Full', M, N, A, LDA, AF, LDA);
  dlacpy('Full', P, N, B, LDB, BF, LDB);

  ANORM = max(dlange('1', M, N, A, LDA, RWORK), UNFL);
  BNORM = max(dlange('1', P, N, B, LDB, RWORK), UNFL);

  // Factorize the matrices A and B in the arrays AF and BF.

  dggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU,
      V, LDV, Q, LDQ, WORK, LWORK, IWORK, INFO);

  // Copy R

  for (I = 1; I <= min(K.value + L.value, M); I++) {
    for (J = I; J <= K.value + L.value; J++) {
      R[I][J] = AF[I][N - K.value - L.value + J];
    }
  }

  if (M - K.value - L.value < 0) {
    for (I = M + 1; I <= K.value + L.value; I++) {
      for (J = I; J <= K.value + L.value; J++) {
        R[I][J] = BF[I - K.value][N - K.value - L.value + J];
      }
    }
  }

  // Compute A:= U'*A*Q - D1*R

  dgemm('No transpose', 'No transpose', M, N, N, ONE, A, LDA, Q, LDQ, ZERO,
      WORK.asMatrix(LDA), LDA);

  dgemm('Transpose', 'No transpose', M, N, M, ONE, U, LDU, WORK.asMatrix(LDA),
      LDA, ZERO, A, LDA);

  for (I = 1; I <= K.value; I++) {
    for (J = I; J <= K.value + L.value; J++) {
      A[I][N - K.value - L.value + J] =
          A[I][N - K.value - L.value + J] - R[I][J];
    }
  }

  for (I = K.value + 1; I <= min(K.value + L.value, M); I++) {
    for (J = I; J <= K.value + L.value; J++) {
      A[I][N - K.value - L.value + J] =
          A[I][N - K.value - L.value + J] - ALPHA[I] * R[I][J];
    }
  }

  // Compute norm( U'*A*Q - D1*R ) / ( max(1,M,N)*norm(A)*ULP ) .

  RESID = dlange('1', M, N, A, LDA, RWORK);

  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, max(M, N))) / ANORM) / ULP;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute B := V'*B*Q - D2*R

  dgemm('No transpose', 'No transpose', P, N, N, ONE, B, LDB, Q, LDQ, ZERO,
      WORK.asMatrix(LDB), LDB);

  dgemm('Transpose', 'No transpose', P, N, P, ONE, V, LDV, WORK.asMatrix(LDB),
      LDB, ZERO, B, LDB);

  for (I = 1; I <= L.value; I++) {
    for (J = I; J <= L.value; J++) {
      B[I][N - L.value + J] = B[I][N - L.value + J] -
          BETA[K.value + I] * R[K.value + I][K.value + J];
    }
  }

  // Compute norm( V'*B*Q - D2*R ) / ( max(P,N)*norm(B)*ULP ) .

  RESID = dlange('1', P, N, B, LDB, RWORK);
  if (BNORM > ZERO) {
    RESULT[2] = ((RESID / max(1, max(P, N))) / BNORM) / ULP;
  } else {
    RESULT[2] = ZERO;
  }

  // Compute I - U'*U

  dlaset('Full', M, M, ZERO, ONE, WORK.asMatrix(LDQ), LDQ);
  dsyrk('Upper', 'Transpose', M, M, -ONE, U, LDU, ONE, WORK.asMatrix(LDU), LDU);

  // Compute norm( I - U'*U ) / ( M * ULP ) .

  RESID = dlansy('1', 'Upper', M, WORK, LDU, RWORK);
  RESULT[3] = (RESID / (max(1, M)).toDouble()) / ULP;

  // Compute I - V'*V

  dlaset('Full', P, P, ZERO, ONE, WORK.asMatrix(LDV), LDV);
  dsyrk('Upper', 'Transpose', P, P, -ONE, V, LDV, ONE, WORK.asMatrix(LDV), LDV);

  // Compute norm( I - V'*V ) / ( P * ULP ) .

  RESID = dlansy('1', 'Upper', P, WORK, LDV, RWORK);
  RESULT[4] = (RESID / (max(1, P)).toDouble()) / ULP;

  // Compute I - Q'*Q

  dlaset('Full', N, N, ZERO, ONE, WORK.asMatrix(LDQ), LDQ);
  dsyrk('Upper', 'Transpose', N, N, -ONE, Q, LDQ, ONE, WORK.asMatrix(LDQ), LDQ);

  // Compute norm( I - Q'*Q ) / ( N * ULP ) .

  RESID = dlansy('1', 'Upper', N, WORK, LDQ, RWORK);
  RESULT[5] = (RESID / (max(1, N)).toDouble()) / ULP;

  // Check sorting

  dcopy(N, ALPHA, 1, WORK, 1);
  for (I = K.value + 1; I <= min(K.value + L.value, M); I++) {
    J = IWORK[I];
    if (I != J) {
      TEMP = WORK[I];
      WORK[I] = WORK[J];
      WORK[J] = TEMP;
    }
  }

  RESULT[6] = ZERO;
  for (I = K.value + 1; I <= min(K.value + L.value, M) - 1; I++) {
    if (WORK[I] < WORK[I + 1]) RESULT[6] = ULPINV;
  }
}
