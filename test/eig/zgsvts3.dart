// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zggsvd3.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zgsvts3(
  final int M,
  final int P,
  final int N,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> BF_,
  final int LDB,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<double> ALPHA_,
  final Array<double> BETA_,
  final Matrix<Complex> R_,
  final int LDR,
  final Array<int> IWORK_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final BF = BF_.having(ld: LDB);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final Q = Q_.having(ld: LDQ);
  final R = R_.having(ld: LDR);
  final IWORK = IWORK_.having();
  final WORK = WORK_.having(length: LWORK);
  final RWORK = RWORK_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL;
  final INFO = Box(0), K = Box(0), L = Box(0);

  ULP = dlamch('Precision');
  ULPINV = ONE / ULP;
  UNFL = dlamch('Safe minimum');

  // Copy the matrix A to the array AF.

  zlacpy('Full', M, N, A, LDA, AF, LDA);
  zlacpy('Full', P, N, B, LDB, BF, LDB);

  ANORM = max(zlange('1', M, N, A, LDA, RWORK), UNFL);
  BNORM = max(zlange('1', P, N, B, LDB, RWORK), UNFL);

  // Factorize the matrices A and B in the arrays AF and BF.

  zggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU,
      V, LDV, Q, LDQ, WORK, LWORK, RWORK, IWORK, INFO);

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

  zgemm('No transpose', 'No transpose', M, N, N, Complex.one, A, LDA, Q, LDQ,
      Complex.zero, WORK.asMatrix(), LDA);

  zgemm('Conjugate transpose', 'No transpose', M, N, M, Complex.one, U, LDU,
      WORK.asMatrix(), LDA, Complex.zero, A, LDA);

  for (I = 1; I <= K.value; I++) {
    for (J = I; J <= K.value + L.value; J++) {
      A[I][N - K.value - L.value + J] =
          A[I][N - K.value - L.value + J] - R[I][J];
    }
  }

  for (I = K.value + 1; I <= min(K.value + L.value, M); I++) {
    for (J = I; J <= K.value + L.value; J++) {
      A[I][N - K.value - L.value + J] =
          A[I][N - K.value - L.value + J] - ALPHA[I].toComplex() * R[I][J];
    }
  }

  // Compute norm( U'*A*Q - D1*R ) / ( max(1,M,N)*norm(A)*ULP ) .

  RESID = zlange('1', M, N, A, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, max(M, N))) / ANORM) / ULP;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute B := V'*B*Q - D2*R

  zgemm('No transpose', 'No transpose', P, N, N, Complex.one, B, LDB, Q, LDQ,
      Complex.zero, WORK.asMatrix(), LDB);

  zgemm('Conjugate transpose', 'No transpose', P, N, P, Complex.one, V, LDV,
      WORK.asMatrix(), LDB, Complex.zero, B, LDB);

  for (I = 1; I <= L.value; I++) {
    for (J = I; J <= L.value; J++) {
      B[I][N - L.value + J] -=
          BETA[K.value + I].toComplex() * R[K.value + I][K.value + J];
    }
  }

  // Compute norm( V'*B*Q - D2*R ) / ( max(P,N)*norm(B)*ULP ) .

  RESID = zlange('1', P, N, B, LDB, RWORK);
  if (BNORM > ZERO) {
    RESULT[2] = ((RESID / max(1, max(P, N))) / BNORM) / ULP;
  } else {
    RESULT[2] = ZERO;
  }

  // Compute I - U'*U

  zlaset('Full', M, M, Complex.zero, Complex.one, WORK.asMatrix(), LDQ);
  zherk('Upper', 'Conjugate transpose', M, M, -ONE, U, LDU, ONE,
      WORK.asMatrix(), LDU);

  // Compute norm( I - U'*U ) / ( M * ULP ) .

  RESID = zlanhe('1', 'Upper', M, WORK.asMatrix(), LDU, RWORK);
  RESULT[3] = (RESID / max(1, M)) / ULP;

  // Compute I - V'*V

  zlaset('Full', P, P, Complex.zero, Complex.one, WORK.asMatrix(), LDV);
  zherk('Upper', 'Conjugate transpose', P, P, -ONE, V, LDV, ONE,
      WORK.asMatrix(), LDV);

  // Compute norm( I - V'*V ) / ( P * ULP ) .

  RESID = zlanhe('1', 'Upper', P, WORK.asMatrix(), LDV, RWORK);
  RESULT[4] = (RESID / max(1, P)) / ULP;

  // Compute I - Q'*Q

  zlaset('Full', N, N, Complex.zero, Complex.one, WORK.asMatrix(), LDQ);
  zherk('Upper', 'Conjugate transpose', N, N, -ONE, Q, LDQ, ONE,
      WORK.asMatrix(), LDQ);

  // Compute norm( I - Q'*Q ) / ( N * ULP ) .

  RESID = zlanhe('1', 'Upper', N, WORK.asMatrix(), LDQ, RWORK);
  RESULT[5] = (RESID / max(1, N)) / ULP;

  // Check sorting

  dcopy(N, ALPHA, 1, RWORK, 1);
  for (I = K.value + 1; I <= min(K.value + L.value, M); I++) {
    J = IWORK[I];
    if (I != J) {
      TEMP = RWORK[I];
      RWORK[I] = RWORK[J];
      RWORK[J] = TEMP;
    }
  }

  RESULT[6] = ZERO;
  for (I = K.value + 1; I <= min(K.value + L.value, M) - 1; I++) {
    if (RWORK[I] < RWORK[I + 1]) RESULT[6] = ULPINV;
  }
}
