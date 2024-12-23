// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dzasum.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlange.dart';

void zbdt01(
  final int M,
  final int N,
  final int KD,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> PT_,
  final int LDPT,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final Q = Q_.having(ld: LDQ);
  final D = D_.having();
  final E = E_.having();
  final PT = PT_.having(ld: LDPT);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, EPS;

  // Quick return if possible

  if (M <= 0 || N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Compute A - Q * B * P**H one column at a time.

  RESID.value = ZERO;
  if (KD != 0) {
    // B is bidiagonal.

    if (KD != 0 && M >= N) {
      // B is upper bidiagonal and M >= N.

      for (J = 1; J <= N; J++) {
        zcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= N - 1; I++) {
          WORK[M + I] =
              D[I].toComplex() * PT[I][J] + E[I].toComplex() * PT[I + 1][J];
        }
        WORK[M + N] = D[N].toComplex() * PT[N][J];
        zgemv('No transpose', M, N, -Complex.one, Q, LDQ, WORK(M + 1), 1,
            Complex.one, WORK, 1);
        RESID.value = max(RESID.value, dzasum(M, WORK, 1));
      }
    } else if (KD < 0) {
      // B is upper bidiagonal and M < N.

      for (J = 1; J <= N; J++) {
        zcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= M - 1; I++) {
          WORK[M + I] =
              D[I].toComplex() * PT[I][J] + E[I].toComplex() * PT[I + 1][J];
        }
        WORK[M + M] = D[M].toComplex() * PT[M][J];
        zgemv('No transpose', M, M, -Complex.one, Q, LDQ, WORK(M + 1), 1,
            Complex.one, WORK, 1);
        RESID.value = max(RESID.value, dzasum(M, WORK, 1));
      }
    } else {
      // B is lower bidiagonal.

      for (J = 1; J <= N; J++) {
        zcopy(M, A(1, J).asArray(), 1, WORK, 1);
        WORK[M + 1] = D[1].toComplex() * PT[1][J];
        for (I = 2; I <= M; I++) {
          WORK[M + I] =
              E[I - 1].toComplex() * PT[I - 1][J] + D[I].toComplex() * PT[I][J];
        }
        zgemv('No transpose', M, M, -Complex.one, Q, LDQ, WORK(M + 1), 1,
            Complex.one, WORK, 1);
        RESID.value = max(RESID.value, dzasum(M, WORK, 1));
      }
    }
  } else {
    // B is diagonal.

    if (M >= N) {
      for (J = 1; J <= N; J++) {
        zcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= N; I++) {
          WORK[M + I] = D[I].toComplex() * PT[I][J];
        }
        zgemv('No transpose', M, N, -Complex.one, Q, LDQ, WORK(M + 1), 1,
            Complex.one, WORK, 1);
        RESID.value = max(RESID.value, dzasum(M, WORK, 1));
      }
    } else {
      for (J = 1; J <= N; J++) {
        zcopy(M, A(1, J).asArray(), 1, WORK, 1);
        for (I = 1; I <= M; I++) {
          WORK[M + I] = D[I].toComplex() * PT[I][J];
        }
        zgemv('No transpose', M, M, -Complex.one, Q, LDQ, WORK(M + 1), 1,
            Complex.one, WORK, 1);
        RESID.value = max(RESID.value, dzasum(M, WORK, 1));
      }
    }
  }

  // Compute norm(A - Q * B * P**H) / ( n * norm(A) * EPS )

  ANORM = zlange('1', M, N, A, LDA, RWORK);
  EPS = dlamch('Precision');

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (ANORM >= RESID.value) {
      RESID.value = (RESID.value / ANORM) / (N * EPS);
    } else {
      if (ANORM < ONE) {
        RESID.value = (min(RESID.value, N * ANORM) / ANORM) / (N * EPS);
      } else {
        RESID.value = min(RESID.value / ANORM, N) / (N * EPS);
      }
    }
  }
}
