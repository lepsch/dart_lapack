// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dpotrs.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

double dla_porcond(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AF_,
  final int LDAF,
  final int CMODE,
  final Array<double> C_,
  final Box<int> INFO,
  final Array<double> WORK_,
  final Array<int> IWORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final C = C_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  int I, J;
  double TMP;
  bool UP;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('DLA_PORCOND', -INFO.value);
    return 0;
  }

  if (N == 0) {
    return 1;
  }
  UP = false;
  if (lsame(UPLO, 'U')) UP = true;

  // Compute the equilibration matrix R such that
  // inv(R)*A*C has unit 1-norm.

  if (UP) {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CMODE == 1) {
        for (J = 1; J <= I; J++) {
          TMP += (A[J][I] * C[J]).abs();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += (A[I][J] * C[J]).abs();
        }
      } else if (CMODE == 0) {
        for (J = 1; J <= I; J++) {
          TMP += A[J][I].abs();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += A[I][J].abs();
        }
      } else {
        for (J = 1; J <= I; J++) {
          TMP += (A[J][I] / C[J]).abs();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += (A[I][J] / C[J]).abs();
        }
      }
      WORK[2 * N + I] = TMP;
    }
  } else {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CMODE == 1) {
        for (J = 1; J <= I; J++) {
          TMP += (A[I][J] * C[J]).abs();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += (A[J][I] * C[J]).abs();
        }
      } else if (CMODE == 0) {
        for (J = 1; J <= I; J++) {
          TMP += A[I][J].abs();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += A[J][I].abs();
        }
      } else {
        for (J = 1; J <= I; J++) {
          TMP += (A[I][J] / C[J]).abs();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += (A[J][I] / C[J]).abs();
        }
      }
      WORK[2 * N + I] = TMP;
    }
  }

  // Estimate the norm of inv(op(A)).

  AINVNM.value = 0.0;

  KASE.value = 0;
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;
    if (KASE.value == 2) {
      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] *= WORK[2 * N + I];
      }

      if (UP) {
        dpotrs('Upper', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
      } else {
        dpotrs('Lower', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
      }

      // Multiply by inv(C).

      if (CMODE == 1) {
        for (I = 1; I <= N; I++) {
          WORK[I] /= C[I];
        }
      } else if (CMODE == -1) {
        for (I = 1; I <= N; I++) {
          WORK[I] *= C[I];
        }
      }
    } else {
      // Multiply by inv(C**T).

      if (CMODE == 1) {
        for (I = 1; I <= N; I++) {
          WORK[I] /= C[I];
        }
      } else if (CMODE == -1) {
        for (I = 1; I <= N; I++) {
          WORK[I] *= C[I];
        }
      }

      if (UP) {
        dpotrs('Upper', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
      } else {
        dpotrs('Lower', N, 1, AF, LDAF, WORK.asMatrix(N), N, INFO);
      }

      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] *= WORK[2 * N + I];
      }
    }
  }

  // Compute the estimate of the reciprocal condition number.

  return AINVNM.value != 0.0 ? 1.0 / AINVNM.value : 0;
}
