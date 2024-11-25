// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacn2.dart';
import 'package:dart_lapack/src/zsytrs.dart';

double zla_syrcond_c(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AF_,
  final int LDAF,
  final Array<int> IPIV_,
  final Array<double> C_,
  final bool CAPPLY,
  final Box<int> INFO,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDAF);
  final C = C_.having();
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  double ANORM, TMP;
  int I, J;
  bool UP, UPPER;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LDAF < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZLA_SYRCOND_C', -INFO.value);
    return 0;
  }
  UP = false;
  if (lsame(UPLO, 'U')) UP = true;

  // Compute norm of op(A)*op2(C).

  ANORM = 0.0;
  if (UP) {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CAPPLY) {
        for (J = 1; J <= I; J++) {
          TMP += A[J][I].cabs1() / C[J];
        }
        for (J = I + 1; J <= N; J++) {
          TMP += A[I][J].cabs1() / C[J];
        }
      } else {
        for (J = 1; J <= I; J++) {
          TMP += A[J][I].cabs1();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += A[I][J].cabs1();
        }
      }
      RWORK[I] = TMP;
      ANORM = max(ANORM, TMP);
    }
  } else {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CAPPLY) {
        for (J = 1; J <= I; J++) {
          TMP += A[I][J].cabs1() / C[J];
        }
        for (J = I + 1; J <= N; J++) {
          TMP += A[J][I].cabs1() / C[J];
        }
      } else {
        for (J = 1; J <= I; J++) {
          TMP += A[I][J].cabs1();
        }
        for (J = I + 1; J <= N; J++) {
          TMP += A[J][I].cabs1();
        }
      }
      RWORK[I] = TMP;
      ANORM = max(ANORM, TMP);
    }
  }

  // Quick return if possible.

  if (N == 0) {
    return 1;
  } else if (ANORM == 0.0) {
    return 0;
  }

  // Estimate the norm of inv(op(A)).

  AINVNM.value = 0.0;

  KASE.value = 0;
  while (true) {
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;
    if (KASE.value == 2) {
      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] *= RWORK[I].toComplex();
      }

      if (UP) {
        zsytrs('U', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      } else {
        zsytrs('L', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      }

      // Multiply by inv(C).

      if (CAPPLY) {
        for (I = 1; I <= N; I++) {
          WORK[I] *= C[I].toComplex();
        }
      }
    } else {
      // Multiply by inv(C**T).

      if (CAPPLY) {
        for (I = 1; I <= N; I++) {
          WORK[I] *= C[I].toComplex();
        }
      }

      if (UP) {
        zsytrs('U', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      } else {
        zsytrs('L', N, 1, AF, LDAF, IPIV, WORK.asMatrix(N), N, INFO);
      }

      // Multiply by R.

      for (I = 1; I <= N; I++) {
        WORK[I] *= RWORK[I].toComplex();
      }
    }
  }

  // Compute the estimate of the reciprocal condition number.

  return AINVNM.value != 0.0 ? 1.0 / AINVNM.value : 0;
}
