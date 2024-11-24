// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbtrs.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

double dla_gbrcond(
  final String TRANS,
  final int N,
  final int KL,
  final int KU,
  final Matrix<double> AB_,
  final int LDAB,
  final Matrix<double> AFB_,
  final int LDAFB,
  final Array<int> IPIV_,
  final int CMODE,
  final Array<double> C_,
  final Box<int> INFO,
  final Array<double> WORK_,
  final Array<int> IWORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAFB);
  final IPIV = IPIV_.having();
  final C = C_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();

  bool NOTRANS;
  int I, J, KD, KE;
  double TMP;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  INFO.value = 0;
  NOTRANS = lsame(TRANS, 'N');
  if (!NOTRANS && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0 || KL > N - 1) {
    INFO.value = -3;
  } else if (KU < 0 || KU > N - 1) {
    INFO.value = -4;
  } else if (LDAB < KL + KU + 1) {
    INFO.value = -6;
  } else if (LDAFB < 2 * KL + KU + 1) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('DLA_GBRCOND', -INFO.value);
    return 0;
  }
  if (N == 0) {
    return 1;
  }

  // Compute the equilibration matrix R such that
  // inv(R)*A*C has unit 1-norm.

  KD = KU + 1;
  KE = KL + 1;
  if (NOTRANS) {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CMODE == 1) {
        for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
          TMP += (AB[KD + I - J][J] * C[J]).abs();
        }
      } else if (CMODE == 0) {
        for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
          TMP += AB[KD + I - J][J].abs();
        }
      } else {
        for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
          TMP += (AB[KD + I - J][J] / C[J]).abs();
        }
      }
      WORK[2 * N + I] = TMP;
    }
  } else {
    for (I = 1; I <= N; I++) {
      TMP = 0.0;
      if (CMODE == 1) {
        for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
          TMP += (AB[KE - I + J][I] * C[J]).abs();
        }
      } else if (CMODE == 0) {
        for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
          TMP += AB[KE - I + J][I].abs();
        }
      } else {
        for (J = max(I - KL, 1); J <= min(I + KU, N); J++) {
          TMP += (AB[KE - I + J][I] / C[J]).abs();
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

      if (NOTRANS) {
        dgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N),
            N, INFO);
      } else {
        dgbtrs('Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N), N,
            INFO);
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

      if (NOTRANS) {
        dgbtrs('Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N), N,
            INFO);
      } else {
        dgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK.asMatrix(N),
            N, INFO);
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
