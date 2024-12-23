// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlassq.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsyequb(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> S_,
  final Box<double> SCOND,
  final Box<double> AMAX,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final S = S_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  const MAX_ITER = 100;
  int I, J, ITER;
  double AVG = 0,
      STD,
      TOL,
      C0,
      C1,
      C2,
      T,
      U,
      SI,
      D,
      BASE,
      SMIN,
      SMAX,
      SMLNUM,
      BIGNUM;
  bool UP;
  final SCALE = Box(0.0), SUMSQ = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;
  if (!(lsame(UPLO, 'U') || lsame(UPLO, 'L'))) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DSYEQUB', -INFO.value);
    return;
  }

  UP = lsame(UPLO, 'U');
  AMAX.value = ZERO;

  // Quick return if possible.

  if (N == 0) {
    SCOND.value = ONE;
    return;
  }

  for (I = 1; I <= N; I++) {
    S[I] = ZERO;
  }

  AMAX.value = ZERO;
  if (UP) {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J - 1; I++) {
        S[I] = max(S[I], A[I][J].abs());
        S[J] = max(S[J], A[I][J].abs());
        AMAX.value = max(AMAX.value, A[I][J].abs());
      }
      S[J] = max(S[J], A[J][J].abs());
      AMAX.value = max(AMAX.value, A[J][J].abs());
    }
  } else {
    for (J = 1; J <= N; J++) {
      S[J] = max(S[J], A[J][J].abs());
      AMAX.value = max(AMAX.value, A[J][J].abs());
      for (I = J + 1; I <= N; I++) {
        S[I] = max(S[I], A[I][J].abs());
        S[J] = max(S[J], A[I][J].abs());
        AMAX.value = max(AMAX.value, A[I][J].abs());
      }
    }
  }
  for (J = 1; J <= N; J++) {
    S[J] = 1.0 / S[J];
  }

  TOL = ONE / sqrt(2.0 * N);

  for (ITER = 1; ITER <= MAX_ITER; ITER++) {
    SCALE.value = 0.0;
    SUMSQ.value = 0.0;
    // beta = |A|s
    for (I = 1; I <= N; I++) {
      WORK[I] = ZERO;
    }
    if (UP) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J - 1; I++) {
          WORK[I] += A[I][J].abs() * S[J];
          WORK[J] += A[I][J].abs() * S[I];
        }
        WORK[J] += A[J][J].abs() * S[J];
      }
    } else {
      for (J = 1; J <= N; J++) {
        WORK[J] += A[J][J].abs() * S[J];
        for (I = J + 1; I <= N; I++) {
          WORK[I] += A[I][J].abs() * S[J];
          WORK[J] += A[I][J].abs() * S[I];
        }
      }
    }

    // avg = s^T beta / n
    AVG = 0.0;
    for (I = 1; I <= N; I++) {
      AVG += S[I] * WORK[I];
    }
    AVG /= N;

    STD = 0.0;
    for (I = N + 1; I <= 2 * N; I++) {
      WORK[I] = S[I - N] * WORK[I - N] - AVG;
    }
    dlassq(N, WORK(N + 1), 1, SCALE, SUMSQ);
    STD = SCALE.value * sqrt(SUMSQ.value / N);

    if (STD < TOL * AVG) break;

    for (I = 1; I <= N; I++) {
      T = A[I][I].abs();
      SI = S[I];
      C2 = (N - 1) * T;
      C1 = (N - 2) * (WORK[I] - T * SI);
      C0 = -(T * SI) * SI + 2 * WORK[I] * SI - N * AVG;
      D = C1 * C1 - 4 * C0 * C2;

      if (D <= 0) {
        INFO.value = -1;
        return;
      }
      SI = -2 * C0 / (C1 + sqrt(D));

      D = SI - S[I];
      U = ZERO;
      if (UP) {
        for (J = 1; J <= I; J++) {
          T = A[J][I].abs();
          U += S[J] * T;
          WORK[J] += D * T;
        }
        for (J = I + 1; J <= N; J++) {
          T = A[I][J].abs();
          U += S[J] * T;
          WORK[J] += D * T;
        }
      } else {
        for (J = 1; J <= I; J++) {
          T = A[I][J].abs();
          U += S[J] * T;
          WORK[J] += D * T;
        }
        for (J = I + 1; J <= N; J++) {
          T = A[J][I].abs();
          U += S[J] * T;
          WORK[J] += D * T;
        }
      }

      AVG += (U + WORK[I]) * D / N;
      S[I] = SI;
    }
  }

  SMLNUM = dlamch('SAFEMIN');
  BIGNUM = ONE / SMLNUM;
  SMIN = BIGNUM;
  SMAX = ZERO;
  T = ONE / sqrt(AVG);
  BASE = dlamch('B');
  U = ONE / log(BASE);
  for (I = 1; I <= N; I++) {
    S[I] = pow(BASE, (U * log(S[I] * T)).toInt()).toDouble();
    SMIN = min(SMIN, S[I]);
    SMAX = max(SMAX, S[I]);
  }
  SCOND.value = max(SMIN, SMLNUM) / min(SMAX, BIGNUM);
}
