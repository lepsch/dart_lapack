import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlassq.dart';

void zsyequb(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> S_,
  final Box<double> SCOND,
  final Box<double> AMAX,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
    xerbla('ZSYEQUB', -INFO.value);
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
        S[I] = max(S[I], A[I][J].cabs1());
        S[J] = max(S[J], A[I][J].cabs1());
        AMAX.value = max(AMAX.value, A[I][J].cabs1());
      }
      S[J] = max(S[J], A[J][J].cabs1());
      AMAX.value = max(AMAX.value, A[J][J].cabs1());
    }
  } else {
    for (J = 1; J <= N; J++) {
      S[J] = max(S[J], A[J][J].cabs1());
      AMAX.value = max(AMAX.value, A[J][J].cabs1());
      for (I = J + 1; I <= N; I++) {
        S[I] = max(S[I], A[I][J].cabs1());
        S[J] = max(S[J], A[I][J].cabs1());
        AMAX.value = max(AMAX.value, A[I][J].cabs1());
      }
    }
  }
  for (J = 1; J <= N; J++) {
    S[J] = 1.0 / S[J];
  }

  TOL = ONE / sqrt(2.0 * N);
  mainLoop:
  for (ITER = 1; ITER <= MAX_ITER; ITER++) {
    SCALE.value = 0.0;
    SUMSQ.value = 0.0;
    // beta = |A|s
    for (I = 1; I <= N; I++) {
      WORK[I] = Complex.zero;
    }
    if (UP) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J - 1; I++) {
          WORK[I] += (A[I][J].cabs1() * S[J]).toComplex();
          WORK[J] += (A[I][J].cabs1() * S[I]).toComplex();
        }
        WORK[J] += (A[J][J].cabs1() * S[J]).toComplex();
      }
    } else {
      for (J = 1; J <= N; J++) {
        WORK[J] += (A[J][J].cabs1() * S[J]).toComplex();
        for (I = J + 1; I <= N; I++) {
          WORK[I] += (A[I][J].cabs1() * S[J]).toComplex();
          WORK[J] += (A[I][J].cabs1() * S[I]).toComplex();
        }
      }
    }

    // avg = s^T beta / n
    AVG = 0.0;
    for (I = 1; I <= N; I++) {
      AVG += S[I] * WORK[I].real;
    }
    AVG /= N;

    STD = 0.0;
    for (I = N + 1; I <= 2 * N; I++) {
      WORK[I] = S[I - N].toComplex() * WORK[I - N] - AVG.toComplex();
    }
    zlassq(N, WORK(N + 1), 1, SCALE, SUMSQ);
    STD = SCALE.value * sqrt(SUMSQ.value / N);

    if (STD < TOL * AVG) break mainLoop;

    for (I = 1; I <= N; I++) {
      T = A[I][I].cabs1();
      SI = S[I];
      C2 = (N - 1) * T;
      C1 = (N - 2) * (WORK[I].real - T * SI);
      C0 = -(T * SI) * SI + 2 * WORK[I].real * SI - N * AVG;
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
          T = A[J][I].cabs1();
          U += S[J] * T;
          WORK[J] += (D * T).toComplex();
        }
        for (J = I + 1; J <= N; J++) {
          T = A[I][J].cabs1();
          U += S[J] * T;
          WORK[J] += (D * T).toComplex();
        }
      } else {
        for (J = 1; J <= I; J++) {
          T = A[I][J].cabs1();
          U += S[J] * T;
          WORK[J] += (D * T).toComplex();
        }
        for (J = I + 1; J <= N; J++) {
          T = (A[J][I].cabs1());
          U += S[J] * T;
          WORK[J] += (D * T).toComplex();
        }
      }

      AVG += (U + WORK[I].real) * D / N;
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
