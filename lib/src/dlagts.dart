import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlagts(
  final int JOB,
  final int N,
  final Array<double> A_,
  final Array<double> B_,
  final Array<double> C_,
  final Array<double> D_,
  final Array<int> IN_,
  final Array<double> Y_,
  final Box<double> TOL,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final B = B_.having();
  final C = C_.having();
  final D = D_.having();
  final IN = IN_.having();
  final Y = Y_.having();
  const ONE = 1.0, ZERO = 0.0;
  int K = 0;
  double ABSAK = 0, AK = 0, BIGNUM = 0, EPS = 0, PERT = 0, SFMIN = 0, TEMP = 0;

  INFO.value = 0;
  if (((JOB).abs() > 2) || (JOB == 0)) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('DLAGTS', -INFO.value);
    return;
  }

  if (N == 0) return;

  EPS = dlamch('Epsilon');
  SFMIN = dlamch('Safe minimum');
  BIGNUM = ONE / SFMIN;

  if (JOB < 0) {
    if (TOL.value <= ZERO) {
      TOL.value = (A[1]).abs();
      if (N > 1) TOL.value = max(TOL.value, max((A[2]).abs(), (B[1]).abs()));
      for (K = 3; K <= N; K++) {
        TOL.value = max(max(TOL.value, (A[K]).abs()),
            max((B[K - 1]).abs(), (D[K - 2]).abs()));
      }
      TOL.value *= EPS;
      if (TOL.value == ZERO) TOL.value = EPS;
    }
  }

  if ((JOB).abs() == 1) {
    for (K = 2; K <= N; K++) {
      if (IN[K - 1] == 0) {
        Y[K] -= C[K - 1] * Y[K - 1];
      } else {
        TEMP = Y[K - 1];
        Y[K - 1] = Y[K];
        Y[K] = TEMP - C[K - 1] * Y[K];
      }
    }
    if (JOB == 1) {
      for (K = N; K >= 1; K--) {
        if (K <= N - 2) {
          TEMP = Y[K] - B[K] * Y[K + 1] - D[K] * Y[K + 2];
        } else if (K == N - 1) {
          TEMP = Y[K] - B[K] * Y[K + 1];
        } else {
          TEMP = Y[K];
        }
        AK = A[K];
        ABSAK = (AK).abs();
        if (ABSAK < ONE) {
          if (ABSAK < SFMIN) {
            if (ABSAK == ZERO || (TEMP).abs() * SFMIN > ABSAK) {
              INFO.value = K;
              return;
            } else {
              TEMP *= BIGNUM;
              AK *= BIGNUM;
            }
          } else if ((TEMP).abs() > ABSAK * BIGNUM) {
            INFO.value = K;
            return;
          }
        }
        Y[K] = TEMP / AK;
      }
    } else {
      for (K = N; K >= 1; K--) {
        if (K <= N - 2) {
          TEMP = Y[K] - B[K] * Y[K + 1] - D[K] * Y[K + 2];
        } else if (K == N - 1) {
          TEMP = Y[K] - B[K] * Y[K + 1];
        } else {
          TEMP = Y[K];
        }
        AK = A[K];
        PERT = sign(TOL.value, AK).toDouble();
        //  }
        while (true) {
          ABSAK = (AK).abs();
          if (ABSAK < ONE) {
            if (ABSAK < SFMIN) {
              if (ABSAK == ZERO || (TEMP).abs() * SFMIN > ABSAK) {
                AK += PERT;
                PERT = 2 * PERT;
                continue;
              } else {
                TEMP *= BIGNUM;
                AK *= BIGNUM;
              }
            } else if ((TEMP).abs() > ABSAK * BIGNUM) {
              AK += PERT;
              PERT = 2 * PERT;
              continue;
            }
          }
          break;
        }
        Y[K] = TEMP / AK;
      }
    }
  } else {
    // Come to here if  JOB = 2 or -2

    if (JOB == 2) {
      for (K = 1; K <= N; K++) {
        if (K >= 3) {
          TEMP = Y[K] - B[K - 1] * Y[K - 1] - D[K - 2] * Y[K - 2];
        } else if (K == 2) {
          TEMP = Y[K] - B[K - 1] * Y[K - 1];
        } else {
          TEMP = Y[K];
        }
        AK = A[K];
        ABSAK = (AK).abs();
        if (ABSAK < ONE) {
          if (ABSAK < SFMIN) {
            if (ABSAK == ZERO || (TEMP).abs() * SFMIN > ABSAK) {
              INFO.value = K;
              return;
            } else {
              TEMP *= BIGNUM;
              AK *= BIGNUM;
            }
          } else if ((TEMP).abs() > ABSAK * BIGNUM) {
            INFO.value = K;
            return;
          }
        }
        Y[K] = TEMP / AK;
      }
    } else {
      for (K = 1; K <= N; K++) {
        if (K >= 3) {
          TEMP = Y[K] - B[K - 1] * Y[K - 1] - D[K - 2] * Y[K - 2];
        } else if (K == 2) {
          TEMP = Y[K] - B[K - 1] * Y[K - 1];
        } else {
          TEMP = Y[K];
        }
        AK = A[K];
        PERT = sign(TOL.value, AK).toDouble();
        while (true) {
          ABSAK = (AK).abs();
          if (ABSAK < ONE) {
            if (ABSAK < SFMIN) {
              if (ABSAK == ZERO || (TEMP).abs() * SFMIN > ABSAK) {
                AK += PERT;
                PERT = 2 * PERT;
                continue;
              } else {
                TEMP *= BIGNUM;
                AK *= BIGNUM;
              }
            } else if ((TEMP).abs() > ABSAK * BIGNUM) {
              AK += PERT;
              PERT = 2 * PERT;
              continue;
            }
          }
          break;
        }
        Y[K] = TEMP / AK;
      }
    }

    for (K = N; K >= 2; K--) {
      if (IN[K - 1] == 0) {
        Y[K - 1] -= C[K - 1] * Y[K];
      } else {
        TEMP = Y[K - 1];
        Y[K - 1] = Y[K];
        Y[K] = TEMP - C[K - 1] * Y[K];
      }
    }
  }
}
