import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgebal(
  final String JOB,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> SCALE_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final SCALE = SCALE_.having();
  const ZERO = 0.0, ONE = 1.0;
  const SCLFAC = 2.0;
  const FACTOR = 0.95;
  bool NOCONV, CANSWAP;
  int I, ICA, IRA, J, K, L;
  double C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, SFMIN2;

  // Test the input parameters

  INFO.value = 0;
  if (!lsame(JOB, 'N') &&
      !lsame(JOB, 'P') &&
      !lsame(JOB, 'S') &&
      !lsame(JOB, 'B')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZGEBAL', -INFO.value);
    return;
  }

  // Quick returns.

  if (N == 0) {
    ILO.value = 1;
    IHI.value = 0;
    return;
  }

  if (lsame(JOB, 'N')) {
    for (I = 1; I <= N; I++) {
      SCALE[I] = ONE;
    }
    ILO.value = 1;
    IHI.value = N;
    return;
  }

  // Permutation to isolate eigenvalues if possible.

  K = 1;
  L = N;

  if (!lsame(JOB, 'S')) {
    // Row and column exchange.

    NOCONV = true;
    while (NOCONV) {
      // Search for rows isolating an eigenvalue and push them down.

      NOCONV = false;
      for (I = L; I >= 1; I--) {
        CANSWAP = true;
        for (J = 1; J <= L; J++) {
          if (I != J && (A[I][J].real != ZERO || A[I][J].imaginary != ZERO)) {
            CANSWAP = false;
            break;
          }
        }

        if (CANSWAP) {
          SCALE[L] = I.toDouble();
          if (I != L) {
            zswap(L, A(1, I).asArray(), 1, A(1, L).asArray(), 1);
            zswap(N - K + 1, A(I, K).asArray(), LDA, A(L, K).asArray(), LDA);
          }
          NOCONV = true;

          if (L == 1) {
            ILO.value = 1;
            IHI.value = 1;
            return;
          }

          L--;
        }
      }
    }

    NOCONV = true;
    while (NOCONV) {
      // Search for columns isolating an eigenvalue and push them left.

      NOCONV = false;
      for (J = K; J <= L; J++) {
        CANSWAP = true;
        for (I = K; I <= L; I++) {
          if (I != J && (A[I][J].real != ZERO || A[I][J].imaginary != ZERO)) {
            CANSWAP = false;
            break;
          }
        }

        if (CANSWAP) {
          SCALE[K] = J.toDouble();
          if (J != K) {
            zswap(L, A(1, J).asArray(), 1, A(1, K).asArray(), 1);
            zswap(N - K + 1, A(J, K).asArray(), LDA, A(K, K).asArray(), LDA);
          }
          NOCONV = true;

          K++;
        }
      }
    }
  }

  // Initialize SCALE for non-permuted submatrix.

  for (I = K; I <= L; I++) {
    SCALE[I] = ONE;
  }

  // If we only had to permute, we are done.

  if (lsame(JOB, 'P')) {
    ILO.value = K;
    IHI.value = L;
    return;
  }

  // Balance the submatrix in rows K to L.

  // Iterative loop for norm reduction.

  SFMIN1 = dlamch('S') / dlamch('P');
  SFMAX1 = ONE / SFMIN1;
  SFMIN2 = SFMIN1 * SCLFAC;
  SFMAX2 = ONE / SFMIN2;

  NOCONV = true;
  while (NOCONV) {
    NOCONV = false;

    for (I = K; I <= L; I++) {
      C = dznrm2(L - K + 1, A(K, I).asArray(), 1);
      R = dznrm2(L - K + 1, A(I, K).asArray(), LDA);
      ICA = izamax(L, A(1, I).asArray(), 1);
      CA = A[ICA][I].abs();
      IRA = izamax(N - K + 1, A(I, K).asArray(), LDA);
      RA = A[I][IRA + K - 1].abs();

      // Guard against zero C or R due to underflow.

      if (C == ZERO || R == ZERO) continue;

      // Exit if NaN to avoid infinite loop

      if (disnan(C + CA + R + RA)) {
        INFO.value = -3;
        xerbla('ZGEBAL', -INFO.value);
        return;
      }

      G = R / SCLFAC;
      F = ONE;
      S = C + R;

      while (
          C < G && max(F, max(C, CA)) < SFMAX2 && min(R, min(G, RA)) > SFMIN2) {
        F *= SCLFAC;
        C *= SCLFAC;
        CA *= SCLFAC;
        R /= SCLFAC;
        G /= SCLFAC;
        RA /= SCLFAC;
      }

      G = C / SCLFAC;

      while (G >= R &&
          max(R, RA) < SFMAX2 &&
          min(min(F, C), min(G, CA)) > SFMIN2) {
        F /= SCLFAC;
        C /= SCLFAC;
        G /= SCLFAC;
        CA /= SCLFAC;
        R *= SCLFAC;
        RA *= SCLFAC;
      }

      // Now balance.

      if ((C + R) >= FACTOR * S) continue;
      if (F < ONE && SCALE[I] < ONE) {
        if (F * SCALE[I] <= SFMIN1) continue;
      }
      if (F > ONE && SCALE[I] > ONE) {
        if (SCALE[I] >= SFMAX1 / F) continue;
      }
      G = ONE / F;
      SCALE[I] *= F;
      NOCONV = true;

      zdscal(N - K + 1, G, A(I, K).asArray(), LDA);
      zdscal(L, F, A(1, I).asArray(), 1);
    }
  }

  ILO.value = K;
  IHI.value = L;
}
