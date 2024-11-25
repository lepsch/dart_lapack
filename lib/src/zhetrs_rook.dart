// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zgeru.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';

void zhetrs_rook(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  const ONE = 1.0;
  bool UPPER;
  int J, K, KP;
  double S;
  Complex AK, AKM1, AKM1K, BK, BKM1, DENOM;

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZHETRS_ROOK', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B, where A = U*D*U**H.

    // First solve U*D*X = B, overwriting B with X.

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = N;
    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(U(K)), where U(K) is the transformation
        // stored in column K of A.

        zgeru(K - 1, NRHS, -Complex.one, A(1, K).asArray(), 1,
            B(K, 1).asArray(), LDB, B(1, 1), LDB);

        // Multiply by the inverse of the diagonal block.

        S = ONE / A[K][K].real;
        zdscal(NRHS, S, B(K, 1).asArray(), LDB);
        K--;
      } else {
        // 2 x 2 diagonal block

        // Interchange rows K and -IPIV(K), then K-1 and -IPIV(K-1)

        KP = -IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        KP = -IPIV[K - 1];
        if (KP != K - 1) {
          zswap(NRHS, B(K - 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(U(K)), where U(K) is the transformation
        // stored in columns K-1 and K of A.

        zgeru(K - 2, NRHS, -Complex.one, A(1, K).asArray(), 1,
            B(K, 1).asArray(), LDB, B(1, 1), LDB);
        zgeru(K - 2, NRHS, -Complex.one, A(1, K - 1).asArray(), 1,
            B(K - 1, 1).asArray(), LDB, B(1, 1), LDB);

        // Multiply by the inverse of the diagonal block.

        AKM1K = A[K - 1][K];
        AKM1 = A[K - 1][K - 1] / AKM1K;
        AK = A[K][K] / AKM1K.conjugate();
        DENOM = AKM1 * AK - Complex.one;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[K - 1][J] / AKM1K;
          BK = B[K][J] / AKM1K.conjugate();
          B[K - 1][J] = (AK * BKM1 - BK) / DENOM;
          B[K][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        K -= 2;
      }
    }

    // Next solve U**H *X = B, overwriting B with X.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Multiply by inv(U**H(K)), where U(K) is the transformation
        // stored in column K of A.

        if (K > 1) {
          zlacgv(NRHS, B(K, 1).asArray(), LDB);
          zgemv('Conjugate transpose', K - 1, NRHS, -Complex.one, B, LDB,
              A(1, K).asArray(), 1, Complex.one, B(K, 1).asArray(), LDB);
          zlacgv(NRHS, B(K, 1).asArray(), LDB);
        }

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K++;
      } else {
        // 2 x 2 diagonal block

        // Multiply by inv(U**H(K+1)), where U(K+1) is the transformation
        // stored in columns K and K+1 of A.

        if (K > 1) {
          zlacgv(NRHS, B(K, 1).asArray(), LDB);
          zgemv('Conjugate transpose', K - 1, NRHS, -Complex.one, B, LDB,
              A(1, K).asArray(), 1, Complex.one, B(K, 1).asArray(), LDB);
          zlacgv(NRHS, B(K, 1).asArray(), LDB);

          zlacgv(NRHS, B(K + 1, 1).asArray(), LDB);
          zgemv(
              'Conjugate transpose',
              K - 1,
              NRHS,
              -Complex.one,
              B,
              LDB,
              A(1, K + 1).asArray(),
              1,
              Complex.one,
              B(K + 1, 1).asArray(),
              LDB);
          zlacgv(NRHS, B(K + 1, 1).asArray(), LDB);
        }

        // Interchange rows K and -IPIV(K), then K+1 and -IPIV(K+1)

        KP = -IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        KP = -IPIV[K + 1];
        if (KP != K + 1) {
          zswap(NRHS, B(K + 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        K += 2;
      }
    }
  } else {
    // Solve A*X = B, where A = L*D*L**H.

    // First solve L*D*X = B, overwriting B with X.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(L(K)), where L(K) is the transformation
        // stored in column K of A.

        if (K < N) {
          zgeru(N - K, NRHS, -Complex.one, A(K + 1, K).asArray(), 1,
              B(K, 1).asArray(), LDB, B(K + 1, 1), LDB);
        }

        // Multiply by the inverse of the diagonal block.

        S = ONE / A[K][K].real;
        zdscal(NRHS, S, B(K, 1).asArray(), LDB);
        K++;
      } else {
        // 2 x 2 diagonal block

        // Interchange rows K and -IPIV(K), then K+1 and -IPIV(K+1)

        KP = -IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        KP = -IPIV[K + 1];
        if (KP != K + 1) {
          zswap(NRHS, B(K + 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(L(K)), where L(K) is the transformation
        // stored in columns K and K+1 of A.

        if (K < N - 1) {
          zgeru(N - K - 1, NRHS, -Complex.one, A(K + 2, K).asArray(), 1,
              B(K, 1).asArray(), LDB, B(K + 2, 1), LDB);
          zgeru(N - K - 1, NRHS, -Complex.one, A(K + 2, K + 1).asArray(), 1,
              B(K + 1, 1).asArray(), LDB, B(K + 2, 1), LDB);
        }

        // Multiply by the inverse of the diagonal block.

        AKM1K = A[K + 1][K];
        AKM1 = A[K][K] / AKM1K.conjugate();
        AK = A[K + 1][K + 1] / AKM1K;
        DENOM = AKM1 * AK - Complex.one;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[K][J] / AKM1K.conjugate();
          BK = B[K + 1][J] / AKM1K;
          B[K][J] = (AK * BKM1 - BK) / DENOM;
          B[K + 1][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        K += 2;
      }
    }

    // Next solve L**H *X = B, overwriting B with X.

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = N;
    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Multiply by inv(L**H(K)), where L(K) is the transformation
        // stored in column K of A.

        if (K < N) {
          zlacgv(NRHS, B(K, 1).asArray(), LDB);
          zgemv(
              'Conjugate transpose',
              N - K,
              NRHS,
              -Complex.one,
              B(K + 1, 1),
              LDB,
              A(K + 1, K).asArray(),
              1,
              Complex.one,
              B(K, 1).asArray(),
              LDB);
          zlacgv(NRHS, B(K, 1).asArray(), LDB);
        }

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K--;
      } else {
        // 2 x 2 diagonal block

        // Multiply by inv(L**H(K-1)), where L(K-1) is the transformation
        // stored in columns K-1 and K of A.

        if (K < N) {
          zlacgv(NRHS, B(K, 1).asArray(), LDB);
          zgemv(
              'Conjugate transpose',
              N - K,
              NRHS,
              -Complex.one,
              B(K + 1, 1),
              LDB,
              A(K + 1, K).asArray(),
              1,
              Complex.one,
              B(K, 1).asArray(),
              LDB);
          zlacgv(NRHS, B(K, 1).asArray(), LDB);

          zlacgv(NRHS, B(K - 1, 1).asArray(), LDB);
          zgemv(
              'Conjugate transpose',
              N - K,
              NRHS,
              -Complex.one,
              B(K + 1, 1),
              LDB,
              A(K + 1, K - 1).asArray(),
              1,
              Complex.one,
              B(K - 1, 1).asArray(),
              LDB);
          zlacgv(NRHS, B(K - 1, 1).asArray(), LDB);
        }

        // Interchange rows K and -IPIV(K), then K-1 and -IPIV(K-1)

        KP = -IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        KP = -IPIV[K - 1];
        if (KP != K - 1) {
          zswap(NRHS, B(K - 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        K -= 2;
      }
    }
  }
}
