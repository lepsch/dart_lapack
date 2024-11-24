// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);

  const ONE = 1.0;
  bool UPPER;
  int J, K, KP;
  double AK, AKM1, AKM1K, BK, BKM1, DENOM;

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
    xerbla('DSYTRS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Solve A*X = B, where A = U*D*U**T.

    // First solve U*D*X = B, overwriting B with X.

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = N;

    // If K < 1, exit from loop.

    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(U(K)), where U(K) is the transformation
        // stored in column K of A.

        dger(K - 1, NRHS, -ONE, A(1, K).asArray(), 1, B(K, 1).asArray(), LDB,
            B(1, 1), LDB);

        // Multiply by the inverse of the diagonal block.

        dscal(NRHS, ONE / A[K][K], B(K, 1).asArray(), LDB);
        K--;
      } else {
        // 2 x 2 diagonal block

        // Interchange rows K-1 and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K - 1) {
          dswap(NRHS, B(K - 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(U(K)), where U(K) is the transformation
        // stored in columns K-1 and K of A.

        dger(K - 2, NRHS, -ONE, A(1, K).asArray(), 1, B(K, 1).asArray(), LDB,
            B(1, 1), LDB);
        dger(K - 2, NRHS, -ONE, A(1, K - 1).asArray(), 1, B(K - 1, 1).asArray(),
            LDB, B(1, 1), LDB);

        // Multiply by the inverse of the diagonal block.

        AKM1K = A[K - 1][K];
        AKM1 = A[K - 1][K - 1] / AKM1K;
        AK = A[K][K] / AKM1K;
        DENOM = AKM1 * AK - ONE;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[K - 1][J] / AKM1K;
          BK = B[K][J] / AKM1K;
          B[K - 1][J] = (AK * BKM1 - BK) / DENOM;
          B[K][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        K -= 2;
      }
    }

    // Next solve U**T *X = B, overwriting B with X.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;

    // If K > N, exit from loop.

    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Multiply by inv(U**T(K)), where U(K) is the transformation
        // stored in column K of A.

        dgemv('Transpose', K - 1, NRHS, -ONE, B, LDB, A(1, K).asArray(), 1, ONE,
            B(K, 1).asArray(), LDB);

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K++;
      } else {
        // 2 x 2 diagonal block

        // Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
        // stored in columns K and K+1 of A.

        dgemv('Transpose', K - 1, NRHS, -ONE, B, LDB, A(1, K).asArray(), 1, ONE,
            B(K, 1).asArray(), LDB);
        dgemv('Transpose', K - 1, NRHS, -ONE, B, LDB, A(1, K + 1).asArray(), 1,
            ONE, B(K + 1, 1).asArray(), LDB);

        // Interchange rows K and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K += 2;
      }
    }
  } else {
    // Solve A*X = B, where A = L*D*L**T.

    // First solve L*D*X = B, overwriting B with X.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;

    // If K > N, exit from loop.

    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(L(K)), where L(K) is the transformation
        // stored in column K of A.

        if (K < N) {
          dger(N - K, NRHS, -ONE, A(K + 1, K).asArray(), 1, B(K, 1).asArray(),
              LDB, B(K + 1, 1), LDB);
        }

        // Multiply by the inverse of the diagonal block.

        dscal(NRHS, ONE / A[K][K], B(K, 1).asArray(), LDB);
        K++;
      } else {
        // 2 x 2 diagonal block

        // Interchange rows K+1 and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K + 1) {
          dswap(NRHS, B(K + 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(L(K)), where L(K) is the transformation
        // stored in columns K and K+1 of A.

        if (K < N - 1) {
          dger(N - K - 1, NRHS, -ONE, A(K + 2, K).asArray(), 1,
              B(K, 1).asArray(), LDB, B(K + 2, 1), LDB);
          dger(N - K - 1, NRHS, -ONE, A(K + 2, K + 1).asArray(), 1,
              B(K + 1, 1).asArray(), LDB, B(K + 2, 1), LDB);
        }

        // Multiply by the inverse of the diagonal block.

        AKM1K = A[K + 1][K];
        AKM1 = A[K][K] / AKM1K;
        AK = A[K + 1][K + 1] / AKM1K;
        DENOM = AKM1 * AK - ONE;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[K][J] / AKM1K;
          BK = B[K + 1][J] / AKM1K;
          B[K][J] = (AK * BKM1 - BK) / DENOM;
          B[K + 1][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        K += 2;
      }
    }

    // Next solve L**T *X = B, overwriting B with X.

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = N;

    // If K < 1, exit from loop.

    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Multiply by inv(L**T(K)), where L(K) is the transformation
        // stored in column K of A.

        if (K < N) {
          dgemv('Transpose', N - K, NRHS, -ONE, B(K + 1, 1), LDB,
              A(K + 1, K).asArray(), 1, ONE, B(K, 1).asArray(), LDB);
        }

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K--;
      } else {
        // 2 x 2 diagonal block

        // Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
        // stored in columns K-1 and K of A.

        if (K < N) {
          dgemv('Transpose', N - K, NRHS, -ONE, B(K + 1, 1), LDB,
              A(K + 1, K).asArray(), 1, ONE, B(K, 1).asArray(), LDB);
          dgemv('Transpose', N - K, NRHS, -ONE, B(K + 1, 1), LDB,
              A(K + 1, K - 1).asArray(), 1, ONE, B(K - 1, 1).asArray(), LDB);
        }

        // Interchange rows K and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K) {
          dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K -= 2;
      }
    }
  }
}
