// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/blas/dtrsm.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsytrs_3(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> E_,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  const ONE = 1.0;
  bool UPPER;
  int I, J, K, KP;
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
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DSYTRS_3', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Begin Upper

    // Solve A*X = B, where A = U*D*U**T.

    // P**T * B

    // Interchange rows K and IPIV(K) of matrix B in the same order
    // that the formation order of IPIV(I) vector for Upper case.

    // (We can do the simple loop over IPIV with decrement -1,
    // since the ABS value of IPIV[I] represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = N; K >= 1; K--) {
      KP = IPIV[K].abs();
      if (KP != K) {
        dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }

    // Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

    dtrsm('L', 'U', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB);

    // Compute D \ B -> B   [ D \ (U \P**T * B) ]

    I = N;
    while (I >= 1) {
      if (IPIV[I] > 0) {
        dscal(NRHS, ONE / A[I][I], B(I, 1).asArray(), LDB);
      } else if (I > 1) {
        AKM1K = E[I];
        AKM1 = A[I - 1][I - 1] / AKM1K;
        AK = A[I][I] / AKM1K;
        DENOM = AKM1 * AK - ONE;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[I - 1][J] / AKM1K;
          BK = B[I][J] / AKM1K;
          B[I - 1][J] = (AK * BKM1 - BK) / DENOM;
          B[I][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        I--;
      }
      I--;
    }

    // Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]

    dtrsm('L', 'U', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB);

    // P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]

    // Interchange rows K and IPIV(K) of matrix B in reverse order
    // from the formation order of IPIV(I) vector for Upper case.

    // (We can do the simple loop over IPIV with increment 1,
    // since the ABS value of IPIV(I) represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = 1; K <= N; K++) {
      KP = IPIV[K].abs();
      if (KP != K) {
        dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }
  } else {
    // Begin Lower

    // Solve A*X = B, where A = L*D*L**T.

    // P**T * B
    // Interchange rows K and IPIV(K) of matrix B in the same order
    // that the formation order of IPIV(I) vector for Lower case.

    // (We can do the simple loop over IPIV with increment 1,
    // since the ABS value of IPIV(I) represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = 1; K <= N; K++) {
      KP = IPIV[K].abs();
      if (KP != K) {
        dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }

    // Compute (L \P**T * B) -> B    [ (L \P**T * B) ]

    dtrsm('L', 'L', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB);

    // Compute D \ B -> B   [ D \ (L \P**T * B) ]

    I = 1;
    while (I <= N) {
      if (IPIV[I] > 0) {
        dscal(NRHS, ONE / A[I][I], B(I, 1).asArray(), LDB);
      } else if (I < N) {
        AKM1K = E[I];
        AKM1 = A[I][I] / AKM1K;
        AK = A[I + 1][I + 1] / AKM1K;
        DENOM = AKM1 * AK - ONE;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[I][J] / AKM1K;
          BK = B[I + 1][J] / AKM1K;
          B[I][J] = (AK * BKM1 - BK) / DENOM;
          B[I + 1][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        I++;
      }
      I++;
    }

    // Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]

    dtrsm('L', 'L', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB);

    // P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]

    // Interchange rows K and IPIV(K) of matrix B in reverse order
    // from the formation order of IPIV(I) vector for Lower case.

    // (We can do the simple loop over IPIV with decrement -1,
    // since the ABS value of IPIV(I) represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = N; K >= 1; K--) {
      KP = IPIV[K].abs();
      if (KP != K) {
        dswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }

    // END Lower
  }
}
