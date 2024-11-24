// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zhetrs_3(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> E_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final E = E_.having();
  const ONE = 1.0;
  bool UPPER;
  int I, J, K, KP;
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
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('ZHETRS_3', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  if (UPPER) {
    // Begin Upper

    // Solve A*X = B, where A = U*D*U**H.

    // P**T * B

    // Interchange rows K and IPIV(K) of matrix B in the same order
    // that the formation order of IPIV(I) vector for Upper case.

    // (We can do the simple loop over IPIV with decrement -1,
    // since the ABS value of IPIV(I) represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = N; K >= 1; K--) {
      KP = IPIV[K].abs();
      if (KP != K) {
        zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }

    // Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

    ztrsm('L', 'U', 'N', 'U', N, NRHS, Complex.one, A, LDA, B, LDB);

    // Compute D \ B -> B   [ D \ (U \P**T * B) ]

    I = N;
    while (I >= 1) {
      if (IPIV[I] > 0) {
        S = ONE / A[I][I].real;
        zdscal(NRHS, S, B(I, 1).asArray(), LDB);
      } else if (I > 1) {
        AKM1K = E[I];
        AKM1 = A[I - 1][I - 1] / AKM1K;
        AK = A[I][I] / AKM1K.conjugate();
        DENOM = AKM1 * AK - Complex.one;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[I - 1][J] / AKM1K;
          BK = B[I][J] / AKM1K.conjugate();
          B[I - 1][J] = (AK * BKM1 - BK) / DENOM;
          B[I][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        I--;
      }
      I--;
    }

    // Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ]

    ztrsm('L', 'U', 'C', 'U', N, NRHS, Complex.one, A, LDA, B, LDB);

    // P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ]

    // Interchange rows K and IPIV(K) of matrix B in reverse order
    // from the formation order of IPIV(I) vector for Upper case.

    // (We can do the simple loop over IPIV with increment 1,
    // since the ABS value of IPIV(I) represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = 1; K <= N; K += 1) {
      KP = IPIV[K].abs();
      if (KP != K) {
        zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }
  } else {
    // Begin Lower

    // Solve A*X = B, where A = L*D*L**H.

    // P**T * B
    // Interchange rows K and IPIV(K) of matrix B in the same order
    // that the formation order of IPIV(I) vector for Lower case.

    // (We can do the simple loop over IPIV with increment 1,
    // since the ABS value of IPIV(I) represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = 1; K <= N; K += 1) {
      KP = IPIV[K].abs();
      if (KP != K) {
        zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }

    // Compute (L \P**T * B) -> B    [ (L \P**T * B) ]

    ztrsm('L', 'L', 'N', 'U', N, NRHS, Complex.one, A, LDA, B, LDB);

    // Compute D \ B -> B   [ D \ (L \P**T * B) ]

    I = 1;
    while (I <= N) {
      if (IPIV[I] > 0) {
        S = ONE / A[I][I].real;
        zdscal(NRHS, S, B(I, 1).asArray(), LDB);
      } else if (I < N) {
        AKM1K = E[I];
        AKM1 = A[I][I] / AKM1K.conjugate();
        AK = A[I + 1][I + 1] / AKM1K;
        DENOM = AKM1 * AK - Complex.one;
        for (J = 1; J <= NRHS; J++) {
          BKM1 = B[I][J] / AKM1K.conjugate();
          BK = B[I + 1][J] / AKM1K;
          B[I][J] = (AK * BKM1 - BK) / DENOM;
          B[I + 1][J] = (AKM1 * BK - BKM1) / DENOM;
        }
        I++;
      }
      I++;
    }

    // Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ]

    ztrsm('L', 'L', 'C', 'U', N, NRHS, Complex.one, A, LDA, B, LDB);

    // P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ]

    // Interchange rows K and IPIV(K) of matrix B in reverse order
    // from the formation order of IPIV(I) vector for Lower case.

    // (We can do the simple loop over IPIV with decrement -1,
    // since the ABS value of IPIV(I) represents the row index
    // of the interchange with row i in both 1x1 and 2x2 pivot cases)

    for (K = N; K >= 1; K--) {
      KP = IPIV[K].abs();
      if (KP != K) {
        zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
      }
    }
  }
}
