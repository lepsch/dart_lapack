// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/blas/ztrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zgtsv.dart';
import 'package:dart_lapack/src/zlacgv.dart';
import 'package:dart_lapack/src/zlacpy.dart';

void zhetrs_aa(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  bool LQUERY, UPPER;
  int K, KP, LWKMIN;

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);
  if (min(N, NRHS) == 0) {
    LWKMIN = 1;
  } else {
    LWKMIN = 3 * N - 2;
  }

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
  } else if (LWORK < LWKMIN && !LQUERY) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZHETRS_AA', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWKMIN.toComplex();
    return;
  }

  // Quick return if possible

  if (min(N, NRHS) == 0) return;

  if (UPPER) {
    // Solve A*X = B, where A = U**H*T*U.

    // 1) Forward substitution with U**H

    if (N > 1) {
      // Pivot, P**T * B -> B

      for (K = 1; K <= N; K++) {
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }

      // Compute U**H \ B -> B    [ (U**H \P**T * B) ]

      ztrsm('L', 'U', 'C', 'U', N - 1, NRHS, Complex.one, A(1, 2), LDA, B(2, 1),
          LDB);
    }

    // 2) Solve with triangular matrix T

    // Compute T \ B -> B   [ T \ (U**H \P**T * B) ]

    zlacpy('F', 1, N, A(1, 1), LDA + 1, WORK(N).asMatrix(), 1);
    if (N > 1) {
      zlacpy('F', 1, N - 1, A(1, 2), LDA + 1, WORK(2 * N).asMatrix(), 1);
      zlacpy('F', 1, N - 1, A(1, 2), LDA + 1, WORK(1).asMatrix(), 1);
      zlacgv(N - 1, WORK(1), 1);
    }
    zgtsv(N, NRHS, WORK(1), WORK(N), WORK(2 * N), B, LDB, INFO);

    // 3) Backward substitution with U

    if (N > 1) {
      // Compute U \ B -> B   [ U \ (T \ (U**H \P**T * B) ) ]

      ztrsm('L', 'U', 'N', 'U', N - 1, NRHS, Complex.one, A(1, 2), LDA, B(2, 1),
          LDB);

      // Pivot, P * B  [ P * (U**H \ (T \ (U \P**T * B) )) ]

      for (K = N; K >= 1; K--) {
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }
    }
  } else {
    // Solve A*X = B, where A = L*T*L**H.

    // 1) Forward substitution with L

    if (N > 1) {
      // Pivot, P**T * B -> B

      for (K = 1; K <= N; K++) {
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }

      // Compute L \ B -> B    [ (L \P**T * B) ]

      ztrsm('L', 'L', 'N', 'U', N - 1, NRHS, Complex.one, A(2, 1), LDA, B(2, 1),
          LDB);
    }

    // 2) Solve with triangular matrix T

    // Compute T \ B -> B   [ T \ (L \P**T * B) ]

    zlacpy('F', 1, N, A(1, 1), LDA + 1, WORK(N).asMatrix(), 1);
    if (N > 1) {
      zlacpy('F', 1, N - 1, A(2, 1), LDA + 1, WORK(1).asMatrix(), 1);
      zlacpy('F', 1, N - 1, A(2, 1), LDA + 1, WORK(2 * N).asMatrix(), 1);
      zlacgv(N - 1, WORK(2 * N), 1);
    }
    zgtsv(N, NRHS, WORK(1), WORK(N), WORK(2 * N), B, LDB, INFO);

    // 3) Backward substitution with L**H

    if (N > 1) {
      // Compute L**H \ B -> B   [ L**H \ (T \ (L \P**T * B) ) ]

      ztrsm('L', 'L', 'C', 'U', N - 1, NRHS, Complex.one, A(2, 1), LDA, B(2, 1),
          LDB);

      // Pivot, P * B  [ P * (L**H \ (T \ (L \P**T * B) )) ]

      for (K = N; K >= 1; K--) {
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
      }
    }
  }
}
