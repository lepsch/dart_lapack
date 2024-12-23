// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zher2.dart';
import 'package:dart_lapack/src/blas/ztrmv.dart';
import 'package:dart_lapack/src/blas/ztrsv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';

void zhegs2(
  final int ITYPE,
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, HALF = 0.5;
  bool UPPER;
  int K;
  double AKK, BKK;
  Complex CT;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZHEGS2', -INFO.value);
    return;
  }

  if (ITYPE == 1) {
    if (UPPER) {
      // Compute inv(U**H)*A*inv(U)

      for (K = 1; K <= N; K++) {
        // Update the upper triangle of A(k:n,k:n)

        AKK = A[K][K].real;
        BKK = B[K][K].real;
        AKK /= pow(BKK, 2);
        A[K][K] = AKK.toComplex();
        if (K < N) {
          zdscal(N - K, ONE / BKK, A(K, K + 1).asArray(), LDA);
          CT = (-HALF * AKK).toComplex();
          zlacgv(N - K, A(K, K + 1).asArray(), LDA);
          zlacgv(N - K, B(K, K + 1).asArray(), LDB);
          zaxpy(N - K, CT, B(K, K + 1).asArray(), LDB, A(K, K + 1).asArray(),
              LDA);
          zher2(UPLO, N - K, -Complex.one, A(K, K + 1).asArray(), LDA,
              B(K, K + 1).asArray(), LDB, A(K + 1, K + 1), LDA);
          zaxpy(N - K, CT, B(K, K + 1).asArray(), LDB, A(K, K + 1).asArray(),
              LDA);
          zlacgv(N - K, B(K, K + 1).asArray(), LDB);
          ztrsv(UPLO, 'Conjugate transpose', 'Non-unit', N - K, B(K + 1, K + 1),
              LDB, A(K, K + 1).asArray(), LDA);
          zlacgv(N - K, A(K, K + 1).asArray(), LDA);
        }
      }
    } else {
      // Compute inv(L)*A*inv(L**H)

      for (K = 1; K <= N; K++) {
        // Update the lower triangle of A(k:n,k:n)

        AKK = A[K][K].real;
        BKK = B[K][K].real;
        AKK /= pow(BKK, 2);
        A[K][K] = AKK.toComplex();
        if (K < N) {
          zdscal(N - K, ONE / BKK, A(K + 1, K).asArray(), 1);
          CT = (-HALF * AKK).toComplex();
          zaxpy(N - K, CT, B(K + 1, K).asArray(), 1, A(K + 1, K).asArray(), 1);
          zher2(UPLO, N - K, -Complex.one, A(K + 1, K).asArray(), 1,
              B(K + 1, K).asArray(), 1, A(K + 1, K + 1), LDA);
          zaxpy(N - K, CT, B(K + 1, K).asArray(), 1, A(K + 1, K).asArray(), 1);
          ztrsv(UPLO, 'No transpose', 'Non-unit', N - K, B(K + 1, K + 1), LDB,
              A(K + 1, K).asArray(), 1);
        }
      }
    }
  } else {
    if (UPPER) {
      // Compute U*A*U**H

      for (K = 1; K <= N; K++) {
        // Update the upper triangle of A(1:k,1:k)

        AKK = A[K][K].real;
        BKK = B[K][K].real;
        ztrmv(UPLO, 'No transpose', 'Non-unit', K - 1, B, LDB,
            A(1, K).asArray(), 1);
        CT = (HALF * AKK).toComplex();
        zaxpy(K - 1, CT, B(1, K).asArray(), 1, A(1, K).asArray(), 1);
        zher2(UPLO, K - 1, Complex.one, A(1, K).asArray(), 1, B(1, K).asArray(),
            1, A, LDA);
        zaxpy(K - 1, CT, B(1, K).asArray(), 1, A(1, K).asArray(), 1);
        zdscal(K - 1, BKK, A(1, K).asArray(), 1);
        A[K][K] = (AKK * pow(BKK, 2)).toComplex();
      }
    } else {
      // Compute L**H *A*L

      for (K = 1; K <= N; K++) {
        // Update the lower triangle of A(1:k,1:k)

        AKK = A[K][K].real;
        BKK = B[K][K].real;
        zlacgv(K - 1, A(K, 1).asArray(), LDA);
        ztrmv(UPLO, 'Conjugate transpose', 'Non-unit', K - 1, B, LDB,
            A(K, 1).asArray(), LDA);
        CT = (HALF * AKK).toComplex();
        zlacgv(K - 1, B(K, 1).asArray(), LDB);
        zaxpy(K - 1, CT, B(K, 1).asArray(), LDB, A(K, 1).asArray(), LDA);
        zher2(UPLO, K - 1, Complex.one, A(K, 1).asArray(), LDA,
            B(K, 1).asArray(), LDB, A, LDA);
        zaxpy(K - 1, CT, B(K, 1).asArray(), LDB, A(K, 1).asArray(), LDA);
        zlacgv(K - 1, B(K, 1).asArray(), LDB);
        zdscal(K - 1, BKK, A(K, 1).asArray(), LDA);
        zlacgv(K - 1, A(K, 1).asArray(), LDA);
        A[K][K] = (AKK * pow(BKK, 2)).toComplex();
      }
    }
  }
}
