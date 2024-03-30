import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsyr2.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/blas/dtrsv.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsygs2(
  final int ITYPE,
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, HALF = 0.5;
  bool UPPER;
  int K;
  double AKK, BKK, CT;

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
    xerbla('DSYGS2', -INFO.value);
    return;
  }

  if (ITYPE == 1) {
    if (UPPER) {
      // Compute inv(U**T)*A*inv(U)

      for (K = 1; K <= N; K++) {
        // Update the upper triangle of A(k:n,k:n)

        AKK = A[K][K];
        BKK = B[K][K];
        AKK /= pow(BKK, 2);
        A[K][K] = AKK;
        if (K < N) {
          dscal(N - K, ONE / BKK, A(K, K + 1).asArray(), LDA);
          CT = -HALF * AKK;
          daxpy(N - K, CT, B(K, K + 1).asArray(), LDB, A(K, K + 1).asArray(),
              LDA);
          dsyr2(UPLO, N - K, -ONE, A(K, K + 1).asArray(), LDA,
              B(K, K + 1).asArray(), LDB, A(K + 1, K + 1), LDA);
          daxpy(N - K, CT, B(K, K + 1).asArray(), LDB, A(K, K + 1).asArray(),
              LDA);
          dtrsv(UPLO, 'Transpose', 'Non-unit', N - K, B(K + 1, K + 1), LDB,
              A(K, K + 1).asArray(), LDA);
        }
      }
    } else {
      // Compute inv(L)*A*inv(L**T)

      for (K = 1; K <= N; K++) {
        // Update the lower triangle of A(k:n,k:n)

        AKK = A[K][K];
        BKK = B[K][K];
        AKK /= pow(BKK, 2);
        A[K][K] = AKK;
        if (K < N) {
          dscal(N - K, ONE / BKK, A(K + 1, K).asArray(), 1);
          CT = -HALF * AKK;
          daxpy(N - K, CT, B(K + 1, K).asArray(), 1, A(K + 1, K).asArray(), 1);
          dsyr2(UPLO, N - K, -ONE, A(K + 1, K).asArray(), 1,
              B(K + 1, K).asArray(), 1, A(K + 1, K + 1), LDA);
          daxpy(N - K, CT, B(K + 1, K).asArray(), 1, A(K + 1, K).asArray(), 1);
          dtrsv(UPLO, 'No transpose', 'Non-unit', N - K, B(K + 1, K + 1), LDB,
              A(K + 1, K).asArray(), 1);
        }
      }
    }
  } else {
    if (UPPER) {
      // Compute U*A*U**T

      for (K = 1; K <= N; K++) {
        // Update the upper triangle of A(1:k,1:k)

        AKK = A[K][K];
        BKK = B[K][K];
        dtrmv(UPLO, 'No transpose', 'Non-unit', K - 1, B, LDB,
            A(1, K).asArray(), 1);
        CT = HALF * AKK;
        daxpy(K - 1, CT, B(1, K).asArray(), 1, A(1, K).asArray(), 1);
        dsyr2(UPLO, K - 1, ONE, A(1, K).asArray(), 1, B(1, K).asArray(), 1, A,
            LDA);
        daxpy(K - 1, CT, B(1, K).asArray(), 1, A(1, K).asArray(), 1);
        dscal(K - 1, BKK, A(1, K).asArray(), 1);
        A[K][K] = AKK * pow(BKK, 2);
      }
    } else {
      // Compute L**T *A*L

      for (K = 1; K <= N; K++) {
        // Update the lower triangle of A(1:k,1:k)

        AKK = A[K][K];
        BKK = B[K][K];
        dtrmv(UPLO, 'Transpose', 'Non-unit', K - 1, B, LDB, A(K, 1).asArray(),
            LDA);
        CT = HALF * AKK;
        daxpy(K - 1, CT, B(K, 1).asArray(), LDB, A(K, 1).asArray(), LDA);
        dsyr2(UPLO, K - 1, ONE, A(K, 1).asArray(), LDA, B(K, 1).asArray(), LDB,
            A, LDA);
        daxpy(K - 1, CT, B(K, 1).asArray(), LDB, A(K, 1).asArray(), LDA);
        dscal(K - 1, BKK, A(K, 1).asArray(), LDA);
        A[K][K] = AKK * pow(BKK, 2);
      }
    }
  }
}
