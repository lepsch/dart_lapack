import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zhemm.dart';
import 'package:lapack/src/blas/zher2k.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhegs2.dart';

void zhegst(
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

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  const HALF = Complex(0.5, 0.0);
  bool UPPER;
  int K, KB, NB;

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
    xerbla('ZHEGST', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Determine the block size for this environment.

  NB = ilaenv(1, 'ZHEGST', UPLO, N, -1, -1, -1);

  if (NB <= 1 || NB >= N) {
    // Use unblocked code

    zhegs2(ITYPE, UPLO, N, A, LDA, B, LDB, INFO);
  } else {
    // Use blocked code

    if (ITYPE == 1) {
      if (UPPER) {
        // Compute inv(U**H)*A*inv(U)

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the upper triangle of A(k:n,k:n)

          zhegs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
          if (K + KB <= N) {
            ztrsm('Left', UPLO, 'Conjugate transpose', 'Non-unit', KB,
                N - K - KB + 1, Complex.one, B(K, K), LDB, A(K, K + KB), LDA);
            zhemm('Left', UPLO, KB, N - K - KB + 1, -HALF, A(K, K), LDA,
                B(K, K + KB), LDB, Complex.one, A(K, K + KB), LDA);
            zher2k(
                UPLO,
                'Conjugate transpose',
                N - K - KB + 1,
                KB,
                -Complex.one,
                A(K, K + KB),
                LDA,
                B(K, K + KB),
                LDB,
                ONE,
                A(K + KB, K + KB),
                LDA);
            zhemm('Left', UPLO, KB, N - K - KB + 1, -HALF, A(K, K), LDA,
                B(K, K + KB), LDB, Complex.one, A(K, K + KB), LDA);
            ztrsm('Right', UPLO, 'No transpose', 'Non-unit', KB, N - K - KB + 1,
                Complex.one, B(K + KB, K + KB), LDB, A(K, K + KB), LDA);
          }
        }
      } else {
        // Compute inv(L)*A*inv(L**H)

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the lower triangle of A(k:n,k:n)

          zhegs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
          if (K + KB <= N) {
            ztrsm(
                'Right',
                UPLO,
                'Conjugate transpose',
                'Non-unit',
                N - K - KB + 1,
                KB,
                Complex.one,
                B(K, K),
                LDB,
                A(K + KB, K),
                LDA);
            zhemm('Right', UPLO, N - K - KB + 1, KB, -HALF, A(K, K), LDA,
                B(K + KB, K), LDB, Complex.one, A(K + KB, K), LDA);
            zher2k(
                UPLO,
                'No transpose',
                N - K - KB + 1,
                KB,
                -Complex.one,
                A(K + KB, K),
                LDA,
                B(K + KB, K),
                LDB,
                ONE,
                A(K + KB, K + KB),
                LDA);
            zhemm('Right', UPLO, N - K - KB + 1, KB, -HALF, A(K, K), LDA,
                B(K + KB, K), LDB, Complex.one, A(K + KB, K), LDA);
            ztrsm('Left', UPLO, 'No transpose', 'Non-unit', N - K - KB + 1, KB,
                Complex.one, B(K + KB, K + KB), LDB, A(K + KB, K), LDA);
          }
        }
      }
    } else {
      if (UPPER) {
        // Compute U*A*U**H

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the upper triangle of A(1:k+kb-1,1:k+kb-1)

          ztrmm('Left', UPLO, 'No transpose', 'Non-unit', K - 1, KB,
              Complex.one, B, LDB, A(1, K), LDA);
          zhemm('Right', UPLO, K - 1, KB, HALF, A(K, K), LDA, B(1, K), LDB,
              Complex.one, A(1, K), LDA);
          zher2k(UPLO, 'No transpose', K - 1, KB, Complex.one, A(1, K), LDA,
              B(1, K), LDB, ONE, A, LDA);
          zhemm('Right', UPLO, K - 1, KB, HALF, A(K, K), LDA, B(1, K), LDB,
              Complex.one, A(1, K), LDA);
          ztrmm('Right', UPLO, 'Conjugate transpose', 'Non-unit', K - 1, KB,
              Complex.one, B(K, K), LDB, A(1, K), LDA);
          zhegs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
        }
      } else {
        // Compute L**H*A*L

        for (K = 1; K <= N; K += NB) {
          KB = min(N - K + 1, NB);

          // Update the lower triangle of A(1:k+kb-1,1:k+kb-1)

          ztrmm('Right', UPLO, 'No transpose', 'Non-unit', KB, K - 1,
              Complex.one, B, LDB, A(K, 1), LDA);
          zhemm('Left', UPLO, KB, K - 1, HALF, A(K, K), LDA, B(K, 1), LDB,
              Complex.one, A(K, 1), LDA);
          zher2k(UPLO, 'Conjugate transpose', K - 1, KB, Complex.one, A(K, 1),
              LDA, B(K, 1), LDB, ONE, A, LDA);
          zhemm('Left', UPLO, KB, K - 1, HALF, A(K, K), LDA, B(K, 1), LDB,
              Complex.one, A(K, 1), LDA);
          ztrmm('Left', UPLO, 'Conjugate transpose', 'Non-unit', KB, K - 1,
              Complex.one, B(K, K), LDB, A(K, 1), LDA);
          zhegs2(ITYPE, UPLO, KB, A(K, K), LDA, B(K, K), LDB, INFO);
        }
      }
    }
  }
}
