import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/dsymv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytri_rook(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int K, KP, KSTEP;
  double AK, AKKP1, AKP1, D, T, TEMP;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DSYTRI_ROOK', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check that the diagonal matrix D is nonsingular.

  if (UPPER) {
    // Upper triangular storage: examine D from bottom to top

    for (INFO.value = N; INFO.value >= 1; INFO.value--) {
      // 10
      if (IPIV[INFO.value] > 0 && A[INFO.value][INFO.value] == ZERO) return;
    } // 10
  } else {
    // Lower triangular storage: examine D from top to bottom.

    for (INFO.value = 1; INFO.value <= N; INFO.value++) {
      // 20
      if (IPIV[INFO.value] > 0 && A[INFO.value][INFO.value] == ZERO) return;
    } // 20
  }
  INFO.value = 0;

  if (UPPER) {
    // Compute inv(A) from the factorization A = U*D*U**T.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Invert the diagonal block.

        A[K][K] = ONE / A[K][K];

        // Compute column K of the inverse.

        if (K > 1) {
          dcopy(K - 1, A(1, K).asArray(), 1, WORK, 1);
          dsymv(UPLO, K - 1, -ONE, A, LDA, WORK, 1, ZERO, A(1, K).asArray(), 1);
          A[K][K] = A[K][K] - ddot(K - 1, WORK, 1, A(1, K).asArray(), 1);
        }
        KSTEP = 1;
      } else {
        // 2 x 2 diagonal block

        // Invert the diagonal block.

        T = (A[K][K + 1]).abs();
        AK = A[K][K] / T;
        AKP1 = A[K + 1][K + 1] / T;
        AKKP1 = A[K][K + 1] / T;
        D = T * (AK * AKP1 - ONE);
        A[K][K] = AKP1 / D;
        A[K + 1][K + 1] = AK / D;
        A[K][K + 1] = -AKKP1 / D;

        // Compute columns K and K+1 of the inverse.

        if (K > 1) {
          dcopy(K - 1, A(1, K).asArray(), 1, WORK, 1);
          dsymv(UPLO, K - 1, -ONE, A, LDA, WORK, 1, ZERO, A(1, K).asArray(), 1);
          A[K][K] = A[K][K] - ddot(K - 1, WORK, 1, A(1, K).asArray(), 1);
          A[K][K + 1] = A[K][K + 1] -
              ddot(K - 1, A(1, K).asArray(), 1, A(1, K + 1).asArray(), 1);
          dcopy(K - 1, A(1, K + 1).asArray(), 1, WORK, 1);
          dsymv(UPLO, K - 1, -ONE, A, LDA, WORK, 1, ZERO, A(1, K + 1).asArray(),
              1);
          A[K + 1][K + 1] =
              A[K + 1][K + 1] - ddot(K - 1, WORK, 1, A(1, K + 1).asArray(), 1);
        }
        KSTEP = 2;
      }

      if (KSTEP == 1) {
        // Interchange rows and columns K and IPIV(K) in the leading
        // submatrix A(1:k+1,1:k+1)

        KP = IPIV[K];
        if (KP != K) {
          if (KP > 1) {
            dswap(KP - 1, A(1, K).asArray(), 1, A(1, KP).asArray(), 1);
          }
          dswap(K - KP - 1, A(KP + 1, K).asArray(), 1, A(KP, KP + 1).asArray(),
              LDA);
          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      } else {
        // Interchange rows and columns K and K+1 with -IPIV(K) and
        // -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1)

        KP = -IPIV[K];
        if (KP != K) {
          if (KP > 1) {
            dswap(KP - 1, A(1, K).asArray(), 1, A(1, KP).asArray(), 1);
          }
          dswap(K - KP - 1, A(KP + 1, K).asArray(), 1, A(KP, KP + 1).asArray(),
              LDA);

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
          TEMP = A[K][K + 1];
          A[K][K + 1] = A[KP][K + 1];
          A[KP][K + 1] = TEMP;
        }

        K = K + 1;
        KP = -IPIV[K];
        if (KP != K) {
          if (KP > 1) {
            dswap(KP - 1, A(1, K).asArray(), 1, A(1, KP).asArray(), 1);
          }
          dswap(K - KP - 1, A(KP + 1, K).asArray(), 1, A(KP, KP + 1).asArray(),
              LDA);
          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      }

      K = K + 1;
    }
  } else {
    // Compute inv(A) from the factorization A = L*D*L**T.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = N;
    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Invert the diagonal block.

        A[K][K] = ONE / A[K][K];

        // Compute column K of the inverse.

        if (K < N) {
          dcopy(N - K, A(K + 1, K).asArray(), 1, WORK, 1);
          dsymv(UPLO, N - K, -ONE, A(K + 1, K + 1), LDA, WORK, 1, ZERO,
              A(K + 1, K).asArray(), 1);
          A[K][K] = A[K][K] - ddot(N - K, WORK, 1, A(K + 1, K).asArray(), 1);
        }
        KSTEP = 1;
      } else {
        // 2 x 2 diagonal block

        // Invert the diagonal block.

        T = (A[K][K - 1]).abs();
        AK = A[K - 1][K - 1] / T;
        AKP1 = A[K][K] / T;
        AKKP1 = A[K][K - 1] / T;
        D = T * (AK * AKP1 - ONE);
        A[K - 1][K - 1] = AKP1 / D;
        A[K][K] = AK / D;
        A[K][K - 1] = -AKKP1 / D;

        // Compute columns K-1 and K of the inverse.

        if (K < N) {
          dcopy(N - K, A(K + 1, K).asArray(), 1, WORK, 1);
          dsymv(UPLO, N - K, -ONE, A(K + 1, K + 1), LDA, WORK, 1, ZERO,
              A(K + 1, K).asArray(), 1);
          A[K][K] = A[K][K] - ddot(N - K, WORK, 1, A(K + 1, K).asArray(), 1);
          A[K][K - 1] = A[K][K - 1] -
              ddot(N - K, A(K + 1, K).asArray(), 1, A(K + 1, K - 1).asArray(),
                  1);
          dcopy(N - K, A(K + 1, K - 1).asArray(), 1, WORK, 1);
          dsymv(UPLO, N - K, -ONE, A(K + 1, K + 1), LDA, WORK, 1, ZERO,
              A(K + 1, K - 1).asArray(), 1);
          A[K - 1][K - 1] = A[K - 1][K - 1] -
              ddot(N - K, WORK, 1, A(K + 1, K - 1).asArray(), 1);
        }
        KSTEP = 2;
      }

      if (KSTEP == 1) {
        // Interchange rows and columns K and IPIV(K) in the trailing
        // submatrix A(k-1:n,k-1:n)

        KP = IPIV[K];
        if (KP != K) {
          if (KP < N) {
            dswap(
                N - KP, A(KP + 1, K).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }
          dswap(KP - K - 1, A(K + 1, K).asArray(), 1, A(KP, K + 1).asArray(),
              LDA);
          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      } else {
        // Interchange rows and columns K and K-1 with -IPIV(K) and
        // -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)

        KP = -IPIV[K];
        if (KP != K) {
          if (KP < N) {
            dswap(
                N - KP, A(KP + 1, K).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }
          dswap(KP - K - 1, A(K + 1, K).asArray(), 1, A(KP, K + 1).asArray(),
              LDA);

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
          TEMP = A[K][K - 1];
          A[K][K - 1] = A[KP][K - 1];
          A[KP][K - 1] = TEMP;
        }

        K = K - 1;
        KP = -IPIV[K];
        if (KP != K) {
          if (KP < N) {
            dswap(
                N - KP, A(KP + 1, K).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }
          dswap(KP - K - 1, A(K + 1, K).asArray(), 1, A(KP, K + 1).asArray(),
              LDA);
          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      }

      K = K - 1;
    }
  }
}
