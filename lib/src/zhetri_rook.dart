import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zhemv.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zhetri_rook(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final IPIV = IPIV_.dim();
  final WORK = WORK_.dim();
  const ONE = 1.0;
  bool UPPER;
  int J, K, KP, KSTEP;
  double AK, AKP1, D, T;
  Complex AKKP1, TEMP;

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
    xerbla('ZHETRI_ROOK', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check that the diagonal matrix D is nonsingular.

  if (UPPER) {
    // Upper triangular storage: examine D from bottom to top

    for (INFO.value = N; INFO.value >= 1; INFO.value--) {
      // 10
      if (IPIV[INFO.value] > 0 && A[INFO.value][INFO.value] == Complex.zero) {
        return;
      }
    } // 10
  } else {
    // Lower triangular storage: examine D from top to bottom.

    for (INFO.value = 1; INFO.value <= N; INFO.value++) {
      // 20
      if (IPIV[INFO.value] > 0 && A[INFO.value][INFO.value] == Complex.zero) {
        return;
      }
    } // 20
  }
  INFO.value = 0;

  if (UPPER) {
    // Compute inv(A) from the factorization A = U*D*U**H.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Invert the diagonal block.

        A[K][K] = (ONE / (A[K][K]).toDouble()).toComplex();

        // Compute column K of the inverse.

        if (K > 1) {
          zcopy(K - 1, A(1, K).asArray(), 1, WORK, 1);
          zhemv(UPLO, K - 1, -Complex.one, A, LDA, WORK, 1, Complex.zero,
              A(1, K).asArray(), 1);
          A[K][K] = A[K][K] -
              zdotc(K - 1, WORK, 1, A(1, K).asArray(), 1).real.toComplex();
        }
        KSTEP = 1;
      } else {
        // 2 x 2 diagonal block

        // Invert the diagonal block.

        T = (A[K][K + 1]).abs();
        AK = (A[K][K]).toDouble() / T;
        AKP1 = (A[K + 1][K + 1]).toDouble() / T;
        AKKP1 = A[K][K + 1] / T.toComplex();
        D = T * (AK * AKP1 - ONE);
        A[K][K] = (AKP1 / D).toComplex();
        A[K + 1][K + 1] = (AK / D).toComplex();
        A[K][K + 1] = -AKKP1 / D.toComplex();

        // Compute columns K and K+1 of the inverse.

        if (K > 1) {
          zcopy(K - 1, A(1, K).asArray(), 1, WORK, 1);
          zhemv(UPLO, K - 1, -Complex.one, A, LDA, WORK, 1, Complex.zero,
              A(1, K).asArray(), 1);
          A[K][K] = A[K][K] -
              zdotc(K - 1, WORK, 1, A(1, K).asArray(), 1).real.toComplex();
          A[K][K + 1] = A[K][K + 1] -
              zdotc(K - 1, A(1, K).asArray(), 1, A(1, K + 1).asArray(), 1);
          zcopy(K - 1, A(1, K + 1).asArray(), 1, WORK, 1);
          zhemv(UPLO, K - 1, -Complex.one, A, LDA, WORK, 1, Complex.zero,
              A(1, K + 1).asArray(), 1);
          A[K + 1][K + 1] = A[K + 1][K + 1] -
              zdotc(K - 1, WORK, 1, A(1, K + 1).asArray(), 1).real.toComplex();
        }
        KSTEP = 2;
      }

      if (KSTEP == 1) {
        // Interchange rows and columns K and IPIV(K) in the leading
        // submatrix A(1:k,1:k)

        KP = IPIV[K];
        if (KP != K) {
          if (KP > 1) {
            zswap(KP - 1, A(1, K).asArray(), 1, A(1, KP).asArray(), 1);
          }

          for (J = KP + 1; J <= K - 1; J++) {
            // 40
            TEMP = A[J][K].conjugate();
            A[J][K] = A[KP][J].conjugate();
            A[KP][J] = TEMP;
          } // 40

          A[KP][K] = A[KP][K].conjugate();

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      } else {
        // Interchange rows and columns K and K+1 with -IPIV(K) and
        // -IPIV(K+1) in the leading submatrix A(k+1:n,k+1:n)

        // (1) Interchange rows and columns K and -IPIV(K)

        KP = -IPIV[K];
        if (KP != K) {
          if (KP > 1) {
            zswap(KP - 1, A(1, K).asArray(), 1, A(1, KP).asArray(), 1);
          }

          for (J = KP + 1; J <= K - 1; J++) {
            // 50
            TEMP = A[J][K].conjugate();
            A[J][K] = A[KP][J].conjugate();
            A[KP][J] = TEMP;
          } // 50

          A[KP][K] = A[KP][K].conjugate();

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;

          TEMP = A[K][K + 1];
          A[K][K + 1] = A[KP][K + 1];
          A[KP][K + 1] = TEMP;
        }

        // (2) Interchange rows and columns K+1 and -IPIV(K+1)

        K = K + 1;
        KP = -IPIV[K];
        if (KP != K) {
          if (KP > 1) {
            zswap(KP - 1, A(1, K).asArray(), 1, A(1, KP).asArray(), 1);
          }

          for (J = KP + 1; J <= K - 1; J++) {
            // 60
            TEMP = A[J][K].conjugate();
            A[J][K] = A[KP][J].conjugate();
            A[KP][J] = TEMP;
          } // 60

          A[KP][K] = A[KP][K].conjugate();

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      }

      K = K + 1;
    } // 70
  } else {
    // Compute inv(A) from the factorization A = L*D*L**H.

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = N;
    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Invert the diagonal block.

        A[K][K] = (ONE / (A[K][K]).toDouble()).toComplex();

        // Compute column K of the inverse.

        if (K < N) {
          zcopy(N - K, A(K + 1, K).asArray(), 1, WORK, 1);
          zhemv(UPLO, N - K, -Complex.one, A(K + 1, K + 1), LDA, WORK, 1,
              Complex.zero, A(K + 1, K).asArray(), 1);
          A[K][K] = A[K][K] -
              zdotc(N - K, WORK, 1, A(K + 1, K).asArray(), 1).real.toComplex();
        }
        KSTEP = 1;
      } else {
        // 2 x 2 diagonal block

        // Invert the diagonal block.

        T = (A[K][K - 1]).abs();
        AK = (A[K - 1][K - 1]).toDouble() / T;
        AKP1 = (A[K][K]).toDouble() / T;
        AKKP1 = A[K][K - 1] / T.toComplex();
        D = T * (AK * AKP1 - ONE);
        A[K - 1][K - 1] = (AKP1 / D).toComplex();
        A[K][K] = (AK / D).toComplex();
        A[K][K - 1] = -AKKP1 / D.toComplex();

        // Compute columns K-1 and K of the inverse.

        if (K < N) {
          zcopy(N - K, A(K + 1, K).asArray(), 1, WORK, 1);
          zhemv(UPLO, N - K, -Complex.one, A(K + 1, K + 1), LDA, WORK, 1,
              Complex.zero, A(K + 1, K).asArray(), 1);
          A[K][K] = A[K][K] -
              zdotc(N - K, WORK, 1, A(K + 1, K).asArray(), 1).real.toComplex();
          A[K][K - 1] = A[K][K - 1] -
              zdotc(N - K, A(K + 1, K).asArray(), 1, A(K + 1, K - 1).asArray(),
                  1);
          zcopy(N - K, A(K + 1, K - 1).asArray(), 1, WORK, 1);
          zhemv(UPLO, N - K, -Complex.one, A(K + 1, K + 1), LDA, WORK, 1,
              Complex.zero, A(K + 1, K - 1).asArray(), 1);
          A[K - 1][K - 1] = A[K - 1][K - 1] -
              zdotc(N - K, WORK, 1, A(K + 1, K - 1).asArray(), 1)
                  .real
                  .toComplex();
        }
        KSTEP = 2;
      }

      if (KSTEP == 1) {
        // Interchange rows and columns K and IPIV(K) in the trailing
        // submatrix A(k:n,k:n)

        KP = IPIV[K];
        if (KP != K) {
          if (KP < N) {
            zswap(
                N - KP, A(KP + 1, K).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }

          for (J = K + 1; J <= KP - 1; J++) {
            // 90
            TEMP = A[J][K].conjugate();
            A[J][K] = A[KP][J].conjugate();
            A[KP][J] = TEMP;
          } // 90

          A[KP][K] = A[KP][K].conjugate();

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      } else {
        // Interchange rows and columns K and K-1 with -IPIV(K) and
        // -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)

        // (1) Interchange rows and columns K and -IPIV(K)

        KP = -IPIV[K];
        if (KP != K) {
          if (KP < N) {
            zswap(
                N - KP, A(KP + 1, K).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }

          for (J = K + 1; J <= KP - 1; J++) {
            // 100
            TEMP = A[J][K].conjugate();
            A[J][K] = A[KP][J].conjugate();
            A[KP][J] = TEMP;
          } // 100

          A[KP][K] = A[KP][K].conjugate();

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;

          TEMP = A[K][K - 1];
          A[K][K - 1] = A[KP][K - 1];
          A[KP][K - 1] = TEMP;
        }

        // (2) Interchange rows and columns K-1 and -IPIV(K-1)

        K = K - 1;
        KP = -IPIV[K];
        if (KP != K) {
          if (KP < N) {
            zswap(
                N - KP, A(KP + 1, K).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }

          for (J = K + 1; J <= KP - 1; J++) {
            // 110
            TEMP = A[J][K].conjugate();
            A[J][K] = A[KP][J].conjugate();
            A[KP][J] = TEMP;
          } // 110

          A[KP][K] = A[KP][K].conjugate();

          TEMP = A[K][K];
          A[K][K] = A[KP][KP];
          A[KP][KP] = TEMP;
        }
      }

      K = K - 1;
    } // 120
  }
}
