import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zsyconv.dart';

void zhetrs2(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool UPPER;
  int I, J, K, KP;
  double S;
  Complex AK, AKM1, AKM1K, BK, BKM1, DENOM;
  final IINFO = Box(0);

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
    xerbla('ZHETRS2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  // Convert A

  zsyconv(UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO);

  if (UPPER) {
    // Solve A*X = B, where A = U*D*U**H.

    // P**T * B
    K = N;
    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block
        // Interchange rows K and IPIV[K].
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K--;
      } else {
        // 2 x 2 diagonal block
        // Interchange rows K-1 and -IPIV[K].
        KP = -IPIV[K];
        if (KP == -IPIV[K - 1]) {
          zswap(NRHS, B(K - 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K -= 2;
      }
    }

    // Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

    ztrsm('L', 'U', 'N', 'U', N, NRHS, Complex.one, A, LDA, B, LDB);

    // Compute D \ B -> B   [ D \ (U \P**T * B) ]

    I = N;
    while (I >= 1) {
      if (IPIV[I] > 0) {
        S = 1.0 / A[I][I].toDouble();
        zdscal(NRHS, S, B(I, 1).asArray(), LDB);
      } else if (I > 1) {
        if (IPIV[I - 1] == IPIV[I]) {
          AKM1K = WORK[I];
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
      }
      I--;
    }

    // Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ]

    ztrsm('L', 'U', 'C', 'U', N, NRHS, Complex.one, A, LDA, B, LDB);

    // P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ]

    K = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block
        // Interchange rows K and IPIV[K].
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K++;
      } else {
        // 2 x 2 diagonal block
        // Interchange rows K-1 and -IPIV[K].
        KP = -IPIV[K];
        if (K < N && KP == -IPIV[K + 1]) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K += 2;
      }
    }
  } else {
    // Solve A*X = B, where A = L*D*L**H.

    // P**T * B
    K = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block
        // Interchange rows K and IPIV[K].
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K++;
      } else {
        // 2 x 2 diagonal block
        // Interchange rows K and -IPIV(K+1).
        KP = -IPIV[K + 1];
        if (KP == -IPIV[K]) {
          zswap(NRHS, B(K + 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K += 2;
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
      } else {
        AKM1K = WORK[I];
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

    K = N;
    while (K >= 1) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block
        // Interchange rows K and IPIV[K].
        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K--;
      } else {
        // 2 x 2 diagonal block
        // Interchange rows K-1 and -IPIV[K].
        KP = -IPIV[K];
        if (K > 1 && KP == -IPIV[K - 1]) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K -= 2;
      }
    }
  }

  // Revert A

  zsyconv(UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO);
}
