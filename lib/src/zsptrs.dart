import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgeru.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zsptrs(
  final String UPLO,
  final int N,
  final int NRHS,
  final Array<Complex> AP_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final IPIV = IPIV_.having();
  final AP = AP_.having();
  final B = B_.having(ld: LDB);
  bool UPPER;
  int J, K, KC, KP;
  Complex AK, AKM1, AKM1K, BK, BKM1, DENOM;

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZSPTRS', -INFO.value);
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
    KC = N * (N + 1) ~/ 2 + 1;
    while (K >= 1) {
      KC = KC - K;
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(U(K)), where U(K) is the transformation
        // stored in column K of A.

        zgeru(K - 1, NRHS, -Complex.one, AP(KC), 1, B(K, 1).asArray(), LDB,
            B(1, 1), LDB);

        // Multiply by the inverse of the diagonal block.

        zscal(NRHS, Complex.one / AP[KC + K - 1], B(K, 1).asArray(), LDB);
        K--;
      } else {
        // 2 x 2 diagonal block

        // Interchange rows K-1 and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K - 1) {
          zswap(NRHS, B(K - 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(U(K)), where U(K) is the transformation
        // stored in columns K-1 and K of A.

        zgeru(K - 2, NRHS, -Complex.one, AP(KC), 1, B(K, 1).asArray(), LDB,
            B(1, 1), LDB);
        zgeru(K - 2, NRHS, -Complex.one, AP(KC - (K - 1)), 1,
            B(K - 1, 1).asArray(), LDB, B(1, 1), LDB);

        // Multiply by the inverse of the diagonal block.

        AKM1K = AP[KC + K - 2];
        AKM1 = AP[KC - 1] / AKM1K;
        AK = AP[KC + K - 1] / AKM1K;
        DENOM = AKM1 * AK - Complex.one;
        for (J = 1; J <= NRHS; J++) {
          // 20
          BKM1 = B[K - 1][J] / AKM1K;
          BK = B[K][J] / AKM1K;
          B[K - 1][J] = (AK * BKM1 - BK) / DENOM;
          B[K][J] = (AKM1 * BK - BKM1) / DENOM;
        } // 20
        KC = KC - K + 1;
        K = K - 2;
      }
    }

    // Next solve U**T*X = B, overwriting B with X.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;
    KC = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Multiply by inv(U**T(K)), where U(K) is the transformation
        // stored in column K of A.

        zgemv('Transpose', K - 1, NRHS, -Complex.one, B, LDB, AP(KC), 1,
            Complex.one, B(K, 1).asArray(), LDB);

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        KC = KC + K;
        K++;
      } else {
        // 2 x 2 diagonal block

        // Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
        // stored in columns K and K+1 of A.

        zgemv('Transpose', K - 1, NRHS, -Complex.one, B, LDB, AP(KC), 1,
            Complex.one, B(K, 1).asArray(), LDB);
        zgemv('Transpose', K - 1, NRHS, -Complex.one, B, LDB, AP(KC + K), 1,
            Complex.one, B(K + 1, 1).asArray(), LDB);

        // Interchange rows K and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        KC = KC + 2 * K + 1;
        K = K + 2;
      }
    }
  } else {
    // Solve A*X = B, where A = L*D*L**T.

    // First solve L*D*X = B, overwriting B with X.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;
    KC = 1;
    while (K <= N) {
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(L(K)), where L(K) is the transformation
        // stored in column K of A.

        if (K < N) {
          zgeru(N - K, NRHS, -Complex.one, AP(KC + 1), 1, B(K, 1).asArray(),
              LDB, B(K + 1, 1), LDB);
        }

        // Multiply by the inverse of the diagonal block.

        zscal(NRHS, Complex.one / AP[KC], B(K, 1).asArray(), LDB);
        KC = KC + N - K + 1;
        K++;
      } else {
        // 2 x 2 diagonal block

        // Interchange rows K+1 and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K + 1) {
          zswap(NRHS, B(K + 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }

        // Multiply by inv(L(K)), where L(K) is the transformation
        // stored in columns K and K+1 of A.

        if (K < N - 1) {
          zgeru(N - K - 1, NRHS, -Complex.one, AP(KC + 2), 1, B(K, 1).asArray(),
              LDB, B(K + 2, 1), LDB);
          zgeru(N - K - 1, NRHS, -Complex.one, AP(KC + N - K + 2), 1,
              B(K + 1, 1).asArray(), LDB, B(K + 2, 1), LDB);
        }

        // Multiply by the inverse of the diagonal block.

        AKM1K = AP[KC + 1];
        AKM1 = AP[KC] / AKM1K;
        AK = AP[KC + N - K + 1] / AKM1K;
        DENOM = AKM1 * AK - Complex.one;
        for (J = 1; J <= NRHS; J++) {
          // 70
          BKM1 = B[K][J] / AKM1K;
          BK = B[K + 1][J] / AKM1K;
          B[K][J] = (AK * BKM1 - BK) / DENOM;
          B[K + 1][J] = (AKM1 * BK - BKM1) / DENOM;
        } // 70
        KC = KC + 2 * (N - K) + 1;
        K = K + 2;
      }
    } // 80

    // Next solve L**T*X = B, overwriting B with X.

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = N;
    KC = N * (N + 1) ~/ 2 + 1;
    while (K >= 1) {
      KC = KC - (N - K + 1);
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Multiply by inv(L**T(K)), where L(K) is the transformation
        // stored in column K of A.

        if (K < N) {
          zgemv('Transpose', N - K, NRHS, -Complex.one, B(K + 1, 1), LDB,
              AP(KC + 1), 1, Complex.one, B(K, 1).asArray(), LDB);
        }

        // Interchange rows K and IPIV(K).

        KP = IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        K--;
      } else {
        // 2 x 2 diagonal block

        // Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
        // stored in columns K-1 and K of A.

        if (K < N) {
          zgemv('Transpose', N - K, NRHS, -Complex.one, B(K + 1, 1), LDB,
              AP(KC + 1), 1, Complex.one, B(K, 1).asArray(), LDB);
          zgemv('Transpose', N - K, NRHS, -Complex.one, B(K + 1, 1), LDB,
              AP(KC - (N - K)), 1, Complex.one, B(K - 1, 1).asArray(), LDB);
        }

        // Interchange rows K and -IPIV(K).

        KP = -IPIV[K];
        if (KP != K) {
          zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
        }
        KC = KC - (N - K + 2);
        K = K - 2;
      }
    } // 100
  }
}
