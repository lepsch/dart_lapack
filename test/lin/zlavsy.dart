// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

void zlavsy(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);
  bool NOUNIT;
  int J, K, KP;
  Complex D11, D12, D21, D22, T1, T2;

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'T')) {
    INFO.value = -2;
  } else if (!lsame(DIAG, 'U') && !lsame(DIAG, 'N')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('ZLAVSY', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  NOUNIT = lsame(DIAG, 'N');
  // ------------------------------------------
  //
  // Compute  B := A * B  (No transpose)
  //
  // ------------------------------------------
  if (lsame(TRANS, 'N')) {
    // Compute  B := U*B
    // where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))

    if (lsame(UPLO, 'U')) {
      // Loop forward applying the transformations.

      K = 1;
      while (K <= N) {
        if (IPIV[K] > 0) {
          // 1 x 1 pivot block

          // Multiply by the diagonal element if forming U * D.

          if (NOUNIT) zscal(NRHS, A[K][K], B(K, 1).asArray(), LDB);

          // Multiply by  P(K) * inv(U(K))  if K > 1.

          if (K > 1) {
            // Apply the transformation.

            zgeru(K - 1, NRHS, Complex.one, A(1, K).asArray(), 1,
                B(K, 1).asArray(), LDB, B(1, 1), LDB);

            // Interchange if P(K) != I.

            KP = IPIV[K];
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          K++;
        } else {
          // 2 x 2 pivot block

          // Multiply by the diagonal block if forming U * D.

          if (NOUNIT) {
            D11 = A[K][K];
            D22 = A[K + 1][K + 1];
            D12 = A[K][K + 1];
            D21 = D12;
            for (J = 1; J <= NRHS; J++) {
              T1 = B[K][J];
              T2 = B[K + 1][J];
              B[K][J] = D11 * T1 + D12 * T2;
              B[K + 1][J] = D21 * T1 + D22 * T2;
            }
          }

          // Multiply by  P(K) * inv(U(K))  if K > 1.

          if (K > 1) {
            // Apply the transformations.

            zgeru(K - 1, NRHS, Complex.one, A(1, K).asArray(), 1,
                B(K, 1).asArray(), LDB, B(1, 1), LDB);
            zgeru(K - 1, NRHS, Complex.one, A(1, K + 1).asArray(), 1,
                B(K + 1, 1).asArray(), LDB, B(1, 1), LDB);

            // Interchange if P(K) != I.

            KP = IPIV[K].abs();
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          K += 2;
        }
      }

      // Compute  B := L*B
      // where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
    } else {
      // Loop backward applying the transformations to B.

      K = N;
      while (K >= 1) {
        // Test the pivot index.  If greater than zero, a 1 x 1
        // pivot was used, otherwise a 2 x 2 pivot was used.

        if (IPIV[K] > 0) {
          // 1 x 1 pivot block:

          // Multiply by the diagonal element if forming L * D.

          if (NOUNIT) zscal(NRHS, A[K][K], B(K, 1).asArray(), LDB);

          // Multiply by  P(K) * inv(L(K))  if K < N.

          if (K != N) {
            KP = IPIV[K];

            // Apply the transformation.

            zgeru(N - K, NRHS, Complex.one, A(K + 1, K).asArray(), 1,
                B(K, 1).asArray(), LDB, B(K + 1, 1), LDB);

            // Interchange if a permutation was applied at the
            // K-th step of the factorization.

            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          K--;
        } else {
          // 2 x 2 pivot block:

          // Multiply by the diagonal block if forming L * D.

          if (NOUNIT) {
            D11 = A[K - 1][K - 1];
            D22 = A[K][K];
            D21 = A[K][K - 1];
            D12 = D21;
            for (J = 1; J <= NRHS; J++) {
              T1 = B[K - 1][J];
              T2 = B[K][J];
              B[K - 1][J] = D11 * T1 + D12 * T2;
              B[K][J] = D21 * T1 + D22 * T2;
            }
          }

          // Multiply by  P(K) * inv(L(K))  if K < N.

          if (K != N) {
            // Apply the transformation.

            zgeru(N - K, NRHS, Complex.one, A(K + 1, K).asArray(), 1,
                B(K, 1).asArray(), LDB, B(K + 1, 1), LDB);
            zgeru(N - K, NRHS, Complex.one, A(K + 1, K - 1).asArray(), 1,
                B(K - 1, 1).asArray(), LDB, B(K + 1, 1), LDB);

            // Interchange if a permutation was applied at the
            // K-th step of the factorization.

            KP = IPIV[K].abs();
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          K -= 2;
        }
      }
    }
    // ----------------------------------------
    //
    // Compute  B := A' * B  (transpose)
    //
    // ----------------------------------------
  } else if (lsame(TRANS, 'T')) {
    // Form  B := U'*B
    // where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
    // and   U' = inv(U'(1))*P(1)* ... *inv(U'(m))*P(m)

    if (lsame(UPLO, 'U')) {
      // Loop backward applying the transformations.

      K = N;
      while (K >= 1) {
        // 1 x 1 pivot block.

        if (IPIV[K] > 0) {
          if (K > 1) {
            // Interchange if P(K) != I.

            KP = IPIV[K];
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformation

            zgemv('Transpose', K - 1, NRHS, Complex.one, B, LDB,
                A(1, K).asArray(), 1, Complex.one, B(K, 1).asArray(), LDB);
          }
          if (NOUNIT) zscal(NRHS, A[K][K], B(K, 1).asArray(), LDB);
          K--;

          // 2 x 2 pivot block.
        } else {
          if (K > 2) {
            // Interchange if P(K) != I.

            KP = IPIV[K].abs();
            if (KP != K - 1) {
              zswap(NRHS, B(K - 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformations

            zgemv('Transpose', K - 2, NRHS, Complex.one, B, LDB,
                A(1, K).asArray(), 1, Complex.one, B(K, 1).asArray(), LDB);
            zgemv(
                'Transpose',
                K - 2,
                NRHS,
                Complex.one,
                B,
                LDB,
                A(1, K - 1).asArray(),
                1,
                Complex.one,
                B(K - 1, 1).asArray(),
                LDB);
          }

          // Multiply by the diagonal block if non-unit.

          if (NOUNIT) {
            D11 = A[K - 1][K - 1];
            D22 = A[K][K];
            D12 = A[K - 1][K];
            D21 = D12;
            for (J = 1; J <= NRHS; J++) {
              T1 = B[K - 1][J];
              T2 = B[K][J];
              B[K - 1][J] = D11 * T1 + D12 * T2;
              B[K][J] = D21 * T1 + D22 * T2;
            }
          }
          K -= 2;
        }
      }

      // Form  B := L'*B
      // where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
      // and   L' = inv(L'(m))*P(m)* ... *inv(L'(1))*P(1)
    } else {
      // Loop forward applying the L-transformations.

      K = 1;
      while (K <= N) {
        // 1 x 1 pivot block

        if (IPIV[K] > 0) {
          if (K < N) {
            // Interchange if P(K) != I.

            KP = IPIV[K];
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformation

            zgemv('Transpose', N - K, NRHS, Complex.one, B(K + 1, 1), LDB,
                A(K + 1, K).asArray(), 1, Complex.one, B(K, 1).asArray(), LDB);
          }
          if (NOUNIT) zscal(NRHS, A[K][K], B(K, 1).asArray(), LDB);
          K++;

          // 2 x 2 pivot block.
        } else {
          if (K < N - 1) {
            // Interchange if P(K) != I.

            KP = IPIV[K].abs();
            if (KP != K + 1) {
              zswap(NRHS, B(K + 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformation

            zgemv(
                'Transpose',
                N - K - 1,
                NRHS,
                Complex.one,
                B(K + 2, 1),
                LDB,
                A(K + 2, K + 1).asArray(),
                1,
                Complex.one,
                B(K + 1, 1).asArray(),
                LDB);
            zgemv('Transpose', N - K - 1, NRHS, Complex.one, B(K + 2, 1), LDB,
                A(K + 2, K).asArray(), 1, Complex.one, B(K, 1).asArray(), LDB);
          }

          // Multiply by the diagonal block if non-unit.

          if (NOUNIT) {
            D11 = A[K][K];
            D22 = A[K + 1][K + 1];
            D21 = A[K + 1][K];
            D12 = D21;
            for (J = 1; J <= NRHS; J++) {
              T1 = B[K][J];
              T2 = B[K + 1][J];
              B[K][J] = D11 * T1 + D12 * T2;
              B[K + 1][J] = D21 * T1 + D22 * T2;
            }
          }
          K += 2;
        }
      }
    }
  }
}
