// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/lapack.dart';

void zlavhp(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final int N,
  final int NRHS,
  final Array<Complex> A_,
  final Array<int> IPIV_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final IPIV = IPIV_.having();
  final B = B_.having(ld: LDB);

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!lsame(TRANS, 'N') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (!lsame(DIAG, 'U') && !lsame(DIAG, 'N')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZLAVHP', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  final NOUNIT = lsame(DIAG, 'N');
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

      var K = 1;
      var KC = 1;
      while (K <= N) {
        // 1 x 1 pivot block

        if (IPIV[K] > 0) {
          // Multiply by the diagonal element if forming U * D.

          if (NOUNIT) zscal(NRHS, A[KC + K - 1], B(K, 1).asArray(), LDB);

          // Multiply by P(K) * inv(U(K))  if K > 1.

          if (K > 1) {
            // Apply the transformation.

            zgeru(K - 1, NRHS, Complex.one, A(KC), 1, B(K, 1).asArray(), LDB,
                B(1, 1), LDB);

            // Interchange if P(K) != I.

            final KP = IPIV[K];
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          KC += K;
          K++;
        } else {
          // 2 x 2 pivot block

          final KCNEXT = KC + K;

          // Multiply by the diagonal block if forming U * D.

          if (NOUNIT) {
            final D11 = A[KCNEXT - 1];
            final D22 = A[KCNEXT + K];
            final D12 = A[KCNEXT + K - 1];
            final D21 = D12.conjugate();
            for (var J = 1; J <= NRHS; J++) {
              final T1 = B[K][J];
              final T2 = B[K + 1][J];
              B[K][J] = D11 * T1 + D12 * T2;
              B[K + 1][J] = D21 * T1 + D22 * T2;
            }
          }

          // Multiply by  P(K) * inv(U(K))  if K > 1.

          if (K > 1) {
            // Apply the transformations.

            zgeru(K - 1, NRHS, Complex.one, A(KC), 1, B(K, 1).asArray(), LDB,
                B(1, 1), LDB);
            zgeru(K - 1, NRHS, Complex.one, A(KCNEXT), 1, B(K + 1, 1).asArray(),
                LDB, B(1, 1), LDB);

            // Interchange if P(K) != I.

            final KP = IPIV[K].abs();
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          KC = KCNEXT + K + 1;
          K += 2;
        }
      }

      // Compute  B := L*B
      // where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
    } else {
      // Loop backward applying the transformations to B.

      var K = N;
      var KC = N * (N + 1) ~/ 2 + 1;
      while (K >= 1) {
        KC -= (N - K + 1);

        // Test the pivot index.  If greater than zero, a 1 x 1
        // pivot was used, otherwise a 2 x 2 pivot was used.

        if (IPIV[K] > 0) {
          // 1 x 1 pivot block:

          // Multiply by the diagonal element if forming L * D.

          if (NOUNIT) zscal(NRHS, A[KC], B(K, 1).asArray(), LDB);

          // Multiply by  P(K) * inv(L(K))  if K < N.

          if (K != N) {
            final KP = IPIV[K];

            // Apply the transformation.

            zgeru(N - K, NRHS, Complex.one, A(KC + 1), 1, B(K, 1).asArray(),
                LDB, B(K + 1, 1), LDB);

            // Interchange if a permutation was applied at the
            // K-th step of the factorization.

            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          K--;
        } else {
          // 2 x 2 pivot block:

          final KCNEXT = KC - (N - K + 2);

          // Multiply by the diagonal block if forming L * D.

          if (NOUNIT) {
            final D11 = A[KCNEXT];
            final D22 = A[KC];
            final D21 = A[KCNEXT + 1];
            final D12 = D21.conjugate();
            for (var J = 1; J <= NRHS; J++) {
              final T1 = B[K - 1][J];
              final T2 = B[K][J];
              B[K - 1][J] = D11 * T1 + D12 * T2;
              B[K][J] = D21 * T1 + D22 * T2;
            }
          }

          // Multiply by  P(K) * inv(L(K))  if K < N.

          if (K != N) {
            // Apply the transformation.

            zgeru(N - K, NRHS, Complex.one, A(KC + 1), 1, B(K, 1).asArray(),
                LDB, B(K + 1, 1), LDB);
            zgeru(N - K, NRHS, Complex.one, A(KCNEXT + 2), 1,
                B(K - 1, 1).asArray(), LDB, B(K + 1, 1), LDB);

            // Interchange if a permutation was applied at the
            // K-th step of the factorization.

            final KP = IPIV[K].abs();
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }
          }
          KC = KCNEXT;
          K -= 2;
        }
      }
    }
    // -------------------------------------------------
    //
    // Compute  B := A^H * B  (conjugate transpose)
    //
    // -------------------------------------------------
  } else {
    // Form  B := U^H*B
    // where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
    // and   U^H = inv(U^H(1))*P(1)* ... *inv(U^H(m))*P(m)

    if (lsame(UPLO, 'U')) {
      // Loop backward applying the transformations.

      var K = N;
      var KC = N * (N + 1) ~/ 2 + 1;
      while (K >= 1) {
        KC -= K;

        // 1 x 1 pivot block.

        if (IPIV[K] > 0) {
          if (K > 1) {
            // Interchange if P(K) != I.

            final KP = IPIV[K];
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformation:
            //    y := y - B' * conjg(x)
            // where x is a column of A and y is a row of B.

            zlacgv(NRHS, B(K, 1).asArray(), LDB);
            zgemv('Conjugate', K - 1, NRHS, Complex.one, B, LDB, A(KC), 1,
                Complex.one, B(K, 1).asArray(), LDB);
            zlacgv(NRHS, B(K, 1).asArray(), LDB);
          }
          if (NOUNIT) zscal(NRHS, A[KC + K - 1], B(K, 1).asArray(), LDB);
          K--;

          // 2 x 2 pivot block.
        } else {
          final KCNEXT = KC - (K - 1);
          if (K > 2) {
            // Interchange if P(K) != I.

            final KP = IPIV[K].abs();
            if (KP != K - 1) {
              zswap(NRHS, B(K - 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformations.

            zlacgv(NRHS, B(K, 1).asArray(), LDB);
            zgemv('Conjugate', K - 2, NRHS, Complex.one, B, LDB, A(KC), 1,
                Complex.one, B(K, 1).asArray(), LDB);
            zlacgv(NRHS, B(K, 1).asArray(), LDB);

            zlacgv(NRHS, B(K - 1, 1).asArray(), LDB);
            zgemv('Conjugate', K - 2, NRHS, Complex.one, B, LDB, A(KCNEXT), 1,
                Complex.one, B(K - 1, 1).asArray(), LDB);
            zlacgv(NRHS, B(K - 1, 1).asArray(), LDB);
          }

          // Multiply by the diagonal block if non-unit.

          if (NOUNIT) {
            final D11 = A[KC - 1];
            final D22 = A[KC + K - 1];
            final D12 = A[KC + K - 2];
            final D21 = D12.conjugate();
            for (var J = 1; J <= NRHS; J++) {
              final T1 = B[K - 1][J];
              final T2 = B[K][J];
              B[K - 1][J] = D11 * T1 + D12 * T2;
              B[K][J] = D21 * T1 + D22 * T2;
            }
          }
          KC = KCNEXT;
          K -= 2;
        }
      }

      // Form  B := L^H*B
      // where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
      // and   L^H = inv(L(m))*P(m)* ... *inv(L(1))*P(1)
    } else {
      // Loop forward applying the L-transformations.

      var K = 1;
      var KC = 1;
      while (K <= N) {
        // 1 x 1 pivot block

        if (IPIV[K] > 0) {
          if (K < N) {
            // Interchange if P(K) != I.

            final KP = IPIV[K];
            if (KP != K) {
              zswap(NRHS, B(K, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformation

            zlacgv(NRHS, B(K, 1).asArray(), LDB);
            zgemv('Conjugate', N - K, NRHS, Complex.one, B(K + 1, 1), LDB,
                A(KC + 1), 1, Complex.one, B(K, 1).asArray(), LDB);
            zlacgv(NRHS, B(K, 1).asArray(), LDB);
          }
          if (NOUNIT) zscal(NRHS, A[KC], B(K, 1).asArray(), LDB);
          KC += N - K + 1;
          K++;

          // 2 x 2 pivot block.
        } else {
          final KCNEXT = KC + N - K + 1;
          if (K < N - 1) {
            // Interchange if P(K) != I.

            final KP = IPIV[K].abs();
            if (KP != K + 1) {
              zswap(NRHS, B(K + 1, 1).asArray(), LDB, B(KP, 1).asArray(), LDB);
            }

            // Apply the transformation

            zlacgv(NRHS, B(K + 1, 1).asArray(), LDB);
            zgemv('Conjugate', N - K - 1, NRHS, Complex.one, B(K + 2, 1), LDB,
                A(KCNEXT + 1), 1, Complex.one, B(K + 1, 1).asArray(), LDB);
            zlacgv(NRHS, B(K + 1, 1).asArray(), LDB);

            zlacgv(NRHS, B(K, 1).asArray(), LDB);
            zgemv('Conjugate', N - K - 1, NRHS, Complex.one, B(K + 2, 1), LDB,
                A(KC + 2), 1, Complex.one, B(K, 1).asArray(), LDB);
            zlacgv(NRHS, B(K, 1).asArray(), LDB);
          }

          // Multiply by the diagonal block if non-unit.

          if (NOUNIT) {
            final D11 = A[KC];
            final D22 = A[KCNEXT];
            final D21 = A[KC + 1];
            final D12 = D21.conjugate();
            for (var J = 1; J <= NRHS; J++) {
              final T1 = B[K][J];
              final T2 = B[K + 1][J];
              B[K][J] = D11 * T1 + D12 * T2;
              B[K + 1][J] = D21 * T1 + D22 * T2;
            }
          }
          KC = KCNEXT + (N - K);
          K += 2;
        }
      }
    }
  }
}
