// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/blas/dsyr.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsytf2(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  const ZERO = 0.0, ONE = 1.0;
  const EIGHT = 8.0, SEVTEN = 17.0;
  bool UPPER;
  int I, IMAX = 0, J, JMAX, K, KK, KP, KSTEP;
  double ABSAKK,
      ALPHA,
      COLMAX,
      D11,
      D12,
      D21,
      D22,
      R1,
      ROWMAX,
      T,
      WK,
      WKM1,
      WKP1;

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
    xerbla('DSYTF2', -INFO.value);
    return;
  }

  // Initialize ALPHA for use in choosing pivot block size.

  ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  if (UPPER) {
    // Factorize A as U*D*U**T using the upper triangle of A

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2

    K = N;
    while (K >= 1) {
      KSTEP = 1;

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used

      ABSAKK = A[K][K].abs();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value.
      // Determine both COLMAX and IMAX.

      if (K > 1) {
        IMAX = idamax(K - 1, A(1, K).asArray(), 1);
        COLMAX = A[IMAX][K].abs();
      } else {
        COLMAX = ZERO;
      }

      if ((max(ABSAKK, COLMAX) == ZERO) || disnan(ABSAKK)) {
        // Column K is zero or underflow, or contains a NaN:
        // set INFO and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
      } else {
        if (ABSAKK >= ALPHA * COLMAX) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          // JMAX is the column-index of the largest off-diagonal
          // element in row IMAX, and ROWMAX is its absolute value

          JMAX = IMAX + idamax(K - IMAX, A(IMAX, IMAX + 1).asArray(), LDA);
          ROWMAX = A[IMAX][JMAX].abs();
          if (IMAX > 1) {
            JMAX = idamax(IMAX - 1, A(1, IMAX).asArray(), 1);
            ROWMAX = max(ROWMAX, A[JMAX][IMAX].abs());
          }

          if (ABSAKK >= ALPHA * COLMAX * (COLMAX / ROWMAX)) {
            // no interchange, use 1-by-1 pivot block

            KP = K;
          } else if (A[IMAX][IMAX].abs() >= ALPHA * ROWMAX) {
            // interchange rows and columns K and IMAX, use 1-by-1
            // pivot block

            KP = IMAX;
          } else {
            // interchange rows and columns K-1 and IMAX, use 2-by-2
            // pivot block

            KP = IMAX;
            KSTEP = 2;
          }
        }

        KK = K - KSTEP + 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the leading
          // submatrix A(1:k,1:k)

          dswap(KP - 1, A(1, KK).asArray(), 1, A(1, KP).asArray(), 1);
          dswap(KK - KP - 1, A(KP + 1, KK).asArray(), 1,
              A(KP, KP + 1).asArray(), LDA);
          T = A[KK][KK];
          A[KK][KK] = A[KP][KP];
          A[KP][KP] = T;
          if (KSTEP == 2) {
            T = A[K - 1][K];
            A[K - 1][K] = A[KP][K];
            A[KP][K] = T;
          }
        }

        // Update the leading submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k now holds

          // W(k) = U(k)*D(k)

          // where U(k) is the k-th column of U

          // Perform a rank-1 update of A(1:k-1,1:k-1) as

          // A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T

          R1 = ONE / A[K][K];
          dsyr(UPLO, K - 1, -R1, A(1, K).asArray(), 1, A, LDA);

          // Store U(k) in column k

          dscal(K - 1, R1, A(1, K).asArray(), 1);
        } else {
          // 2-by-2 pivot block D(k): columns k and k-1 now hold

          // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

          // where U(k) and U(k-1) are the k-th and (k-1)-th columns
          // of U

          // Perform a rank-2 update of A(1:k-2,1:k-2) as

          // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
          //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T

          if (K > 2) {
            D12 = A[K - 1][K];
            D22 = A[K - 1][K - 1] / D12;
            D11 = A[K][K] / D12;
            T = ONE / (D11 * D22 - ONE);
            D12 = T / D12;

            for (J = K - 2; J >= 1; J--) {
              WKM1 = D12 * (D11 * A[J][K - 1] - A[J][K]);
              WK = D12 * (D22 * A[J][K] - A[J][K - 1]);
              for (I = J; I >= 1; I--) {
                A[I][J] -= A[I][K] * WK + A[I][K - 1] * WKM1;
              }
              A[J][K] = WK;
              A[J][K - 1] = WKM1;
            }
          }
        }
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -KP;
        IPIV[K - 1] = -KP;
      }

      // Decrease K and return to the start of the main loop

      K -= KSTEP;
    }
  } else {
    // Factorize A as L*D*L**T using the lower triangle of A

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2

    K = 1;
    while (K <= N) {
      KSTEP = 1;

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used

      ABSAKK = A[K][K].abs();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value.
      // Determine both COLMAX and IMAX.

      if (K < N) {
        IMAX = K + idamax(N - K, A(K + 1, K).asArray(), 1);
        COLMAX = A[IMAX][K].abs();
      } else {
        COLMAX = ZERO;
      }

      if ((max(ABSAKK, COLMAX) == ZERO) || disnan(ABSAKK)) {
        // Column K is zero or underflow, or contains a NaN:
        // set INFO and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
      } else {
        if (ABSAKK >= ALPHA * COLMAX) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          // JMAX is the column-index of the largest off-diagonal
          // element in row IMAX, and ROWMAX is its absolute value

          JMAX = K - 1 + idamax(IMAX - K, A(IMAX, K).asArray(), LDA);
          ROWMAX = A[IMAX][JMAX].abs();
          if (IMAX < N) {
            JMAX = IMAX + idamax(N - IMAX, A(IMAX + 1, IMAX).asArray(), 1);
            ROWMAX = max(ROWMAX, A[JMAX][IMAX].abs());
          }

          if (ABSAKK >= ALPHA * COLMAX * (COLMAX / ROWMAX)) {
            // no interchange, use 1-by-1 pivot block

            KP = K;
          } else if (A[IMAX][IMAX].abs() >= ALPHA * ROWMAX) {
            // interchange rows and columns K and IMAX, use 1-by-1
            // pivot block

            KP = IMAX;
          } else {
            // interchange rows and columns K+1 and IMAX, use 2-by-2
            // pivot block

            KP = IMAX;
            KSTEP = 2;
          }
        }

        KK = K + KSTEP - 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the trailing
          // submatrix A(k:n,k:n)

          if (KP < N) {
            dswap(
                N - KP, A(KP + 1, KK).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }
          dswap(KP - KK - 1, A(KK + 1, KK).asArray(), 1,
              A(KP, KK + 1).asArray(), LDA);
          T = A[KK][KK];
          A[KK][KK] = A[KP][KP];
          A[KP][KP] = T;
          if (KSTEP == 2) {
            T = A[K + 1][K];
            A[K + 1][K] = A[KP][K];
            A[KP][K] = T;
          }
        }

        // Update the trailing submatrix

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k now holds

          // W(k) = L(k)*D(k)

          // where L(k) is the k-th column of L

          if (K < N) {
            // Perform a rank-1 update of A(k+1:n,k+1:n) as

            // A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T

            D11 = ONE / A[K][K];
            dsyr(UPLO, N - K, -D11, A(K + 1, K).asArray(), 1, A(K + 1, K + 1),
                LDA);

            // Store L(k) in column K

            dscal(N - K, D11, A(K + 1, K).asArray(), 1);
          }
        } else {
          // 2-by-2 pivot block D(k)

          if (K < N - 1) {
            // Perform a rank-2 update of A(k+2:n,k+2:n) as

            // A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))**T

            // where L(k) and L(k+1) are the k-th and (k+1)-th
            // columns of L

            D21 = A[K + 1][K];
            D11 = A[K + 1][K + 1] / D21;
            D22 = A[K][K] / D21;
            T = ONE / (D11 * D22 - ONE);
            D21 = T / D21;

            for (J = K + 2; J <= N; J++) {
              WK = D21 * (D11 * A[J][K] - A[J][K + 1]);
              WKP1 = D21 * (D22 * A[J][K + 1] - A[J][K]);

              for (I = J; I <= N; I++) {
                A[I][J] -= A[I][K] * WK + A[I][K + 1] * WKP1;
              }

              A[J][K] = WK;
              A[J][K + 1] = WKP1;
            }
          }
        }
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -KP;
        IPIV[K + 1] = -KP;
      }

      // Increase K and return to the start of the main loop

      K += KSTEP;
    }
  }
}
