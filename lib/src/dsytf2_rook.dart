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
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsytf2_rook(
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
  bool UPPER, DONE;
  int I, IMAX = 0, J, JMAX = 0, ITEMP, K, KK, KP = 0, KSTEP, P, II;
  double ABSAKK,
      ALPHA,
      COLMAX,
      D11,
      D12,
      D21,
      D22,
      ROWMAX,
      DTEMP,
      T,
      WK,
      WKM1,
      WKP1,
      SFMIN;

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
    xerbla('DSYTF2_ROOK', -INFO.value);
    return;
  }

  // Initialize ALPHA for use in choosing pivot block size.

  ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  // Compute machine safe minimum

  SFMIN = dlamch('S');

  if (UPPER) {
    // Factorize A as U*D*U**T using the upper triangle of A

    // K is the main loop index, decreasing from N to 1 in steps of
    // 1 or 2

    K = N;
    while (K >= 1) {
      KSTEP = 1;
      P = K;

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

      if ((max(ABSAKK, COLMAX) == ZERO)) {
        // Column K is zero or underflow: set INFO and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
      } else {
        // Test for interchange

        // Equivalent to testing for (used to handle NaN and Inf)
        // ABSAKK >= ALPHA*COLMAX

        if (!(ABSAKK < ALPHA * COLMAX)) {
          // no interchange,
          // use 1-by-1 pivot block

          KP = K;
        } else {
          DONE = false;

          // Loop until pivot found

          do {
            // Begin pivot search loop body

            // JMAX is the column-index of the largest off-diagonal
            // element in row IMAX, and ROWMAX is its absolute value.
            // Determine both ROWMAX and JMAX.

            if (IMAX != K) {
              JMAX = IMAX + idamax(K - IMAX, A(IMAX, IMAX + 1).asArray(), LDA);
              ROWMAX = A[IMAX][JMAX].abs();
            } else {
              ROWMAX = ZERO;
            }

            if (IMAX > 1) {
              ITEMP = idamax(IMAX - 1, A(1, IMAX).asArray(), 1);
              DTEMP = A[ITEMP][IMAX].abs();
              if (DTEMP > ROWMAX) {
                ROWMAX = DTEMP;
                JMAX = ITEMP;
              }
            }

            // Equivalent to testing for (used to handle NaN and Inf)
            // ABS( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX

            if (!(A[IMAX][IMAX].abs() < ALPHA * ROWMAX)) {
              // interchange rows and columns K and IMAX,
              // use 1-by-1 pivot block

              KP = IMAX;
              DONE = true;

              // Equivalent to testing for ROWMAX == COLMAX,
              // used to handle NaN and Inf
            } else if ((P == JMAX) || (ROWMAX <= COLMAX)) {
              // interchange rows and columns K+1 and IMAX,
              // use 2-by-2 pivot block

              KP = IMAX;
              KSTEP = 2;
              DONE = true;
            } else {
              // Pivot NOT found, set variables and repeat

              P = IMAX;
              COLMAX = ROWMAX;
              IMAX = JMAX;
            }

            // End pivot search loop body
          } while (!DONE);
        }

        // Swap TWO rows and TWO columns

        // First swap

        if ((KSTEP == 2) && (P != K)) {
          // Interchange rows and column K and P in the leading
          // submatrix A(1:k,1:k) if we have a 2-by-2 pivot

          if (P > 1) dswap(P - 1, A(1, K).asArray(), 1, A(1, P).asArray(), 1);
          if (P < (K - 1)) {
            dswap(K - P - 1, A(P + 1, K).asArray(), 1, A(P, P + 1).asArray(),
                LDA);
          }
          T = A[K][K];
          A[K][K] = A[P][P];
          A[P][P] = T;
        }

        // Second swap

        KK = K - KSTEP + 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the leading
          // submatrix A(1:k,1:k)

          if (KP > 1) {
            dswap(KP - 1, A(1, KK).asArray(), 1, A(1, KP).asArray(), 1);
          }
          if ((KK > 1) && (KP < (KK - 1))) {
            dswap(KK - KP - 1, A(KP + 1, KK).asArray(), 1,
                A(KP, KP + 1).asArray(), LDA);
          }
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

          if (K > 1) {
            // Perform a rank-1 update of A(1:k-1,1:k-1) and
            // store U(k) in column k

            if (A[K][K].abs() >= SFMIN) {
              // Perform a rank-1 update of A(1:k-1,1:k-1) as
              // A := A - U(k)*D(k)*U(k)**T
              //    = A - W(k)*1/D(k)*W(k)**T

              D11 = ONE / A[K][K];
              dsyr(UPLO, K - 1, -D11, A(1, K).asArray(), 1, A, LDA);

              // Store U(k) in column k

              dscal(K - 1, D11, A(1, K).asArray(), 1);
            } else {
              // Store L(k) in column K

              D11 = A[K][K];
              for (II = 1; II <= K - 1; II++) {
                A[II][K] /= D11;
              }

              // Perform a rank-1 update of A(k+1:n,k+1:n) as
              // A := A - U(k)*D(k)*U(k)**T
              //    = A - W(k)*(1/D(k))*W(k)**T
              //    = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

              dsyr(UPLO, K - 1, -D11, A(1, K).asArray(), 1, A, LDA);
            }
          }
        } else {
          // 2-by-2 pivot block D(k): columns k and k-1 now hold

          // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

          // where U(k) and U(k-1) are the k-th and (k-1)-th columns
          // of U

          // Perform a rank-2 update of A(1:k-2,1:k-2) as

          // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
          //    = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T

          // and store L(k) and L(k+1) in columns k and k+1

          if (K > 2) {
            D12 = A[K - 1][K];
            D22 = A[K - 1][K - 1] / D12;
            D11 = A[K][K] / D12;
            T = ONE / (D11 * D22 - ONE);

            for (J = K - 2; J >= 1; J--) {
              WKM1 = T * (D11 * A[J][K - 1] - A[J][K]);
              WK = T * (D22 * A[J][K] - A[J][K - 1]);

              for (I = J; I >= 1; I--) {
                A[I][J] =
                    A[I][J] - (A[I][K] / D12) * WK - (A[I][K - 1] / D12) * WKM1;
              }

              // Store U(k) and U(k-1) in cols k and k-1 for row J

              A[J][K] = WK / D12;
              A[J][K - 1] = WKM1 / D12;
            }
          }
        }
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -P;
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
      P = K;

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

      if ((max(ABSAKK, COLMAX) == ZERO)) {
        // Column K is zero or underflow: set INFO and continue

        if (INFO.value == 0) INFO.value = K;
        KP = K;
      } else {
        // Test for interchange

        // Equivalent to testing for (used to handle NaN and Inf)
        // ABSAKK >= ALPHA*COLMAX

        if (!(ABSAKK < ALPHA * COLMAX)) {
          // no interchange, use 1-by-1 pivot block

          KP = K;
        } else {
          DONE = false;

          // Loop until pivot found

          do {
            // Begin pivot search loop body

            // JMAX is the column-index of the largest off-diagonal
            // element in row IMAX, and ROWMAX is its absolute value.
            // Determine both ROWMAX and JMAX.

            if (IMAX != K) {
              JMAX = K - 1 + idamax(IMAX - K, A(IMAX, K).asArray(), LDA);
              ROWMAX = A[IMAX][JMAX].abs();
            } else {
              ROWMAX = ZERO;
            }

            if (IMAX < N) {
              ITEMP = IMAX + idamax(N - IMAX, A(IMAX + 1, IMAX).asArray(), 1);
              DTEMP = A[ITEMP][IMAX].abs();
              if (DTEMP > ROWMAX) {
                ROWMAX = DTEMP;
                JMAX = ITEMP;
              }
            }

            // Equivalent to testing for (used to handle NaN and Inf)
            // ABS( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX

            if (!(A[IMAX][IMAX].abs() < ALPHA * ROWMAX)) {
              // interchange rows and columns K and IMAX,
              // use 1-by-1 pivot block

              KP = IMAX;
              DONE = true;

              // Equivalent to testing for ROWMAX == COLMAX,
              // used to handle NaN and Inf
            } else if ((P == JMAX) || (ROWMAX <= COLMAX)) {
              // interchange rows and columns K+1 and IMAX,
              // use 2-by-2 pivot block

              KP = IMAX;
              KSTEP = 2;
              DONE = true;
            } else {
              // Pivot NOT found, set variables and repeat

              P = IMAX;
              COLMAX = ROWMAX;
              IMAX = JMAX;
            }

            // End pivot search loop body
          } while (!DONE);
        }

        // Swap TWO rows and TWO columns

        // First swap

        if ((KSTEP == 2) && (P != K)) {
          // Interchange rows and column K and P in the trailing
          // submatrix A(k:n,k:n) if we have a 2-by-2 pivot

          if (P < N) {
            dswap(N - P, A(P + 1, K).asArray(), 1, A(P + 1, P).asArray(), 1);
          }
          if (P > (K + 1)) {
            dswap(P - K - 1, A(K + 1, K).asArray(), 1, A(P, K + 1).asArray(),
                LDA);
          }
          T = A[K][K];
          A[K][K] = A[P][P];
          A[P][P] = T;
        }

        // Second swap

        KK = K + KSTEP - 1;
        if (KP != KK) {
          // Interchange rows and columns KK and KP in the trailing
          // submatrix A(k:n,k:n)

          if (KP < N) {
            dswap(
                N - KP, A(KP + 1, KK).asArray(), 1, A(KP + 1, KP).asArray(), 1);
          }
          if ((KK < N) && (KP > (KK + 1))) {
            dswap(KP - KK - 1, A(KK + 1, KK).asArray(), 1,
                A(KP, KK + 1).asArray(), LDA);
          }
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
            // Perform a rank-1 update of A(k+1:n,k+1:n) and
            // store L(k) in column k

            if (A[K][K].abs() >= SFMIN) {
              // Perform a rank-1 update of A(k+1:n,k+1:n) as
              // A := A - L(k)*D(k)*L(k)**T
              //    = A - W(k)*(1/D(k))*W(k)**T

              D11 = ONE / A[K][K];
              dsyr(UPLO, N - K, -D11, A(K + 1, K).asArray(), 1, A(K + 1, K + 1),
                  LDA);

              // Store L(k) in column k

              dscal(N - K, D11, A(K + 1, K).asArray(), 1);
            } else {
              // Store L(k) in column k

              D11 = A[K][K];
              for (II = K + 1; II <= N; II++) {
                A[II][K] /= D11;
              }

              // Perform a rank-1 update of A(k+1:n,k+1:n) as
              // A := A - L(k)*D(k)*L(k)**T
              //    = A - W(k)*(1/D(k))*W(k)**T
              //    = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

              dsyr(UPLO, N - K, -D11, A(K + 1, K).asArray(), 1, A(K + 1, K + 1),
                  LDA);
            }
          }
        } else {
          // 2-by-2 pivot block D(k): columns k and k+1 now hold

          // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

          // where L(k) and L(k+1) are the k-th and (k+1)-th columns
          // of L

          // Perform a rank-2 update of A(k+2:n,k+2:n) as

          // A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
          //    = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T

          // and store L(k) and L(k+1) in columns k and k+1

          if (K < N - 1) {
            D21 = A[K + 1][K];
            D11 = A[K + 1][K + 1] / D21;
            D22 = A[K][K] / D21;
            T = ONE / (D11 * D22 - ONE);

            for (J = K + 2; J <= N; J++) {
              // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J

              WK = T * (D11 * A[J][K] - A[J][K + 1]);
              WKP1 = T * (D22 * A[J][K + 1] - A[J][K]);

              // Perform a rank-2 update of A(k+2:n,k+2:n)

              for (I = J; I <= N; I++) {
                A[I][J] =
                    A[I][J] - (A[I][K] / D21) * WK - (A[I][K + 1] / D21) * WKP1;
              }

              // Store L(k) and L(k+1) in cols k and k+1 for row J

              A[J][K] = WK / D21;
              A[J][K + 1] = WKP1 / D21;
            }
          }
        }
      }

      // Store details of the interchanges in IPIV

      if (KSTEP == 1) {
        IPIV[K] = KP;
      } else {
        IPIV[K] = -P;
        IPIV[K + 1] = -KP;
      }

      // Increase K and return to the start of the main loop

      K += KSTEP;
    }
  }
}
