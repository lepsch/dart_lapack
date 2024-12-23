// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/izamax.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlasyf_rk(
  final String UPLO,
  final int N,
  final int NB,
  final Box<int> KB,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> E_,
  final Array<int> IPIV_,
  final Matrix<Complex> W_,
  final int LDW,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();
  final W = W_.having(ld: LDW);

  const ZERO = 0.0, ONE = 1.0;
  const EIGHT = 8.0, SEVTEN = 17.0;
  bool DONE;
  int IMAX = 0,
      ITEMP,
      J,
      JB,
      JJ,
      JMAX = 0,
      K,
      KK,
      KW,
      KKW,
      KP = 0,
      KSTEP,
      P,
      II;
  double ABSAKK, ALPHA, COLMAX, ROWMAX, SFMIN, DTEMP;
  Complex D11, D12, D21, D22, R1, T;

  INFO.value = 0;

  // Initialize ALPHA for use in choosing pivot block size.
  ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

  // Compute machine safe minimum
  SFMIN = dlamch('S');

  if (lsame(UPLO, 'U')) {
    // Factorize the trailing columns of A using the upper triangle
    // of A and working backwards, and compute the matrix W = U12*D
    // for use in updating A11

    // Initialize the first entry of array E, where superdiagonal
    // elements of D are stored
    E[1] = Complex.zero;

    // K is the main loop index, decreasing from N in steps of 1 or 2
    K = N;
    while (true) {
      // KW is the column of W which corresponds to column K of A
      KW = NB + K - N;

      // Exit from loop
      if ((K <= N - NB + 1 && NB < N) || K < 1) break;

      KSTEP = 1;
      P = K;

      // Copy column K of A to column KW of W and update it
      zcopy(K, A(1, K).asArray(), 1, W(1, KW).asArray(), 1);
      if (K < N) {
        zgemv('No transpose', K, N - K, -Complex.one, A(1, K + 1), LDA,
            W(K, KW + 1).asArray(), LDW, Complex.one, W(1, KW).asArray(), 1);
      }

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used
      ABSAKK = W[K][KW].cabs1();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value.
      // Determine both COLMAX and IMAX.
      if (K > 1) {
        IMAX = izamax(K - 1, W(1, KW).asArray(), 1);
        COLMAX = W[IMAX][KW].cabs1();
      } else {
        COLMAX = ZERO;
      }

      if (max(ABSAKK, COLMAX) == ZERO) {
        // Column K is zero or underflow: set INFO and continue
        if (INFO.value == 0) INFO.value = K;
        KP = K;
        zcopy(K, W(1, KW).asArray(), 1, A(1, K).asArray(), 1);

        // Set E[K] to zero
        if (K > 1) E[K] = Complex.zero;
      } else {
        // Test for interchange

        // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
        // (used to handle NaN and Inf)
        if (!(ABSAKK < ALPHA * COLMAX)) {
          // no interchange, use 1-by-1 pivot block
          KP = K;
        } else {
          DONE = false;

          // Loop until pivot found
          do {
            // Begin pivot search loop body

            // Copy column IMAX to column KW-1 of W and update it
            zcopy(IMAX, A(1, IMAX).asArray(), 1, W(1, KW - 1).asArray(), 1);
            zcopy(K - IMAX, A(IMAX, IMAX + 1).asArray(), LDA,
                W(IMAX + 1, KW - 1).asArray(), 1);

            if (K < N) {
              zgemv(
                  'No transpose',
                  K,
                  N - K,
                  -Complex.one,
                  A(1, K + 1),
                  LDA,
                  W(IMAX, KW + 1).asArray(),
                  LDW,
                  Complex.one,
                  W(1, KW - 1).asArray(),
                  1);
            }

            // JMAX is the column-index of the largest off-diagonal
            // element in row IMAX, and ROWMAX is its absolute value.
            // Determine both ROWMAX and JMAX.
            if (IMAX != K) {
              JMAX = IMAX + izamax(K - IMAX, W(IMAX + 1, KW - 1).asArray(), 1);
              ROWMAX = W[JMAX][KW - 1].cabs1();
            } else {
              ROWMAX = ZERO;
            }

            if (IMAX > 1) {
              ITEMP = izamax(IMAX - 1, W(1, KW - 1).asArray(), 1);
              DTEMP = W[ITEMP][KW - 1].cabs1();
              if (DTEMP > ROWMAX) {
                ROWMAX = DTEMP;
                JMAX = ITEMP;
              }
            }

            // Equivalent to testing for
            // CABS1( W[IMAX][KW-1] ) >= ALPHA*ROWMAX
            // (used to handle NaN and Inf)
            if (!(W[IMAX][KW - 1].cabs1() < ALPHA * ROWMAX)) {
              // interchange rows and columns K and IMAX,
              // use 1-by-1 pivot block
              KP = IMAX;

              // copy column KW-1 of W to column KW of W
              zcopy(K, W(1, KW - 1).asArray(), 1, W(1, KW).asArray(), 1);

              DONE = true;

              // Equivalent to testing for ROWMAX == COLMAX,
              // (used to handle NaN and Inf)
            } else if ((P == JMAX) || (ROWMAX <= COLMAX)) {
              // interchange rows and columns K-1 and IMAX,
              // use 2-by-2 pivot block
              KP = IMAX;
              KSTEP = 2;
              DONE = true;
            } else {
              // Pivot not found: set params and repeat
              P = IMAX;
              COLMAX = ROWMAX;
              IMAX = JMAX;

              // Copy updated JMAXth (next IMAXth) column to Kth of W
              zcopy(K, W(1, KW - 1).asArray(), 1, W(1, KW).asArray(), 1);
            }

            // End pivot search loop body
          } while (!DONE);
        }

        KK = K - KSTEP + 1;

        // KKW is the column of W which corresponds to column KK of A
        KKW = NB + KK - N;

        if ((KSTEP == 2) && (P != K)) {
          // Copy non-updated column K to column P
          zcopy(K - P, A(P + 1, K).asArray(), 1, A(P, P + 1).asArray(), LDA);
          zcopy(P, A(1, K).asArray(), 1, A(1, P).asArray(), 1);

          // Interchange rows K and P in last N-K+1 columns of A
          // and last N-K+2 columns of W
          zswap(N - K + 1, A(K, K).asArray(), LDA, A(P, K).asArray(), LDA);
          zswap(N - KK + 1, W(K, KKW).asArray(), LDW, W(P, KKW).asArray(), LDW);
        }

        // Updated column KP is already stored in column KKW of W
        if (KP != KK) {
          // Copy non-updated column KK to column KP
          A[KP][K] = A[KK][K];
          zcopy(K - 1 - KP, A(KP + 1, KK).asArray(), 1, A(KP, KP + 1).asArray(),
              LDA);
          zcopy(KP, A(1, KK).asArray(), 1, A(1, KP).asArray(), 1);

          // Interchange rows KK and KP in last N-KK+1 columns
          // of A and W
          zswap(N - KK + 1, A(KK, KK).asArray(), LDA, A(KP, KK).asArray(), LDA);
          zswap(
              N - KK + 1, W(KK, KKW).asArray(), LDW, W(KP, KKW).asArray(), LDW);
        }

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column KW of W now holds
          //
          // W(k) = U(k)*D(k)
          //
          // where U(k) is the k-th column of U

          // Store U(k) in column k of A
          zcopy(K, W(1, KW).asArray(), 1, A(1, K).asArray(), 1);
          if (K > 1) {
            if (A[K][K].cabs1() >= SFMIN) {
              R1 = Complex.one / A[K][K];
              zscal(K - 1, R1, A(1, K).asArray(), 1);
            } else if (A[K][K] != Complex.zero) {
              for (II = 1; II <= K - 1; II++) {
                A[II][K] /= A[K][K];
              }
            }

            // Store the superdiagonal element of D in array E
            E[K] = Complex.zero;
          }
        } else {
          // 2-by-2 pivot block D(k): columns KW and KW-1 of W now
          // hold
          //
          // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
          //
          // where U(k) and U(k-1) are the k-th and (k-1)-th columns
          // of U
          if (K > 2) {
            // Store U(k) and U(k-1) in columns k and k-1 of A
            D12 = W[K - 1][KW];
            D11 = W[K][KW] / D12;
            D22 = W[K - 1][KW - 1] / D12;
            T = Complex.one / (D11 * D22 - Complex.one);
            for (J = 1; J <= K - 2; J++) {
              A[J][K - 1] = T * ((D11 * W[J][KW - 1] - W[J][KW]) / D12);
              A[J][K] = T * ((D22 * W[J][KW] - W[J][KW - 1]) / D12);
            }
          }

          // Copy diagonal elements of D(K) to A,
          // copy superdiagonal element of D(K) to E(K) and
          // ZERO out superdiagonal entry of A
          A[K - 1][K - 1] = W[K - 1][KW - 1];
          A[K - 1][K] = Complex.zero;
          A[K][K] = W[K][KW];
          E[K] = W[K - 1][KW];
          E[K - 1] = Complex.zero;
        }

        // End column K is nonsingular
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

    // Update the upper triangle of A11 (= A(1:k,1:k)) as
    //
    // A11 := A11 - U12*D*U12**T = A11 - U12*W**T

    // computing blocks of NB columns at a time
    for (J = ((K - 1) ~/ NB) * NB + 1; -NB < 0 ? J >= 1 : J <= 1; J += -NB) {
      JB = min(NB, K - J + 1);

      // Update the upper triangle of the diagonal block
      for (JJ = J; JJ <= J + JB - 1; JJ++) {
        zgemv('No transpose', JJ - J + 1, N - K, -Complex.one, A(J, K + 1), LDA,
            W(JJ, KW + 1).asArray(), LDW, Complex.one, A(J, JJ).asArray(), 1);
      }

      // Update the rectangular superdiagonal block
      if (J >= 2) {
        zgemm('No transpose', 'Transpose', J - 1, JB, N - K, -Complex.one,
            A(1, K + 1), LDA, W(J, KW + 1), LDW, Complex.one, A(1, J), LDA);
      }
    }

    // Set KB to the number of columns factorized
    KB.value = N - K;
  } else {
    // Factorize the leading columns of A using the lower triangle
    // of A and working forwards, and compute the matrix W = L21*D
    // for use in updating A22

    // Initialize the unused last entry of the subdiagonal array E.
    E[N] = Complex.zero;

    // K is the main loop index, increasing from 1 in steps of 1 or 2
    K = 1;
    while (!((K >= NB && NB < N) || K > N)) {
      KSTEP = 1;
      P = K;

      // Copy column K of A to column K of W and update it

      zcopy(N - K + 1, A(K, K).asArray(), 1, W(K, K).asArray(), 1);
      if (K > 1) {
        zgemv('No transpose', N - K + 1, K - 1, -Complex.one, A(K, 1), LDA,
            W(K, 1).asArray(), LDW, Complex.one, W(K, K).asArray(), 1);
      }

      // Determine rows and columns to be interchanged and whether
      // a 1-by-1 or 2-by-2 pivot block will be used
      ABSAKK = W[K][K].cabs1();

      // IMAX is the row-index of the largest off-diagonal element in
      // column K, and COLMAX is its absolute value.
      // Determine both COLMAX and IMAX.
      if (K < N) {
        IMAX = K + izamax(N - K, W(K + 1, K).asArray(), 1);
        COLMAX = W[IMAX][K].cabs1();
      } else {
        COLMAX = ZERO;
      }

      if (max(ABSAKK, COLMAX) == ZERO) {
        // Column K is zero or underflow: set INFO and continue
        if (INFO.value == 0) INFO.value = K;
        KP = K;
        zcopy(N - K + 1, W(K, K).asArray(), 1, A(K, K).asArray(), 1);

        // Set E[K] to zero
        if (K < N) E[K] = Complex.zero;
      } else {
        // Test for interchange

        // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
        // (used to handle NaN and Inf)
        if (!(ABSAKK < ALPHA * COLMAX)) {
          // no interchange, use 1-by-1 pivot block
          KP = K;
        } else {
          DONE = false;

          // Loop until pivot found
          do {
            // Begin pivot search loop body

            // Copy column IMAX to column K+1 of W and update it
            zcopy(
                IMAX - K, A(IMAX, K).asArray(), LDA, W(K, K + 1).asArray(), 1);
            zcopy(N - IMAX + 1, A(IMAX, IMAX).asArray(), 1,
                W(IMAX, K + 1).asArray(), 1);
            if (K > 1) {
              zgemv(
                  'No transpose',
                  N - K + 1,
                  K - 1,
                  -Complex.one,
                  A(K, 1),
                  LDA,
                  W(IMAX, 1).asArray(),
                  LDW,
                  Complex.one,
                  W(K, K + 1).asArray(),
                  1);
            }

            // JMAX is the column-index of the largest off-diagonal
            // element in row IMAX, and ROWMAX is its absolute value.
            // Determine both ROWMAX and JMAX.
            if (IMAX != K) {
              JMAX = K - 1 + izamax(IMAX - K, W(K, K + 1).asArray(), 1);
              ROWMAX = W[JMAX][K + 1].cabs1();
            } else {
              ROWMAX = ZERO;
            }

            if (IMAX < N) {
              ITEMP = IMAX + izamax(N - IMAX, W(IMAX + 1, K + 1).asArray(), 1);
              DTEMP = W[ITEMP][K + 1].cabs1();
              if (DTEMP > ROWMAX) {
                ROWMAX = DTEMP;
                JMAX = ITEMP;
              }
            }

            // Equivalent to testing for
            // CABS1( W( IMAX, K+1 ) ) >= ALPHA*ROWMAX
            // (used to handle NaN and Inf)
            if (!(W[IMAX][K + 1].cabs1() < ALPHA * ROWMAX)) {
              // interchange rows and columns K and IMAX,
              // use 1-by-1 pivot block
              KP = IMAX;

              // copy column K+1 of W to column K of W
              zcopy(N - K + 1, W(K, K + 1).asArray(), 1, W(K, K).asArray(), 1);

              DONE = true;

              // Equivalent to testing for ROWMAX == COLMAX,
              // (used to handle NaN and Inf)
            } else if ((P == JMAX) || (ROWMAX <= COLMAX)) {
              // interchange rows and columns K+1 and IMAX,
              // use 2-by-2 pivot block
              KP = IMAX;
              KSTEP = 2;
              DONE = true;
            } else {
              // Pivot not found: set params and repeat
              P = IMAX;
              COLMAX = ROWMAX;
              IMAX = JMAX;

              // Copy updated JMAXth (next IMAXth) column to Kth of W
              zcopy(N - K + 1, W(K, K + 1).asArray(), 1, W(K, K).asArray(), 1);
            }

            // End pivot search loop body
          } while (!DONE);
        }

        KK = K + KSTEP - 1;

        if ((KSTEP == 2) && (P != K)) {
          // Copy non-updated column K to column P
          zcopy(P - K, A(K, K).asArray(), 1, A(P, K).asArray(), LDA);
          zcopy(N - P + 1, A(P, K).asArray(), 1, A(P, P).asArray(), 1);

          // Interchange rows K and P in first K columns of A
          // and first K+1 columns of W
          zswap(K, A(K, 1).asArray(), LDA, A(P, 1).asArray(), LDA);
          zswap(KK, W(K, 1).asArray(), LDW, W(P, 1).asArray(), LDW);
        }

        // Updated column KP is already stored in column KK of W
        if (KP != KK) {
          // Copy non-updated column KK to column KP
          A[KP][K] = A[KK][K];
          zcopy(KP - K - 1, A(K + 1, KK).asArray(), 1, A(KP, K + 1).asArray(),
              LDA);
          zcopy(N - KP + 1, A(KP, KK).asArray(), 1, A(KP, KP).asArray(), 1);

          // Interchange rows KK and KP in first KK columns of A and W
          zswap(KK, A(KK, 1).asArray(), LDA, A(KP, 1).asArray(), LDA);
          zswap(KK, W(KK, 1).asArray(), LDW, W(KP, 1).asArray(), LDW);
        }

        if (KSTEP == 1) {
          // 1-by-1 pivot block D(k): column k of W now holds
          //
          // W(k) = L(k)*D(k)
          //
          // where L(k) is the k-th column of L

          // Store L(k) in column k of A
          zcopy(N - K + 1, W(K, K).asArray(), 1, A(K, K).asArray(), 1);
          if (K < N) {
            if (A[K][K].cabs1() >= SFMIN) {
              R1 = Complex.one / A[K][K];
              zscal(N - K, R1, A(K + 1, K).asArray(), 1);
            } else if (A[K][K] != Complex.zero) {
              for (II = K + 1; II <= N; II++) {
                A[II][K] /= A[K][K];
              }
            }

            // Store the subdiagonal element of D in array E
            E[K] = Complex.zero;
          }
        } else {
          // 2-by-2 pivot block D(k): columns k and k+1 of W now hold
          //
          // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
          //
          // where L(k) and L(k+1) are the k-th and (k+1)-th columns
          // of L
          if (K < N - 1) {
            // Store L(k) and L(k+1) in columns k and k+1 of A
            D21 = W[K + 1][K];
            D11 = W[K + 1][K + 1] / D21;
            D22 = W[K][K] / D21;
            T = Complex.one / (D11 * D22 - Complex.one);
            for (J = K + 2; J <= N; J++) {
              A[J][K] = T * ((D11 * W[J][K] - W[J][K + 1]) / D21);
              A[J][K + 1] = T * ((D22 * W[J][K + 1] - W[J][K]) / D21);
            }
          }

          // Copy diagonal elements of D(K) to A,
          // copy subdiagonal element of D(K) to E(K) and
          // ZERO out subdiagonal entry of A
          A[K][K] = W[K][K];
          A[K + 1][K] = Complex.zero;
          A[K + 1][K + 1] = W[K + 1][K + 1];
          E[K] = W[K + 1][K];
          E[K + 1] = Complex.zero;
        }

        // End column K is nonsingular
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

    // Update the lower triangle of A22 (= A(k:n,k:n)) as
    //
    // A22 := A22 - L21*D*L21**T = A22 - L21*W**T
    //
    // computing blocks of NB columns at a time

    for (J = K; NB < 0 ? J >= N : J <= N; J += NB) {
      JB = min(NB, N - J + 1);

      // Update the lower triangle of the diagonal block
      for (JJ = J; JJ <= J + JB - 1; JJ++) {
        zgemv('No transpose', J + JB - JJ, K - 1, -Complex.one, A(JJ, 1), LDA,
            W(JJ, 1).asArray(), LDW, Complex.one, A(JJ, JJ).asArray(), 1);
      }

      // Update the rectangular subdiagonal block
      if (J + JB <= N) {
        zgemm(
            'No transpose',
            'Transpose',
            N - J - JB + 1,
            JB,
            K - 1,
            -Complex.one,
            A(J + JB, 1),
            LDA,
            W(J, 1),
            LDW,
            Complex.one,
            A(J + JB, J),
            LDA);
      }
    }

    // Set KB to the number of columns factorized
    KB.value = K - 1;
  }
}
