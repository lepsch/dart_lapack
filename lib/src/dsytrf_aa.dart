// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlasyf_aa.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dsytrf_aa(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  const ONE = 1.0;
  bool LQUERY, UPPER;
  int J, LWKMIN, LWKOPT;
  int NB, MJ, NJ, K1, K2, J1, J2, J3, JB;
  double ALPHA;

  // Determine the block size
  NB = ilaenv(1, 'DSYTRF_AA', UPLO, N, -1, -1, -1);

  // Test the input parameters.
  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1);
  if (N <= 1) {
    LWKMIN = 1;
    LWKOPT = 1;
  } else {
    LWKMIN = 2 * N;
    LWKOPT = (NB + 1) * N;
  }

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < LWKMIN && !LQUERY) {
    INFO.value = -7;
  }

  if (INFO.value == 0) {
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DSYTRF_AA', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return;
  if (N == 0) {
    return;
  }
  IPIV[1] = 1;
  if (N == 1) {
    return;
  }

  // Adjust block size based on the workspace size
  if (LWORK < (1 + NB) * N) {
    NB = (LWORK - N) ~/ N;
  }

  if (UPPER) {
    // Factorize A as U**T*D*U using the upper triangle of A

    // Copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N))
    dcopy(N, A(1, 1).asArray(), LDA, WORK(1), 1);

    // J is the main loop index, increasing from 1 to N in steps of
    // JB, where JB is the number of columns factorized by DLASYF;
    // JB is either NB, or N-J+1 for the last block
    J = 0;
    while (J < N) {
      // each step of the main loop
      //  J is the last column of the previous panel
      //  J1 is the first column of the current panel
      //  K1 identifies if the previous column of the panel has been
      //   explicitly stored, e.g., K1=1 for the first panel, and
      //   K1=0 for the rest
      J1 = J + 1;
      JB = min(N - J1 + 1, NB);
      K1 = max(1, J) - J;

      // Panel factorization
      dlasyf_aa(UPLO, 2 - K1, N - J, JB, A(max(1, J), J + 1), LDA, IPIV(J + 1),
          WORK.asMatrix(N), N, WORK(N * NB + 1));

      // Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
      for (J2 = J + 2; J2 <= min(N, J + JB + 1); J2++) {
        IPIV[J2] += J;
        if ((J2 != IPIV[J2]) && ((J1 - K1) > 2)) {
          dswap(
              J1 - K1 - 2, A(1, J2).asArray(), 1, A(1, IPIV[J2]).asArray(), 1);
        }
      }
      J += JB;

      // Trailing submatrix update, where
      //  the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and
      //  WORK stores the current block of the auxiriarly matrix H
      if (J < N) {
        // If first panel and JB=1 (NB=1), then nothing to do
        if (J1 > 1 || JB > 1) {
          // Merge rank-1 update with BLAS-3 update
          ALPHA = A[J][J + 1];
          A[J][J + 1] = ONE;
          dcopy(N - J, A(J - 1, J + 1).asArray(), LDA,
              WORK((J + 1 - J1 + 1) + JB * N), 1);
          dscal(N - J, ALPHA, WORK((J + 1 - J1 + 1) + JB * N), 1);

          // K1 identifies if the previous column of the panel has been
          //  explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
          //  while K1=0 and K2=1 for the rest
          if (J1 > 1) {
            // Not first panel
            K2 = 1;
          } else {
            // First panel
            K2 = 0;

            // First update skips the first column
            JB--;
          }

          for (J2 = J + 1; NB < 0 ? J2 >= N : J2 <= N; J2 += NB) {
            NJ = min(NB, N - J2 + 1);

            // Update (J2, J2) diagonal block with DGEMV
            J3 = J2;
            for (MJ = NJ - 1; MJ >= 1; MJ--) {
              dgemv(
                  'No transpose',
                  MJ,
                  JB + 1,
                  -ONE,
                  WORK(J3 - J1 + 1 + K1 * N).asMatrix(N),
                  N,
                  A(J1 - K2, J3).asArray(),
                  1,
                  ONE,
                  A(J3, J3).asArray(),
                  LDA);
              J3++;
            }

            // Update off-diagonal block of J2-th block row with DGEMM
            dgemm(
                'Transpose',
                'Transpose',
                NJ,
                N - J3 + 1,
                JB + 1,
                -ONE,
                A(J1 - K2, J2),
                LDA,
                WORK(J3 - J1 + 1 + K1 * N).asMatrix(N),
                N,
                ONE,
                A(J2, J3),
                LDA);
          }

          // Recover T( J, J+1 )
          A[J][J + 1] = ALPHA;
        }

        // WORK(J+1, 1) stores H(J+1, 1)
        dcopy(N - J, A(J + 1, J + 1).asArray(), LDA, WORK(1), 1);
      }
    }
  } else {
    // Factorize A as L*D*L**T using the lower triangle of A

    // copy first column A(1:N, 1) into H(1:N, 1)
    //  (stored in WORK(1:N))
    dcopy(N, A(1, 1).asArray(), 1, WORK(1), 1);

    // J is the main loop index, increasing from 1 to N in steps of
    // JB, where JB is the number of columns factorized by DLASYF;
    // JB is either NB, or N-J+1 for the last block
    J = 0;
    while (J < N) {
      // each step of the main loop
      //  J is the last column of the previous panel
      //  J1 is the first column of the current panel
      //  K1 identifies if the previous column of the panel has been
      //   explicitly stored, e.g., K1=1 for the first panel, and
      //   K1=0 for the rest
      J1 = J + 1;
      JB = min(N - J1 + 1, NB);
      K1 = max(1, J) - J;

      // Panel factorization
      dlasyf_aa(UPLO, 2 - K1, N - J, JB, A(J + 1, max(1, J)), LDA, IPIV(J + 1),
          WORK.asMatrix(N), N, WORK(N * NB + 1));

      // Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
      for (J2 = J + 2; J2 <= min(N, J + JB + 1); J2++) {
        IPIV[J2] += J;
        if ((J2 != IPIV[J2]) && ((J1 - K1) > 2)) {
          dswap(J1 - K1 - 2, A(J2, 1).asArray(), LDA, A(IPIV[J2], 1).asArray(),
              LDA);
        }
      }
      J += JB;

      // Trailing submatrix update, where
      //   A(J2+1, J1-1) stores L(J2+1, J1) and
      //   WORK(J2+1, 1) stores H(J2+1, 1)
      if (J < N) {
        // if first panel and JB=1 (NB=1), then nothing to do
        if (J1 > 1 || JB > 1) {
          // Merge rank-1 update with BLAS-3 update
          ALPHA = A[J + 1][J];
          A[J + 1][J] = ONE;
          dcopy(N - J, A(J + 1, J - 1).asArray(), 1,
              WORK((J + 1 - J1 + 1) + JB * N), 1);
          dscal(N - J, ALPHA, WORK((J + 1 - J1 + 1) + JB * N), 1);

          // K1 identifies if the previous column of the panel has been
          //  explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
          //  while K1=0 and K2=1 for the rest
          if (J1 > 1) {
            // Not first panel
            K2 = 1;
          } else {
            // First panel
            K2 = 0;

            // First update skips the first column
            JB--;
          }

          for (J2 = J + 1; NB < 0 ? J2 >= N : J2 <= N; J2 += NB) {
            NJ = min(NB, N - J2 + 1);

            // Update (J2, J2) diagonal block with DGEMV
            J3 = J2;
            for (MJ = NJ - 1; MJ >= 1; MJ--) {
              dgemv(
                  'No transpose',
                  MJ,
                  JB + 1,
                  -ONE,
                  WORK(J3 - J1 + 1 + K1 * N).asMatrix(N),
                  N,
                  A(J3, J1 - K2).asArray(),
                  LDA,
                  ONE,
                  A(J3, J3).asArray(),
                  1);
              J3++;
            }

            // Update off-diagonal block in J2-th block column with DGEMM
            dgemm(
                'No transpose',
                'Transpose',
                N - J3 + 1,
                NJ,
                JB + 1,
                -ONE,
                WORK(J3 - J1 + 1 + K1 * N).asMatrix(N),
                N,
                A(J2, J1 - K2),
                LDA,
                ONE,
                A(J3, J2),
                LDA);
          }

          // Recover T( J+1, J )
          A[J + 1][J] = ALPHA;
        }

        // WORK(J+1, 1) stores H(J+1, 1)
        dcopy(N - J, A(J + 1, J + 1).asArray(), 1, WORK(1), 1);
      }
    }
  }

  WORK[1] = LWKOPT.toDouble();
}
