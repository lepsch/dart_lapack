// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/dtrsm.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbtrf.dart';
import 'package:lapack/src/dgetrf.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dsygst.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytrf_aa_2stage(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TB_,
  final int LTB,
  final Array<int> IPIV_,
  final Array<int> IPIV2_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TB = TB_.having();
  final IPIV = IPIV_.having();
  final IPIV2 = IPIV2_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool UPPER, TQUERY, WQUERY;
  int I, J, K, I1, I2, TD;
  int LDTB, NB, KB, JB, NT;
  double PIV;
  final IINFO = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  WQUERY = (LWORK == -1);
  TQUERY = (LTB == -1);
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LTB < max(1, 4 * N) && !TQUERY) {
    INFO.value = -6;
  } else if (LWORK < max(1, N) && !WQUERY) {
    INFO.value = -10;
  }

  if (INFO.value != 0) {
    xerbla('DSYTRF_AA_2STAGE', -INFO.value);
    return;
  }

  // Answer the query

  NB = ilaenv(1, 'DSYTRF_AA_2STAGE', UPLO, N, -1, -1, -1);
  if (INFO.value == 0) {
    if (TQUERY) {
      TB[1] = max(1, (3 * NB + 1) * N).toDouble();
    }
    if (WQUERY) {
      WORK[1] = max(1, N * NB).toDouble();
    }
  }
  if (TQUERY || WQUERY) {
    return;
  }

  // Quick return;

  if (N == 0) {
    return;
  }

  // Determine the number of the block size

  LDTB = LTB ~/ N;
  if (LDTB < 3 * NB + 1) {
    NB = (LDTB - 1) ~/ 3;
  }
  if (LWORK < NB * N) {
    NB = LWORK ~/ N;
  }

  // Determine the number of the block columns

  NT = (N + NB - 1) ~/ NB;
  TD = 2 * NB;
  KB = min(NB, N);

  // Initialize vectors/matrices

  for (J = 1; J <= KB; J++) {
    IPIV[J] = J;
  }

  // Save NB

  TB[1] = NB.toDouble();

  if (UPPER) {
    // .....................................................
    // Factorize A as U**T*D*U using the upper triangle of A
    // .....................................................

    for (J = 0; J <= NT - 1; J++) {
      // Generate Jth column of W and H

      KB = min(NB, N - J * NB);
      for (I = 1; I <= J - 1; I++) {
        if (I == 1) {
          // H(I,J) = T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
          if (I == (J - 1)) {
            JB = NB + KB;
          } else {
            JB = 2 * NB;
          }
          dgemm(
              'NoTranspose',
              'NoTranspose',
              NB,
              KB,
              JB,
              ONE,
              TB(TD + 1 + (I * NB) * LDTB).asMatrix(LDTB - 1),
              LDTB - 1,
              A((I - 1) * NB + 1, J * NB + 1),
              LDA,
              ZERO,
              WORK(I * NB + 1).asMatrix(N),
              N);
        } else {
          // H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
          if (I == J - 1) {
            JB = 2 * NB + KB;
          } else {
            JB = 3 * NB;
          }
          dgemm(
              'NoTranspose',
              'NoTranspose',
              NB,
              KB,
              JB,
              ONE,
              TB(TD + NB + 1 + ((I - 1) * NB) * LDTB).asMatrix(LDTB - 1),
              LDTB - 1,
              A((I - 2) * NB + 1, J * NB + 1),
              LDA,
              ZERO,
              WORK(I * NB + 1).asMatrix(N),
              N);
        }
      }

      // Compute T(J,J)

      dlacpy('Upper', KB, KB, A(J * NB + 1, J * NB + 1), LDA,
          TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1), LDTB - 1);
      if (J > 1) {
        // T(J,J) = U(1:J,J)'*H(1:J)
        dgemm(
            'Transpose',
            'NoTranspose',
            KB,
            KB,
            (J - 1) * NB,
            -ONE,
            A(1, J * NB + 1),
            LDA,
            WORK(NB + 1).asMatrix(N),
            N,
            ONE,
            TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1);
        // T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
        dgemm(
            'Transpose',
            'NoTranspose',
            KB,
            NB,
            KB,
            ONE,
            A((J - 1) * NB + 1, J * NB + 1),
            LDA,
            TB(TD + NB + 1 + ((J - 1) * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1,
            ZERO,
            WORK(1).asMatrix(N),
            N);
        dgemm(
            'NoTranspose',
            'NoTranspose',
            KB,
            KB,
            NB,
            -ONE,
            WORK(1).asMatrix(N),
            N,
            A((J - 2) * NB + 1, J * NB + 1),
            LDA,
            ONE,
            TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1);
      }
      if (J > 0) {
        dsygst(1, 'Upper', KB, TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1, A((J - 1) * NB + 1, J * NB + 1), LDA, IINFO);
      }

      // Expand T(J,J) into full format

      for (I = 1; I <= KB; I++) {
        for (K = I + 1; K <= KB; K++) {
          TB[TD + (K - I) + 1 + (J * NB + I - 1) * LDTB] =
              TB[TD - (K - (I + 1)) + (J * NB + K - 1) * LDTB];
        }
      }

      if (J < NT - 1) {
        if (J > 0) {
          // Compute H(J,J)

          if (J == 1) {
            dgemm(
                'NoTranspose',
                'NoTranspose',
                KB,
                KB,
                KB,
                ONE,
                TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
                LDTB - 1,
                A((J - 1) * NB + 1, J * NB + 1),
                LDA,
                ZERO,
                WORK(J * NB + 1).asMatrix(N),
                N);
          } else {
            dgemm(
                'NoTranspose',
                'NoTranspose',
                KB,
                KB,
                NB + KB,
                ONE,
                TB(TD + NB + 1 + ((J - 1) * NB) * LDTB).asMatrix(LDTB - 1),
                LDTB - 1,
                A((J - 2) * NB + 1, J * NB + 1),
                LDA,
                ZERO,
                WORK(J * NB + 1).asMatrix(N),
                N);
          }

          // Update with the previous column

          dgemm(
              'Transpose',
              'NoTranspose',
              NB,
              N - (J + 1) * NB,
              J * NB,
              -ONE,
              WORK(NB + 1).asMatrix(N),
              N,
              A(1, (J + 1) * NB + 1),
              LDA,
              ONE,
              A(J * NB + 1, (J + 1) * NB + 1),
              LDA);
        }

        // Copy panel to workspace to call DGETRF

        for (K = 1; K <= NB; K++) {
          dcopy(N - (J + 1) * NB, A(J * NB + K, (J + 1) * NB + 1).asArray(),
              LDA, WORK(1 + (K - 1) * N), 1);
        }

        // Factorize panel

        dgetrf(N - (J + 1) * NB, NB, WORK.asMatrix(N), N,
            IPIV((J + 1) * NB + 1), IINFO);
        // IF (IINFO != 0 && INFO == 0) THEN
        //    INFO = IINFO+(J+1)*NB
        // END IF

        // Copy panel back

        for (K = 1; K <= NB; K++) {
          dcopy(N - (J + 1) * NB, WORK(1 + (K - 1) * N), 1,
              A(J * NB + K, (J + 1) * NB + 1).asArray(), LDA);
        }

        // Compute T(J+1, J), zero out for GEMM update

        KB = min(NB, N - (J + 1) * NB);
        dlaset('Full', KB, NB, ZERO, ZERO,
            TB(TD + NB + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1), LDTB - 1);
        dlacpy('Upper', KB, NB, WORK.asMatrix(N), N,
            TB(TD + NB + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1), LDTB - 1);
        if (J > 0) {
          dtrsm(
              'R',
              'U',
              'N',
              'U',
              KB,
              NB,
              ONE,
              A((J - 1) * NB + 1, J * NB + 1),
              LDA,
              TB(TD + NB + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
              LDTB - 1);
        }

        // Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
        // updates

        for (K = 1; K <= NB; K++) {
          for (I = 1; I <= KB; I++) {
            TB[TD - NB + K - I + 1 + (J * NB + NB + I - 1) * LDTB] =
                TB[TD + NB + I - K + 1 + (J * NB + K - 1) * LDTB];
          }
        }
        dlaset(
            'Lower', KB, NB, ZERO, ONE, A(J * NB + 1, (J + 1) * NB + 1), LDA);

        // Apply pivots to trailing submatrix of A

        for (K = 1; K <= KB; K++) {
          // > Adjust ipiv
          IPIV[(J + 1) * NB + K] += (J + 1) * NB;

          I1 = (J + 1) * NB + K;
          I2 = IPIV[(J + 1) * NB + K];
          if (I1 != I2) {
            // > Apply pivots to previous columns of L
            dswap(K - 1, A((J + 1) * NB + 1, I1).asArray(), 1,
                A((J + 1) * NB + 1, I2).asArray(), 1);
            // > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
            if (I2 > (I1 + 1)) {
              dswap(I2 - I1 - 1, A(I1, I1 + 1).asArray(), LDA,
                  A(I1 + 1, I2).asArray(), 1);
            }
            // > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
            if (I2 < N) {
              dswap(N - I2, A(I1, I2 + 1).asArray(), LDA,
                  A(I2, I2 + 1).asArray(), LDA);
            }
            // > Swap A(I1, I1) with A(I2, I2)
            PIV = A[I1][I1];
            A[I1][I1] = A[I2][I2];
            A[I2][I2] = PIV;
            // > Apply pivots to previous columns of L
            if (J > 0) {
              dswap(J * NB, A(1, I1).asArray(), 1, A(1, I2).asArray(), 1);
            }
          }
        }
      }
    }
  } else {
    // .....................................................
    // Factorize A as L*D*L**T using the lower triangle of A
    // .....................................................

    for (J = 0; J <= NT - 1; J++) {
      // Generate Jth column of W and H

      KB = min(NB, N - J * NB);
      for (I = 1; I <= J - 1; I++) {
        if (I == 1) {
          // H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)'
          if (I == J - 1) {
            JB = NB + KB;
          } else {
            JB = 2 * NB;
          }
          dgemm(
              'NoTranspose',
              'Transpose',
              NB,
              KB,
              JB,
              ONE,
              TB(TD + 1 + (I * NB) * LDTB).asMatrix(LDTB - 1),
              LDTB - 1,
              A(J * NB + 1, (I - 1) * NB + 1),
              LDA,
              ZERO,
              WORK(I * NB + 1).asMatrix(N),
              N);
        } else {
          // H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
          if (I == J - 1) {
            JB = 2 * NB + KB;
          } else {
            JB = 3 * NB;
          }
          dgemm(
              'NoTranspose',
              'Transpose',
              NB,
              KB,
              JB,
              ONE,
              TB(TD + NB + 1 + ((I - 1) * NB) * LDTB).asMatrix(LDTB - 1),
              LDTB - 1,
              A(J * NB + 1, (I - 2) * NB + 1),
              LDA,
              ZERO,
              WORK(I * NB + 1).asMatrix(N),
              N);
        }
      }

      // Compute T(J,J)

      dlacpy('Lower', KB, KB, A(J * NB + 1, J * NB + 1), LDA,
          TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1), LDTB - 1);
      if (J > 1) {
        // T(J,J) = L(J,1:J)*H(1:J)
        dgemm(
            'NoTranspose',
            'NoTranspose',
            KB,
            KB,
            (J - 1) * NB,
            -ONE,
            A(J * NB + 1, 1),
            LDA,
            WORK(NB + 1).asMatrix(N),
            N,
            ONE,
            TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1);
        // T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
        dgemm(
            'NoTranspose',
            'NoTranspose',
            KB,
            NB,
            KB,
            ONE,
            A(J * NB + 1, (J - 1) * NB + 1),
            LDA,
            TB(TD + NB + 1 + ((J - 1) * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1,
            ZERO,
            WORK(1).asMatrix(N),
            N);
        dgemm(
            'NoTranspose',
            'Transpose',
            KB,
            KB,
            NB,
            -ONE,
            WORK(1).asMatrix(N),
            N,
            A(J * NB + 1, (J - 2) * NB + 1),
            LDA,
            ONE,
            TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1);
      }
      if (J > 0) {
        dsygst(1, 'Lower', KB, TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
            LDTB - 1, A(J * NB + 1, (J - 1) * NB + 1), LDA, IINFO);
      }

      // Expand T(J,J) into full format

      for (I = 1; I <= KB; I++) {
        for (K = I + 1; K <= KB; K++) {
          TB[TD - (K - (I + 1)) + (J * NB + K - 1) * LDTB] =
              TB[TD + (K - I) + 1 + (J * NB + I - 1) * LDTB];
        }
      }

      if (J < NT - 1) {
        if (J > 0) {
          // Compute H(J,J)

          if (J == 1) {
            dgemm(
                'NoTranspose',
                'Transpose',
                KB,
                KB,
                KB,
                ONE,
                TB(TD + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
                LDTB - 1,
                A(J * NB + 1, (J - 1) * NB + 1),
                LDA,
                ZERO,
                WORK(J * NB + 1).asMatrix(N),
                N);
          } else {
            dgemm(
                'NoTranspose',
                'Transpose',
                KB,
                KB,
                NB + KB,
                ONE,
                TB(TD + NB + 1 + ((J - 1) * NB) * LDTB).asMatrix(LDTB - 1),
                LDTB - 1,
                A(J * NB + 1, (J - 2) * NB + 1),
                LDA,
                ZERO,
                WORK(J * NB + 1).asMatrix(N),
                N);
          }

          // Update with the previous column

          dgemm(
              'NoTranspose',
              'NoTranspose',
              N - (J + 1) * NB,
              NB,
              J * NB,
              -ONE,
              A((J + 1) * NB + 1, 1),
              LDA,
              WORK(NB + 1).asMatrix(N),
              N,
              ONE,
              A((J + 1) * NB + 1, J * NB + 1),
              LDA);
        }

        // Factorize panel

        dgetrf(N - (J + 1) * NB, NB, A((J + 1) * NB + 1, J * NB + 1), LDA,
            IPIV((J + 1) * NB + 1), IINFO);
        // IF (IINFO != 0 && INFO == 0) THEN
        //    INFO = IINFO+(J+1)*NB
        // END IF

        // Compute T(J+1, J), zero out for GEMM update

        KB = min(NB, N - (J + 1) * NB);
        dlaset('Full', KB, NB, ZERO, ZERO,
            TB(TD + NB + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1), LDTB - 1);
        dlacpy('Upper', KB, NB, A((J + 1) * NB + 1, J * NB + 1), LDA,
            TB(TD + NB + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1), LDTB - 1);
        if (J > 0) {
          dtrsm(
              'R',
              'L',
              'T',
              'U',
              KB,
              NB,
              ONE,
              A(J * NB + 1, (J - 1) * NB + 1),
              LDA,
              TB(TD + NB + 1 + (J * NB) * LDTB).asMatrix(LDTB - 1),
              LDTB - 1);
        }

        // Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
        // updates

        for (K = 1; K <= NB; K++) {
          for (I = 1; I <= KB; I++) {
            TB[TD - NB + K - I + 1 + (J * NB + NB + I - 1) * LDTB] =
                TB[TD + NB + I - K + 1 + (J * NB + K - 1) * LDTB];
          }
        }
        dlaset(
            'Upper', KB, NB, ZERO, ONE, A((J + 1) * NB + 1, J * NB + 1), LDA);

        // Apply pivots to trailing submatrix of A

        for (K = 1; K <= KB; K++) {
          // > Adjust ipiv
          IPIV[(J + 1) * NB + K] += (J + 1) * NB;

          I1 = (J + 1) * NB + K;
          I2 = IPIV[(J + 1) * NB + K];
          if (I1 != I2) {
            // > Apply pivots to previous columns of L
            dswap(K - 1, A(I1, (J + 1) * NB + 1).asArray(), LDA,
                A(I2, (J + 1) * NB + 1).asArray(), LDA);
            // > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
            if (I2 > (I1 + 1)) {
              dswap(I2 - I1 - 1, A(I1 + 1, I1).asArray(), 1,
                  A(I2, I1 + 1).asArray(), LDA);
            }
            // > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
            if (I2 < N) {
              dswap(N - I2, A(I2 + 1, I1).asArray(), 1, A(I2 + 1, I2).asArray(),
                  1);
            }
            // > Swap A(I1, I1) with A(I2, I2)
            PIV = A[I1][I1];
            A[I1][I1] = A[I2][I2];
            A[I2][I2] = PIV;
            // > Apply pivots to previous columns of L
            if (J > 0) {
              dswap(J * NB, A(I1, 1).asArray(), LDA, A(I2, 1).asArray(), LDA);
            }
          }
        }

        // Apply pivots to previous columns of L

        // CALL DLASWP( J*NB, A( 1, 1 ), LDA,
        // $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
      }
    }
  }

  // Factor the band matrix
  dgbtrf(N, N, NB, NB, TB.asMatrix(LDTB), LDTB, IPIV2, INFO);
}
