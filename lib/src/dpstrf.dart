import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dpstf2.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpstrf(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> PIV_,
  final Box<int> RANK,
  final double TOL,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final PIV = PIV_.having(length: N);
  final WORK = WORK_.having(length: 2 * N);
  const ONE = 1.0, ZERO = 0.0;
  double AJJ, DSTOP, DTEMP;
  int I, ITEMP, J, JB, K, NB, PVT;
  bool UPPER;

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
    xerbla('DPSTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Get block size

  NB = ilaenv(1, 'DPOTRF', UPLO, N, -1, -1, -1);
  if (NB <= 1 || NB >= N) {
    // Use unblocked code

    dpstf2(UPLO, N, A(1, 1), LDA, PIV, RANK, TOL, WORK, INFO);
    return;
  } else {
    // Initialize PIV

    for (I = 1; I <= N; I++) {
      PIV[I] = I;
    }

    // Compute stopping value

    PVT = 1;
    AJJ = A[PVT][PVT];
    for (I = 2; I <= N; I++) {
      if (A[I][I] > AJJ) {
        PVT = I;
        AJJ = A[PVT][PVT];
      }
    }
    if (AJJ <= ZERO || disnan(AJJ)) {
      RANK.value = 0;
      INFO.value = 1;
      return;
    }

    // Compute stopping value if not supplied

    if (TOL < ZERO) {
      DSTOP = N * dlamch('Epsilon') * AJJ;
    } else {
      DSTOP = TOL;
    }

    if (UPPER) {
      // Compute the Cholesky factorization P**T * A * P = U**T * U

      for (K = 1; K <= N; K += NB) {
        // Account for last block not being NB wide

        JB = min(NB, N - K + 1);

        // Set relevant part of first half of WORK to zero,
        // holds dot products

        for (I = K; I <= N; I++) {
          WORK[I] = 0;
        }

        for (J = K; J <= K + JB - 1; J++) {
          // Find pivot, test for exit, else swap rows and columns
          // Update dot products, compute possible pivots which are
          // stored in the second half of WORK

          for (I = J; I <= N; I++) {
            if (J > K) {
              WORK[I] += pow(A[J - 1][I], 2);
            }
            WORK[N + I] = A[I][I] - WORK[I];
          }

          if (J > 1) {
            ITEMP = WORK.maxloc(N + J, 2 * N, dim: 1);
            PVT = ITEMP + J - 1;
            AJJ = WORK[N + PVT];
            if (AJJ <= DSTOP || disnan(AJJ)) {
              A[J][J] = AJJ;

              // Rank is the number of steps completed.  Set INFO = 1 to signal
              // that the factorization cannot be used to solve a system.

              RANK.value = J - 1;
              INFO.value = 1;
              return;
            }
          }

          if (J != PVT) {
            // Pivot OK, so can now swap pivot rows and columns

            A[PVT][PVT] = A[J][J];
            dswap(J - 1, A(1, J).asArray(), 1, A(1, PVT).asArray(), 1);
            if (PVT < N) {
              dswap(N - PVT, A(J, PVT + 1).asArray(), LDA,
                  A(PVT, PVT + 1).asArray(), LDA);
            }
            dswap(PVT - J - 1, A(J, J + 1).asArray(), LDA,
                A(J + 1, PVT).asArray(), 1);

            // Swap dot products and PIV

            DTEMP = WORK[J];
            WORK[J] = WORK[PVT];
            WORK[PVT] = DTEMP;
            ITEMP = PIV[PVT];
            PIV[PVT] = PIV[J];
            PIV[J] = ITEMP;
          }

          AJJ = sqrt(AJJ);
          A[J][J] = AJJ;

          // Compute elements J+1:N of row J.

          if (J < N) {
            dgemv('Trans', J - K, N - J, -ONE, A(K, J + 1), LDA,
                A(K, J).asArray(), 1, ONE, A(J, J + 1).asArray(), LDA);
            dscal(N - J, ONE / AJJ, A(J, J + 1).asArray(), LDA);
          }
        }

        // Update trailing matrix, J already incremented

        if (K + JB <= N) {
          dsyrk('Upper', 'Trans', N - J + 1, JB, -ONE, A(K, J), LDA, ONE,
              A(J, J), LDA);
        }
      }
    } else {
      // Compute the Cholesky factorization P**T * A * P = L * L**T

      for (K = 1; K <= N; K += NB) {
        // Account for last block not being NB wide

        JB = min(NB, N - K + 1);

        // Set relevant part of first half of WORK to zero,
        // holds dot products

        for (I = K; I <= N; I++) {
          WORK[I] = 0;
        }

        for (J = K; J <= K + JB - 1; J++) {
          // Find pivot, test for exit, else swap rows and columns
          // Update dot products, compute possible pivots which are
          // stored in the second half of WORK

          for (I = J; I <= N; I++) {
            if (J > K) {
              WORK[I] += pow(A[I][J - 1], 2);
            }
            WORK[N + I] = A[I][I] - WORK[I];
          }

          if (J > 1) {
            ITEMP = WORK.maxloc(N + J, 2 * N, dim: 1);
            PVT = ITEMP + J - 1;
            AJJ = WORK[N + PVT];
            if (AJJ <= DSTOP || disnan(AJJ)) {
              A[J][J] = AJJ;

              // Rank is the number of steps completed.  Set INFO = 1 to signal
              // that the factorization cannot be used to solve a system.

              RANK.value = J - 1;
              INFO.value = 1;
              return;
            }
          }

          if (J != PVT) {
            // Pivot OK, so can now swap pivot rows and columns

            A[PVT][PVT] = A[J][J];
            dswap(J - 1, A(J, 1).asArray(), LDA, A(PVT, 1).asArray(), LDA);
            if (PVT < N) {
              dswap(N - PVT, A(PVT + 1, J).asArray(), 1,
                  A(PVT + 1, PVT).asArray(), 1);
            }
            dswap(PVT - J - 1, A(J + 1, J).asArray(), 1,
                A(PVT, J + 1).asArray(), LDA);

            // Swap dot products and PIV

            DTEMP = WORK[J];
            WORK[J] = WORK[PVT];
            WORK[PVT] = DTEMP;
            ITEMP = PIV[PVT];
            PIV[PVT] = PIV[J];
            PIV[J] = ITEMP;
          }

          AJJ = sqrt(AJJ);
          A[J][J] = AJJ;

          // Compute elements J+1:N of column J.

          if (J < N) {
            dgemv('No Trans', N - J, J - K, -ONE, A(J + 1, K), LDA,
                A(J, K).asArray(), LDA, ONE, A(J + 1, J).asArray(), 1);
            dscal(N - J, ONE / AJJ, A(J + 1, J).asArray(), 1);
          }
        }

        // Update trailing matrix, J already incremented

        if (K + JB <= N) {
          dsyrk('Lower', 'No Trans', N - J + 1, JB, -ONE, A(J, K), LDA, ONE,
              A(J, J), LDA);
        }
      }
    }
  }

  // Ran to completion, A has full rank
  // Rank is the number of steps completed.

  RANK.value = N;
}
