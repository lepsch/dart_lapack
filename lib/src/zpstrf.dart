import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zpstf2.dart';

void zpstrf(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> PIV_,
  final Box<int> RANK,
  final double TOL,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final PIV = PIV_.having(length: N);
  final WORK = WORK_.having(length: 2 * N);

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  Complex ZTEMP;
  double AJJ = 0, DSTOP, DTEMP;
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
    xerbla('ZPSTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Get block size

  NB = ilaenv(1, 'ZPOTRF', UPLO, N, -1, -1, -1);
  if (NB <= 1 || NB >= N) {
    // Use unblocked code

    zpstf2(UPLO, N, A(1, 1), LDA, PIV, RANK, TOL, WORK, INFO);
    return;
  }

  // Initialize PIV

  for (I = 1; I <= N; I++) {
    PIV[I] = I;
  }

  // Compute stopping value

  for (I = 1; I <= N; I++) {
    WORK[I] = A[I][I].real;
  }
  PVT = WORK.maxloc(1, N);
  AJJ = A[PVT][PVT].real;
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
    // Compute the Cholesky factorization P**T * A * P = U**H * U

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
            WORK[I] =
                WORK[I] + (A[J - 1][I].conjugate() * A[J - 1][I]).real;
          }
          WORK[N + I] = A[I][I].real - WORK[I];
        }

        if (J > 1) {
          ITEMP = WORK.maxloc(N + J, 2 * N);
          PVT = ITEMP + J - 1;
          AJJ = WORK[N + PVT];
          if (AJJ <= DSTOP || disnan(AJJ)) {
            A[J][J] = AJJ.toComplex();
            // Rank is the number of steps completed.  Set INFO.value = 1 to signal
            // that the factorization cannot be used to solve a system.

            RANK.value = J - 1;
            INFO.value = 1;
            return;
          }
        }

        if (J != PVT) {
          // Pivot OK, so can now swap pivot rows and columns

          A[PVT][PVT] = A[J][J];
          zswap(J - 1, A(1, J).asArray(), 1, A(1, PVT).asArray(), 1);
          if (PVT < N) {
            zswap(N - PVT, A(J, PVT + 1).asArray(), LDA,
                A(PVT, PVT + 1).asArray(), LDA);
          }
          for (I = J + 1; I <= PVT - 1; I++) {
            ZTEMP = A[J][I].conjugate();
            A[J][I] = A[I][PVT].conjugate();
            A[I][PVT] = ZTEMP;
          }
          A[J][PVT] = A[J][PVT].conjugate();

          // Swap dot products and PIV

          DTEMP = WORK[J];
          WORK[J] = WORK[PVT];
          WORK[PVT] = DTEMP;
          ITEMP = PIV[PVT];
          PIV[PVT] = PIV[J];
          PIV[J] = ITEMP;
        }

        AJJ = sqrt(AJJ);
        A[J][J] = AJJ.toComplex();

        // Compute elements J+1:N of row J.

        if (J < N) {
          zlacgv(J - 1, A(1, J).asArray(), 1);
          zgemv('Trans', J - K, N - J, -Complex.one, A(K, J + 1), LDA,
              A(K, J).asArray(), 1, Complex.one, A(J, J + 1).asArray(), LDA);
          zlacgv(J - 1, A(1, J).asArray(), 1);
          zdscal(N - J, ONE / AJJ, A(J, J + 1).asArray(), LDA);
        }
      }

      // Update trailing matrix, J already incremented

      if (K + JB <= N) {
        zherk('Upper', 'Conj Trans', N - J + 1, JB, -ONE, A(K, J), LDA, ONE,
            A(J, J), LDA);
      }
    }
  } else {
    // Compute the Cholesky factorization P**T * A * P = L * L**H

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
            WORK[I] =
                WORK[I] + (A[I][J - 1].conjugate() * A[I][J - 1]).real;
          }
          WORK[N + I] = A[I][I].real - WORK[I];
        }

        if (J > 1) {
          ITEMP = WORK.maxloc(N + J, 2 * N);
          PVT = ITEMP + J - 1;
          AJJ = WORK[N + PVT];
          if (AJJ <= DSTOP || disnan(AJJ)) {
            A[J][J] = AJJ.toComplex();
            // Rank is the number of steps completed.  Set INFO.value = 1 to signal
            // that the factorization cannot be used to solve a system.

            RANK.value = J - 1;
            INFO.value = 1;
            return;
          }
        }

        if (J != PVT) {
          // Pivot OK, so can now swap pivot rows and columns

          A[PVT][PVT] = A[J][J];
          zswap(J - 1, A(J, 1).asArray(), LDA, A(PVT, 1).asArray(), LDA);
          if (PVT < N) {
            zswap(N - PVT, A(PVT + 1, J).asArray(), 1,
                A(PVT + 1, PVT).asArray(), 1);
          }
          for (I = J + 1; I <= PVT - 1; I++) {
            ZTEMP = A[I][J].conjugate();
            A[I][J] = A[PVT][I].conjugate();
            A[PVT][I] = ZTEMP;
          }
          A[PVT][J] = A[PVT][J].conjugate();

          // Swap dot products and PIV

          DTEMP = WORK[J];
          WORK[J] = WORK[PVT];
          WORK[PVT] = DTEMP;
          ITEMP = PIV[PVT];
          PIV[PVT] = PIV[J];
          PIV[J] = ITEMP;
        }

        AJJ = sqrt(AJJ);
        A[J][J] = AJJ.toComplex();

        // Compute elements J+1:N of column J.

        if (J < N) {
          zlacgv(J - 1, A(J, 1).asArray(), LDA);
          zgemv('No Trans', N - J, J - K, -Complex.one, A(J + 1, K), LDA,
              A(J, K).asArray(), LDA, Complex.one, A(J + 1, J).asArray(), 1);
          zlacgv(J - 1, A(J, 1).asArray(), LDA);
          zdscal(N - J, ONE / AJJ, A(J + 1, J).asArray(), 1);
        }
      }

      // Update trailing matrix, J already incremented

      if (K + JB <= N) {
        zherk('Lower', 'No Trans', N - J + 1, JB, -ONE, A(J, K), LDA, ONE,
            A(J, J), LDA);
      }
    }
  }

  // Ran to completion, A has full rank

  RANK.value = N;
}
