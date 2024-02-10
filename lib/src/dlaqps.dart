import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaqps(
  final int M,
  final int N,
  final int OFFSET,
  final int NB,
  final Box<int> KB,
  final Matrix<double> A,
  final int LDA,
  final Array<int> JPVT,
  final Array<double> TAU,
  final Array<double> VN1,
  final Array<double> VN2,
  final Array<double> AUXV,
  final Matrix<double> F,
  final int LDF,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int ITEMP, J, K, LASTRK, LSTICC, PVT, RK;
  double AKK, TEMP, TEMP2, TOL3Z;

  LASTRK = min(M, N + OFFSET);
  LSTICC = 0;
  K = 0;
  TOL3Z = sqrt(dlamch('Epsilon'));

  // Beginning of while loop.

  while ((K < NB) && (LSTICC == 0)) {
    K = K + 1;
    RK = OFFSET + K;

    // Determine ith pivot column and swap if necessary

    PVT = (K - 1) + idamax(N - K + 1, VN1(K), 1);
    if (PVT != K) {
      dswap(M, A(1, PVT).asArray(), 1, A(1, K).asArray(), 1);
      dswap(K - 1, F(PVT, 1).asArray(), LDF, F(K, 1).asArray(), LDF);
      ITEMP = JPVT[PVT];
      JPVT[PVT] = JPVT[K];
      JPVT[K] = ITEMP;
      VN1[PVT] = VN1[K];
      VN2[PVT] = VN2[K];
    }

    // Apply previous Householder reflectors to column K:
    // A[RK:M][K] := A[RK:M][K] - A[RK:M][1:K-1]*F[K,1:K-1]**T.

    if (K > 1) {
      dgemv('No transpose', M - RK + 1, K - 1, -ONE, A(RK, 1), LDA,
          F(K, 1).asArray(), LDF, ONE, A(RK, K).asArray(), 1);
    }

    // Generate elementary reflector H(k).

    if (RK < M) {
      dlarfg(M - RK + 1, A.box(RK, K), A(RK + 1, K).asArray(), 1, TAU.box(K));
    } else {
      dlarfg(1, A.box(RK, K), A(RK, K).asArray(), 1, TAU.box(K));
    }

    AKK = A[RK][K];
    A[RK][K] = ONE;

    // Compute Kth column of F:

    // Compute  F[K+1:N,K] := tau[K]*A[RK:M][K+1:N]**T*A[RK:M][K].

    if (K < N) {
      dgemv('Transpose', M - RK + 1, N - K, TAU[K], A(RK, K + 1), LDA,
          A(RK, K).asArray(), 1, ZERO, F(K + 1, K).asArray(), 1);
    }

    // Padding F[1:K,K] with zeros.

    for (J = 1; J <= K; J++) {
      F[J][K] = ZERO;
    }

    // Incremental updating of F:
    // F[1:N,K] := F[1:N,K] - tau[K]*F[1:N,1:K-1]*A[RK:M][1:K-1]**T
    // *A[RK:M][K].

    if (K > 1) {
      dgemv('Transpose', M - RK + 1, K - 1, -TAU[K], A(RK, 1), LDA,
          A(RK, K).asArray(), 1, ZERO, AUXV(1), 1);

      dgemv('No transpose', N, K - 1, ONE, F(1, 1), LDF, AUXV(1), 1, ONE,
          F(1, K).asArray(), 1);
    }

    // Update the current row of A:
    // A[RK][K+1:N] := A[RK][K+1:N] - A[RK][1:K]*F[K+1:N,1:K]**T.

    if (K < N) {
      dgemv('No transpose', N - K, K, -ONE, F(K + 1, 1), LDF,
          A(RK, 1).asArray(), LDA, ONE, A(RK, K + 1).asArray(), LDA);
    }

    // Update partial column norms.

    if (RK < LASTRK) {
      for (J = K + 1; J <= N; J++) {
        if (VN1[J] != ZERO) {
          // NOTE: The following 4 lines follow from the analysis in
          // Lapack Working Note 176.

          TEMP = (A[RK][J]).abs() / VN1[J];
          TEMP = max(ZERO, (ONE + TEMP) * (ONE - TEMP));
          TEMP2 = TEMP * (VN1[J] / pow((VN2[J]), 2));
          if (TEMP2 <= TOL3Z) {
            VN2[J] = LSTICC.toDouble();
            LSTICC = J;
          } else {
            VN1[J] = VN1[J] * sqrt(TEMP);
          }
        }
      }
    }

    A[RK][K] = AKK;
  } // End of while loop.
  KB.value = K;
  RK = OFFSET + KB.value;

  // Apply the block reflector to the rest of the matrix:
  // A[OFFSET+KB.value+1:M][KB.value+1:N] := A[OFFSET+KB.value+1:M][KB.value+1:N] -
  // A[OFFSET+KB.value+1:M][1:KB.value]*F[KB.value+1:N,1:KB.value]**T.

  if (KB.value < min(N, M - OFFSET)) {
    dgemm(
        'No transpose',
        'Transpose',
        M - RK,
        N - KB.value,
        KB.value,
        -ONE,
        A(RK + 1, 1),
        LDA,
        F(KB.value + 1, 1),
        LDF,
        ONE,
        A(RK + 1, KB.value + 1),
        LDA);
  }

  // Recomputation of difficult columns.

  while (LSTICC > 0) {
    ITEMP = VN2[LSTICC].round();
    VN1[LSTICC] = dnrm2(M - RK, A(RK + 1, LSTICC).asArray(), 1);

    // NOTE: The computation of VN1[ LSTICC ] relies on the fact that
    // SNRM2 does not fail on vectors with norm below the value of
    // sqrt(dlamch('S'))

    VN2[LSTICC] = VN1[LSTICC];
    LSTICC = ITEMP;
  }
}
