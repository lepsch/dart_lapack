import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/nint.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlarfg.dart';

void zlaqps(
  final int M,
  final int N,
  final int OFFSET,
  final int NB,
  final Box<int> KB,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> JPVT_,
  final Array<Complex> TAU_,
  final Array<double> VN1_,
  final Array<double> VN2_,
  final Array<Complex> AUXV_,
  final Matrix<Complex> F_,
  final int LDF,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final F = F_.having(ld: LDF);
  final JPVT = JPVT_.having();
  final TAU = TAU_.having();
  final VN1 = VN1_.having();
  final VN2 = VN2_.having();
  final AUXV = AUXV_.having();
  const ZERO = 0.0, ONE = 1.0;
  int ITEMP, J, K, LASTRK, LSTICC, PVT, RK;
  double TEMP, TEMP2, TOL3Z;
  Complex AKK;

  LASTRK = min(M, N + OFFSET);
  LSTICC = 0;
  K = 0;
  TOL3Z = sqrt(dlamch('Epsilon'));

  // Beginning of while loop.

  while ((K < NB) && (LSTICC == 0)) {
    K++;
    RK = OFFSET + K;

    // Determine ith pivot column and swap if necessary

    PVT = (K - 1) + idamax(N - K + 1, VN1(K), 1);
    if (PVT != K) {
      zswap(M, A(1, PVT).asArray(), 1, A(1, K).asArray(), 1);
      zswap(K - 1, F(PVT, 1).asArray(), LDF, F(K, 1).asArray(), LDF);
      ITEMP = JPVT[PVT];
      JPVT[PVT] = JPVT[K];
      JPVT[K] = ITEMP;
      VN1[PVT] = VN1[K];
      VN2[PVT] = VN2[K];
    }

    // Apply previous Householder reflectors to column K:
    // A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**H.

    if (K > 1) {
      for (J = 1; J <= K - 1; J++) {
        // 20
        F[K][J] = F[K][J].conjugate();
      } // 20
      zgemv('No transpose', M - RK + 1, K - 1, -Complex.one, A(RK, 1), LDA,
          F(K, 1).asArray(), LDF, Complex.one, A(RK, K).asArray(), 1);
      for (J = 1; J <= K - 1; J++) {
        // 30
        F[K][J] = F[K][J].conjugate();
      } // 30
    }

    // Generate elementary reflector H(k).

    if (RK < M) {
      zlarfg(M - RK + 1, A(RK, K), A(RK + 1, K).asArray(), 1, TAU(K));
    } else {
      zlarfg(1, A(RK, K), A(RK, K).asArray(), 1, TAU(K));
    }

    AKK = A[RK][K];
    A[RK][K] = Complex.one;

    // Compute Kth column of F:

    // Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).

    if (K < N) {
      zgemv('Conjugate transpose', M - RK + 1, N - K, TAU[K], A(RK, K + 1), LDA,
          A(RK, K).asArray(), 1, Complex.zero, F(K + 1, K).asArray(), 1);
    }

    // Padding F(1:K,K) with zeros.

    for (J = 1; J <= K; J++) {
      // 40
      F[J][K] = Complex.zero;
    } // 40

    // Incremental updating of F:
    // F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**H
    //             *A(RK:M,K).

    if (K > 1) {
      zgemv('Conjugate transpose', M - RK + 1, K - 1, -TAU[K], A(RK, 1), LDA,
          A(RK, K).asArray(), 1, Complex.zero, AUXV(1), 1);

      zgemv('No transpose', N, K - 1, Complex.one, F(1, 1), LDF, AUXV(1), 1,
          Complex.one, F(1, K).asArray(), 1);
    }

    // Update the current row of A:
    // A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H.

    if (K < N) {
      zgemm('No transpose', 'Conjugate transpose', 1, N - K, K, -Complex.one,
          A(RK, 1), LDA, F(K + 1, 1), LDF, Complex.one, A(RK, K + 1), LDA);
    }

    // Update partial column norms.

    if (RK < LASTRK) {
      for (J = K + 1; J <= N; J++) {
        // 50
        if (VN1[J] != ZERO) {
          // NOTE: The following 4 lines follow from the analysis in
          // Lapack Working Note 176.

          TEMP = A[RK][J].abs() / VN1[J];
          TEMP = max(ZERO, (ONE + TEMP) * (ONE - TEMP));
          TEMP2 = TEMP * pow(VN1[J] / VN2[J], 2);
          if (TEMP2 <= TOL3Z) {
            VN2[J] = LSTICC.toDouble();
            LSTICC = J;
          } else {
            VN1[J] = VN1[J] * sqrt(TEMP);
          }
        }
      } // 50
    }

    A[RK][K] = AKK;
  }
  KB.value = K;
  RK = OFFSET + KB.value;

  // Apply the block reflector to the rest of the matrix:
  // A(OFFSET+KB.value+1:M,KB.value+1:N) := A(OFFSET+KB.value+1:M,KB.value+1:N) -
  //                     A(OFFSET+KB.value+1:M,1:KB.value)*F(KB.value+1:N,1:KB.value)**H.

  if (KB.value < min(N, M - OFFSET)) {
    zgemm(
        'No transpose',
        'Conjugate transpose',
        M - RK,
        N - KB.value,
        KB.value,
        -Complex.one,
        A(RK + 1, 1),
        LDA,
        F(KB.value + 1, 1),
        LDF,
        Complex.one,
        A(RK + 1, KB.value + 1),
        LDA);
  }

  // Recomputation of difficult columns.

  while (LSTICC > 0) {
    ITEMP = nint(VN2[LSTICC]);
    VN1[LSTICC] = dznrm2(M - RK, A(RK + 1, LSTICC).asArray(), 1);

    // NOTE: The computation of VN1[ LSTICC ] relies on the fact that
    // SNRM2 does not fail on vectors with norm below the value of
    // sqrt(dlamch('S'))

    VN2[LSTICC] = VN1[LSTICC];
    LSTICC = ITEMP;
  }
}
