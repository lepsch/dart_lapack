import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeqr2p.dart';
import 'package:lapack/src/zlarfb.dart';
import 'package:lapack/src/zlarft.dart';

void zgeqrfp(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  bool LQUERY;
  int I, IB, IWS, K, LDWORK = 0, LWKMIN, LWKOPT, NB, NBMIN, NX;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  NB = ilaenv(1, 'ZGEQRF', ' ', M, N, -1, -1);
  K = min(M, N);
  if (K == 0) {
    LWKMIN = 1;
    LWKOPT = 1;
  } else {
    LWKMIN = N;
    LWKOPT = N * NB;
  }
  WORK[1] = LWKOPT.toComplex();

  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (LWORK < LWKMIN && !LQUERY) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZGEQRFP', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (K == 0) {
    WORK[1] = Complex.one;
    return;
  }

  NBMIN = 2;
  NX = 0;
  IWS = LWKMIN;
  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'ZGEQRF', ' ', M, N, -1, -1));
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = N;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        NB = LWORK ~/ LDWORK;
        NBMIN = max(2, ilaenv(2, 'ZGEQRF', ' ', M, N, -1, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code initially

    for (I = 1; NB < 0 ? I >= K - NX : I <= K - NX; I += NB) {
      // 10
      IB = min(K - I + 1, NB);

      // Compute the QR factorization of the current block
      // A(i:m,i:i+ib-1)

      zgeqr2p(M - I + 1, IB, A(I, I), LDA, TAU(I), WORK, IINFO);
      if (I + IB <= N) {
        // Form the triangular factor of the block reflector
        // H = H(i) H(i+1) . . . H(i+ib-1)

        zlarft('Forward', 'Columnwise', M - I + 1, IB, A(I, I), LDA, TAU(I),
            WORK.asMatrix(LDWORK), LDWORK);

        // Apply H**H to A(i:m,i+ib:n) from the left

        zlarfb(
            'Left',
            'Conjugate transpose',
            'Forward',
            'Columnwise',
            M - I + 1,
            N - I - IB + 1,
            IB,
            A(I, I),
            LDA,
            WORK.asMatrix(LDWORK),
            LDWORK,
            A(I, I + IB),
            LDA,
            WORK(IB + 1).asMatrix(LDWORK),
            LDWORK);
      }
    } // 10
  } else {
    I = 1;
  }

  // Use unblocked code to factor the last or only block.

  if (I <= K) zgeqr2p(M - I + 1, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO);

  WORK[1] = IWS.toComplex();
}
