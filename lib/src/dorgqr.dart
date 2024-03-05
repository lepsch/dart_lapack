import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/dlarfb.dart';
import 'package:lapack/src/dlarft.dart';
import 'package:lapack/src/dorg2r.dart';

void dorgqr(
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  bool LQUERY;
  int I, IB, IWS, J, KI = 0, KK, L, LDWORK = 0, LWKOPT, NB, NBMIN, NX;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  NB = ilaenv(1, 'DORGQR', ' ', M, N, K, -1);
  LWKOPT = max(1, N) * NB;
  WORK[1] = LWKOPT.toDouble();
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0 || N > M) {
    INFO.value = -2;
  } else if (K < 0 || K > N) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('DORGQR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N <= 0) {
    WORK[1] = 1;
    return;
  }

  NBMIN = 2;
  NX = 0;
  IWS = N;
  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'DORGQR', ' ', M, N, K, -1));
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = N;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        NB = LWORK ~/ LDWORK;
        NBMIN = max(2, ilaenv(2, 'DORGQR', ' ', M, N, K, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code after the last block.
    // The first kk columns are handled by the block method.

    KI = ((K - NX - 1) ~/ NB) * NB;
    KK = min(K, KI + NB);

    // Set A(1:kk,kk+1:n) to zero.

    for (J = KK + 1; J <= N; J++) {
      for (I = 1; I <= KK; I++) {
        A[I][J] = ZERO;
      }
    }
  } else {
    KK = 0;
  }

  // Use unblocked code for the last or only block.

  if (KK < N) {
    dorg2r(M - KK, N - KK, K - KK, A(KK + 1, KK + 1), LDA, TAU(KK + 1), WORK,
        IINFO);
  }

  if (KK > 0) {
    // Use blocked code

    for (I = KI + 1; -NB < 0 ? I >= 1 : I <= 1; I += -NB) {
      IB = min(NB, K - I + 1);
      if (I + IB <= N) {
        // Form the triangular factor of the block reflector
        // H = H(i) H(i+1) . . . H(i+ib-1)

        dlarft('Forward', 'Columnwise', M - I + 1, IB, A(I, I), LDA, TAU(I),
            WORK.asMatrix(LDWORK), LDWORK);

        // Apply H to A(i:m,i+ib:n) from the left

        dlarfb(
            'Left',
            'No transpose',
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

      // Apply H to rows i:m of current block

      dorg2r(M - I + 1, IB, IB, A(I, I), LDA, TAU(I), WORK, IINFO);

      // Set rows 1:i-1 of current block to zero

      for (J = I; J <= I + IB - 1; J++) {
        for (L = 1; L <= I - 1; L++) {
          A[L][J] = ZERO;
        }
      }
    }
  }

  WORK[1] = IWS.toDouble();
}
