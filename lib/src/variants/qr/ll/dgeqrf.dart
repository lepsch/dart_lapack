import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqr2.dart';
import 'package:lapack/src/dlarfb.dart';
import 'package:lapack/src/dlarft.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/droundup_lwork.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgeqrf(
  final int M,
  final int N,
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
  bool LQUERY;
  int I, IB, IWS, J, K, LWKOPT, NB, NBMIN, NX, LBWORK, NT, LLWORK;
  final IINFO = Box(0);

  INFO.value = 0;
  NBMIN = 2;
  NX = 0;
  IWS = N;
  K = min(M, N);
  NB = ilaenv(1, 'DGEQRF', ' ', M, N, -1, -1);

  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'DGEQRF', ' ', M, N, -1, -1));
  }

  // Get NT, the size of the very last T, which is the left-over from in-between K-NX and K to K, eg.:
  //
  //        NB=3     2NB=6       K=10
  //        |        |           |
  //  1--2--3--4--5--6--7--8--9--10
  //              |     \________/
  //           K-NX=5      NT=4
  //
  // So here 4 x 4 is the last T stored in the workspace

  NT = K - ((K - NX) / NB).ceil() * NB;

  // optimal workspace = space for dlarfb + space for normal T's + space for the last T

  LLWORK = max(max((N - M) * K, (N - M) * NB), max(K * NB, NB * NB));
  LLWORK = (LLWORK / NB).ceil();

  if (K == 0) {
    LBWORK = 0;
    LWKOPT = 1;
    WORK[1] = LWKOPT.toDouble();
  } else if (NT > NB) {
    LBWORK = K - NT;

    // Optimal workspace for dlarfb = max(1,N)*NT

    LWKOPT = (LBWORK + LLWORK) * NB;
    WORK[1] = droundup_lwork(LWKOPT + NT * NT);
  } else {
    LBWORK = (K / NB).ceil() * NB;
    LWKOPT = (LBWORK + LLWORK - NB) * NB;
    WORK[1] = droundup_lwork(LWKOPT);
  }

  // Test the input arguments

  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (!LQUERY) {
    if (LWORK <= 0 || (M > 0 && LWORK < max(1, N))) INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DGEQRF', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (K == 0) {
    return;
  }

  if (NB > 1 && NB < K) {
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      if (NT <= NB) {
        IWS = (LBWORK + LLWORK - NB) * NB;
      } else {
        IWS = (LBWORK + LLWORK) * NB + NT * NT;
      }

      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        if (NT <= NB) {
          NB = LWORK ~/ (LLWORK + (LBWORK - NB));
        } else {
          NB = (LWORK - NT * NT) ~/ (LBWORK + LLWORK);
        }
        NBMIN = max(2, ilaenv(2, 'DGEQRF', ' ', M, N, -1, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code initially

    for (I = 1; I <= K - NX; I += NB) {
      IB = min(K - I + 1, NB);

      // Update the current column using old T's

      for (J = 1; J <= I - NB; J += NB) {
        // Apply H' to A(J:M,I:I+IB-1) from the left

        dlarfb(
            'Left',
            'Transpose',
            'Forward',
            'Columnwise',
            M - J + 1,
            IB,
            NB,
            A(J, J),
            LDA,
            WORK(J).asMatrix(LBWORK),
            LBWORK,
            A(J, I),
            LDA,
            WORK(LBWORK * NB + NT * NT + 1).asMatrix(IB),
            IB);
      }

      // Compute the QR factorization of the current block
      // A(I:M,I:I+IB-1)

      dgeqr2(M - I + 1, IB, A(I, I), LDA, TAU(I),
          WORK(LBWORK * NB + NT * NT + 1), IINFO);

      if (I + IB <= N) {
        // Form the triangular factor of the block reflector
        // H = H(i) H(i+1) . . . H(i+ib-1)

        dlarft('Forward', 'Columnwise', M - I + 1, IB, A(I, I), LDA, TAU(I),
            WORK(I).asMatrix(LBWORK), LBWORK);
      }
    }
  } else {
    I = 1;
  }

  // Use unblocked code to factor the last or only block.

  if (I <= K) {
    if (I != 1) {
      for (J = 1; J <= I - NB; J += NB) {
        // Apply H' to A(J:M,I:K) from the left

        dlarfb(
            'Left',
            'Transpose',
            'Forward',
            'Columnwise',
            M - J + 1,
            K - I + 1,
            NB,
            A(J, J),
            LDA,
            WORK(J).asMatrix(LBWORK),
            LBWORK,
            A(J, I),
            LDA,
            WORK(LBWORK * NB + NT * NT + 1).asMatrix(K - I + 1),
            K - I + 1);
      }
      dgeqr2(M - I + 1, K - I + 1, A(I, I), LDA, TAU(I),
          WORK(LBWORK * NB + NT * NT + 1), IINFO);
    } else {
      // Use unblocked code to factor the last or only block.

      dgeqr2(M - I + 1, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO);
    }
  }

  // Apply update to the column M+1:N when N > M

  if (M < N && I != 1) {
    // Form the last triangular factor of the block reflector
    // H = H(i) H(i+1) . . . H(i+ib-1)

    if (NT <= NB) {
      dlarft('Forward', 'Columnwise', M - I + 1, K - I + 1, A(I, I), LDA,
          TAU(I), WORK(I).asMatrix(LBWORK), LBWORK);
    } else {
      dlarft('Forward', 'Columnwise', M - I + 1, K - I + 1, A(I, I), LDA,
          TAU(I), WORK(LBWORK * NB + 1).asMatrix(NT), NT);
    }

    // Apply H' to A(1:M,M+1:N) from the left

    for (J = 1; J <= K - NX; J += NB) {
      IB = min(K - J + 1, NB);
      dlarfb(
          'Left',
          'Transpose',
          'Forward',
          'Columnwise',
          M - J + 1,
          N - M,
          IB,
          A(J, J),
          LDA,
          WORK(J).asMatrix(LBWORK),
          LBWORK,
          A(J, M + 1),
          LDA,
          WORK(LBWORK * NB + NT * NT + 1).asMatrix(N - M),
          N - M);
    }

    if (NT <= NB) {
      dlarfb(
          'Left',
          'Transpose',
          'Forward',
          'Columnwise',
          M - J + 1,
          N - M,
          K - J + 1,
          A(J, J),
          LDA,
          WORK(J).asMatrix(LBWORK),
          LBWORK,
          A(J, M + 1),
          LDA,
          WORK(LBWORK * NB + NT * NT + 1).asMatrix(N - M),
          N - M);
    } else {
      dlarfb(
          'Left',
          'Transpose',
          'Forward',
          'Columnwise',
          M - J + 1,
          N - M,
          K - J + 1,
          A(J, J),
          LDA,
          WORK(LBWORK * NB + 1).asMatrix(NT),
          NT,
          A(J, M + 1),
          LDA,
          WORK(LBWORK * NB + NT * NT + 1).asMatrix(N - M),
          N - M);
    }
  }

  WORK[1] = droundup_lwork(IWS);
}
