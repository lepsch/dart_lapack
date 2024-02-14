import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgelq2.dart';
import 'package:lapack/src/dlarfb.dart';
import 'package:lapack/src/dlarft.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgelqf(
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
  final A = A_.dim(LDA);
  final TAU = TAU_.dim();
  final WORK = WORK_.dim();
  bool LQUERY;
  int I, IB, IWS, K, LDWORK = 0, LWKOPT, NB, NBMIN, NX;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  K = min(M, N);
  NB = ilaenv(1, 'DGELQF', ' ', M, N, -1, -1);
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (!LQUERY) {
    if (LWORK <= 0 || (N > 0 && LWORK < max(1, M))) INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DGELQF', -INFO.value);
    return;
  } else if (LQUERY) {
    if (K == 0) {
      LWKOPT = 1;
    } else {
      LWKOPT = M * NB;
    }
    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // Quick return if possible

  if (K == 0) {
    WORK[1] = 1;
    return;
  }

  NBMIN = 2;
  NX = 0;
  IWS = M;
  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'DGELQF', ' ', M, N, -1, -1));
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = M;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        NB = LWORK ~/ LDWORK;
        NBMIN = max(2, ilaenv(2, 'DGELQF', ' ', M, N, -1, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code initially

    for (I = 1; NB < 0 ? I >= K - NX : I <= K - NX; I += NB) {

      IB = min(K - I + 1, NB);

      // Compute the LQ factorization of the current block
      // A(i:i+ib-1,i:n)

      dgelq2(IB, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO);
      if (I + IB <= M) {
        // Form the triangular factor of the block reflector
        // H = H(i) H(i+1) . . . H(i+ib-1)

        dlarft('Forward', 'Rowwise', N - I + 1, IB, A(I, I), LDA, TAU(I), WORK.asMatrix(LDWORK),
            LDWORK);

        // Apply H to A(i+ib:m,i:n) from the right

        dlarfb(
            'Right',
            'No transpose',
            'Forward',
            'Rowwise',
            M - I - IB + 1,
            N - I + 1,
            IB,
            A(I, I),
            LDA,
            WORK.asMatrix(LDWORK),
            LDWORK,
            A(I + IB, I),
            LDA,
            WORK(IB + 1).asMatrix(LDWORK),
            LDWORK);
      }
    }
  } else {
    I = 1;
  }

  // Use unblocked code to factor the last or only block.

  if (I <= K) dgelq2(M - I + 1, N - I + 1, A(I, I), LDA, TAU(I), WORK, IINFO);

  WORK[1] = IWS.toDouble();
}
