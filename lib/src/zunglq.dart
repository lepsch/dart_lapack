// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarfb.dart';
import 'package:lapack/src/zlarft.dart';
import 'package:lapack/src/zungl2.dart';

void zunglq(
  final int M,
  final int N,
  final int K,
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
  int I, IB, IWS, J, KI = 0, KK, L, LDWORK = 0, LWKOPT, NB, NBMIN, NX;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  NB = ilaenv(1, 'ZUNGLQ', ' ', M, N, K, -1);
  LWKOPT = max(1, M) * NB;
  WORK[1] = LWKOPT.toComplex();
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < M) {
    INFO.value = -2;
  } else if (K < 0 || K > M) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LWORK < max(1, M) && !LQUERY) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZUNGLQ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M <= 0) {
    WORK[1] = Complex.one;
    return;
  }

  NBMIN = 2;
  NX = 0;
  IWS = M;
  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'ZUNGLQ', ' ', M, N, K, -1));
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = M;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        NB = LWORK ~/ LDWORK;
        NBMIN = max(2, ilaenv(2, 'ZUNGLQ', ' ', M, N, K, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code after the last block.
    // The first kk rows are handled by the block method.

    KI = ((K - NX - 1) ~/ NB) * NB;
    KK = min(K, KI + NB);

    // Set A(kk+1:m,1:kk) to zero.

    for (J = 1; J <= KK; J++) {
      for (I = KK + 1; I <= M; I++) {
        A[I][J] = Complex.zero;
      }
    }
  } else {
    KK = 0;
  }

  // Use unblocked code for the last or only block.

  if (KK < M) {
    zungl2(M - KK, N - KK, K - KK, A(KK + 1, KK + 1), LDA, TAU(KK + 1), WORK,
        IINFO);
  }

  if (KK > 0) {
    // Use blocked code

    for (I = KI + 1; I >= 1; I -= NB) {
      IB = min(NB, K - I + 1);
      if (I + IB <= M) {
        // Form the triangular factor of the block reflector
        // H = H(i) H(i+1) . . . H(i+ib-1)

        zlarft('Forward', 'Rowwise', N - I + 1, IB, A(I, I), LDA, TAU(I),
            WORK.asMatrix(LDWORK), LDWORK);

        // Apply H**H to A(i+ib:m,i:n) from the right

        zlarfb(
            'Right',
            'Conjugate transpose',
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

      // Apply H**H to columns i:n of current block

      zungl2(IB, N - I + 1, IB, A(I, I), LDA, TAU(I), WORK, IINFO);

      // Set columns 1:i-1 of current block to zero

      for (J = 1; J <= I - 1; J++) {
        for (L = I; L <= I + IB - 1; L++) {
          A[L][J] = Complex.zero;
        }
      }
    }
  }

  WORK[1] = IWS.toComplex();
}
