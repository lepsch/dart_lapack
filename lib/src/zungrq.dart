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
import 'package:lapack/src/zungr2.dart';

void zungrq(
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
  int I, IB, II, IWS, J, KK, L, LDWORK = 0, LWKOPT, NB = 0, NBMIN, NX;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < M) {
    INFO.value = -2;
  } else if (K < 0 || K > M) {
    INFO.value = -3;
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  }

  if (INFO.value == 0) {
    if (M <= 0) {
      LWKOPT = 1;
    } else {
      NB = ilaenv(1, 'ZUNGRQ', ' ', M, N, K, -1);
      LWKOPT = M * NB;
    }
    WORK[1] = LWKOPT.toComplex();

    if (LWORK < max(1, M) && !LQUERY) {
      INFO.value = -8;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZUNGRQ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M <= 0) {
    return;
  }

  NBMIN = 2;
  NX = 0;
  IWS = M;
  if (NB > 1 && NB < K) {
    // Determine when to cross over from blocked to unblocked code.

    NX = max(0, ilaenv(3, 'ZUNGRQ', ' ', M, N, K, -1));
    if (NX < K) {
      // Determine if workspace is large enough for blocked code.

      LDWORK = M;
      IWS = LDWORK * NB;
      if (LWORK < IWS) {
        // Not enough workspace to use optimal NB:  reduce NB and
        // determine the minimum value of NB.

        NB = LWORK ~/ LDWORK;
        NBMIN = max(2, ilaenv(2, 'ZUNGRQ', ' ', M, N, K, -1));
      }
    }
  }

  if (NB >= NBMIN && NB < K && NX < K) {
    // Use blocked code after the first block.
    // The last kk rows are handled by the block method.

    KK = min(K, ((K - NX + NB - 1) ~/ NB) * NB);

    // Set A(1:m-kk,n-kk+1:n) to zero.

    for (J = N - KK + 1; J <= N; J++) {
      for (I = 1; I <= M - KK; I++) {
        A[I][J] = Complex.zero;
      }
    }
  } else {
    KK = 0;
  }

  // Use unblocked code for the first or only block.

  zungr2(M - KK, N - KK, K - KK, A, LDA, TAU, WORK, IINFO);

  if (KK > 0) {
    // Use blocked code

    for (I = K - KK + 1; I <= K; I += NB) {
      IB = min(NB, K - I + 1);
      II = M - K + I;
      if (II > 1) {
        // Form the triangular factor of the block reflector
        // H = H(i+ib-1) . . . H(i+1) H(i)

        zlarft('Backward', 'Rowwise', N - K + I + IB - 1, IB, A(II, 1), LDA,
            TAU(I), WORK.asMatrix(), LDWORK);

        // Apply H**H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right

        zlarfb(
            'Right',
            'Conjugate transpose',
            'Backward',
            'Rowwise',
            II - 1,
            N - K + I + IB - 1,
            IB,
            A(II, 1),
            LDA,
            WORK.asMatrix(),
            LDWORK,
            A,
            LDA,
            WORK(IB + 1).asMatrix(),
            LDWORK);
      }

      // Apply H**H to columns 1:n-k+i+ib-1 of current block

      zungr2(IB, N - K + I + IB - 1, IB, A(II, 1), LDA, TAU(I), WORK, IINFO);

      // Set columns n-k+i+ib:n of current block to zero

      for (L = N - K + I + IB; L <= N; L++) {
        for (J = II; J <= II + IB - 1; J++) {
          A[J][L] = Complex.zero;
        }
      }
    }
  }

  WORK[1] = IWS.toComplex();
}
