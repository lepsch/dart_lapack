// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zgemv.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/blas/ztrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/ztrtri.dart';

void zgetri(
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  bool LQUERY;
  int I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, NBMIN, NN;

  // Test the input parameters.

  INFO.value = 0;
  NB = ilaenv(1, 'ZGETRI', ' ', N, -1, -1, -1);
  LWKOPT = max(1, N * NB);
  WORK[1] = LWKOPT.toComplex();
  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -1;
  } else if (LDA < max(1, N)) {
    INFO.value = -3;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZGETRI', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
  // and the inverse is not computed.

  ztrtri('Upper', 'Non-unit', N, A, LDA, INFO);
  if (INFO.value > 0) return;

  NBMIN = 2;
  LDWORK = N;
  if (NB > 1 && NB < N) {
    IWS = max(LDWORK * NB, 1);
    if (LWORK < IWS) {
      NB = LWORK ~/ LDWORK;
      NBMIN = max(2, ilaenv(2, 'ZGETRI', ' ', N, -1, -1, -1));
    }
  } else {
    IWS = N;
  }

  // Solve the equation inv(A)*L = inv(U) for inv(A).

  if (NB < NBMIN || NB >= N) {
    // Use unblocked code.

    for (J = N; J >= 1; J--) {
      // Copy current column of L to WORK and replace with zeros.

      for (I = J + 1; I <= N; I++) {
        WORK[I] = A[I][J];
        A[I][J] = Complex.zero;
      }

      // Compute current column of inv(A).

      if (J < N) {
        zgemv('No transpose', N, N - J, -Complex.one, A(1, J + 1), LDA,
            WORK(J + 1), 1, Complex.one, A(1, J).asArray(), 1);
      }
    }
  } else {
    // Use blocked code.

    NN = ((N - 1) ~/ NB) * NB + 1;
    for (J = NN; J >= 1; J -= NB) {
      JB = min(NB, N - J + 1);

      // Copy current block column of L to WORK and replace with
      // zeros.

      for (JJ = J; JJ <= J + JB - 1; JJ++) {
        for (I = JJ + 1; I <= N; I++) {
          WORK[I + (JJ - J) * LDWORK] = A[I][JJ];
          A[I][JJ] = Complex.zero;
        }
      }

      // Compute current block column of inv(A).

      if (J + JB <= N) {
        zgemm(
            'No transpose',
            'No transpose',
            N,
            JB,
            N - J - JB + 1,
            -Complex.one,
            A(1, J + JB),
            LDA,
            WORK(J + JB).asMatrix(LDWORK),
            LDWORK,
            Complex.one,
            A(1, J),
            LDA);
      }
      ztrsm('Right', 'Lower', 'No transpose', 'Unit', N, JB, Complex.one,
          WORK(J).asMatrix(LDWORK), LDWORK, A(1, J), LDA);
    }
  }

  // Apply column interchanges.

  for (J = N - 1; J >= 1; J--) {
    JP = IPIV[J];
    if (JP != J) zswap(N, A(1, J).asArray(), 1, A(1, JP).asArray(), 1);
  }

  WORK[1] = IWS.toComplex();
}
