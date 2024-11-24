// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zungqr.dart';

void zunghr(
  final int N,
  final int ILO,
  final int IHI,
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
  int I, J, LWKOPT = 0, NB, NH;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  NH = IHI - ILO;
  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -1;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -2;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LWORK < max(1, NH) && !LQUERY) {
    INFO.value = -8;
  }

  if (INFO.value == 0) {
    NB = ilaenv(1, 'ZUNGQR', ' ', NH, NH, NH, -1);
    LWKOPT = max(1, NH) * NB;
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZUNGHR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    WORK[1] = Complex.one;
    return;
  }

  // Shift the vectors which define the elementary reflectors one
  // column to the right, and set the first ilo and the last n-ihi
  // rows and columns to those of the unit matrix

  for (J = IHI; J >= ILO + 1; J--) {
    for (I = 1; I <= J - 1; I++) {
      A[I][J] = Complex.zero;
    }
    for (I = J + 1; I <= IHI; I++) {
      A[I][J] = A[I][J - 1];
    }
    for (I = IHI + 1; I <= N; I++) {
      A[I][J] = Complex.zero;
    }
  }
  for (J = 1; J <= ILO; J++) {
    for (I = 1; I <= N; I++) {
      A[I][J] = Complex.zero;
    }
    A[J][J] = Complex.one;
  }
  for (J = IHI + 1; J <= N; J++) {
    for (I = 1; I <= N; I++) {
      A[I][J] = Complex.zero;
    }
    A[J][J] = Complex.one;
  }

  if (NH > 0) {
    // Generate Q(ilo+1:ihi,ilo+1:ihi)

    zungqr(NH, NH, NH, A(ILO + 1, ILO + 1), LDA, TAU(ILO), WORK, LWORK, IINFO);
  }
  WORK[1] = LWKOPT.toComplex();
}
