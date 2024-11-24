// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/lapack.dart';

void zlahilb(
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> WORK_,
  final Box<int> INFO,
  final String PATH,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA, lastd: N);
  final X = X_.having(ld: LDX, lastd: NRHS);
  final B = B_.having(ld: LDB, lastd: NRHS);
  final WORK = WORK_.having(length: N);
  const NMAX_EXACT = 6, NMAX_APPROX = 11, SIZE_D = 8;

  // d's are generated from random permutation of those eight elements.
  // Complex d1(8), d2(8), invd1(8), invd2(8);
  final D1 = Array.fromList([
    (-1, 0),
    (0, 1),
    (-1, -1),
    (0, -1),
    (1, 0),
    (-1, 1),
    (1, 1),
    (1, -1)
  ].toComplexList());
  final D2 = Array.fromList([
    (-1, 0),
    (0, -1),
    (-1, 1),
    (0, 1),
    (1, 0),
    (-1, -1),
    (1, -1),
    (1, 1)
  ].toComplexList());
  final INVD1 = Array.fromList([
    (-1, 0),
    (0, -1),
    (-.5, .5),
    (0, 1),
    (1, 0),
    (-.5, -.5),
    (.5, -.5),
    (.5, .5)
  ].toComplexList());
  final INVD2 = Array.fromList([
    (-1, 0),
    (0, 1),
    (-.5, -.5),
    (0, -1),
    (1, 0),
    (-.5, .5),
    (.5, .5),
    (.5, -.5)
  ].toComplexList());
  final C2 = PATH.substring(1, 3);

  // Test the input arguments

  INFO.value = 0;
  if (N < 0 || N > NMAX_APPROX) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDA < N) {
    INFO.value = -4;
  } else if (LDX < N) {
    INFO.value = -6;
  } else if (LDB < N) {
    INFO.value = -8;
  }
  if (INFO.value < 0) {
    xerbla('ZLAHILB', -INFO.value);
    return;
  }
  if (N > NMAX_EXACT) {
    INFO.value = 1;
  }

  // Compute M = the LCM of the integers [1, 2*N-1].  The largest
  // reasonable N is small enough that integers suffice (up to N = 11).
  var M = 1;
  for (var I = 2; I <= (2 * N - 1); I++) {
    var TM = M;
    var TI = I;
    var R = (TM % TI);
    while (R != 0) {
      TM = TI;
      TI = R;
      R = (TM % TI);
    }
    M = (M ~/ TI) * I;
  }

  // Generate the scaled Hilbert matrix in A
  // If we are testing SY routines,
  //    take D1_i = D2_i, else, D1_i = D2_i*
  if (lsamen(2, C2, 'SY')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= N; I++) {
        A[I][J] = D1[(J % SIZE_D) + 1] *
            (M / (I + J - 1)).toComplex() *
            D1[(I % SIZE_D) + 1];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= N; I++) {
        A[I][J] = D1[(J % SIZE_D) + 1] *
            (M / (I + J - 1)).toComplex() *
            D2[(I % SIZE_D) + 1];
      }
    }
  }

  // Generate matrix B as simply the first NRHS columns of M * the
  // identity.
  final TMP = M.toComplex();
  zlaset('Full', N, NRHS, Complex.zero, TMP, B, LDB);

  // Generate the true solutions in X.  Because B = the first NRHS
  // columns of M*I, the true solutions are just the first NRHS columns
  // of the inverse Hilbert matrix.
  WORK[1] = N.toDouble();
  for (var J = 2; J <= N; J++) {
    WORK[J] = (((WORK[J - 1] / (J - 1)) * (J - 1 - N)) / (J - 1)) * (N + J - 1);
  }

  // If we are testing SY routines,
  //       take D1_i = D2_i, else, D1_i = D2_i*
  if (lsamen(2, C2, 'SY')) {
    for (var J = 1; J <= NRHS; J++) {
      for (var I = 1; I <= N; I++) {
        X[I][J] = INVD1[(J % SIZE_D) + 1] *
            ((WORK[I] * WORK[J]) / (I + J - 1)).toComplex() *
            INVD1[(I % SIZE_D) + 1];
      }
    }
  } else {
    for (var J = 1; J <= NRHS; J++) {
      for (var I = 1; I <= N; I++) {
        X[I][J] = INVD2[(J % SIZE_D) + 1] *
            ((WORK[I] * WORK[J]) / (I + J - 1)).toComplex() *
            INVD1[(I % SIZE_D) + 1];
      }
    }
  }
}
