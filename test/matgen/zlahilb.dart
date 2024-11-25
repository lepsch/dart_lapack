// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlaset.dart';

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
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having(length: N);
  int TM, TI, R;
  int M;
  int I, J;
  Complex TMP;
  String C2;
  const NMAX_EXACT = 6, NMAX_APPROX = 11, SIZE_D = 8;

  // D's are generated from random permutation of those eight elements.
  // final D1=Array<Complex>(8), D2=Array<Complex>(8), INVD1=Array<Complex>(8), INVD2=Array<Complex>(8);
  final D1 = Array.fromList([
    -Complex.one, Complex(0.0, 1.0), Complex(-1.0, -1.0), Complex(0.0, -1.0), //
    Complex.one, Complex(-1.0, 1.0), Complex(1.0, 1.0), Complex(1.0, -1.0)
  ]);
  final D2 = Array.fromList([
    -Complex.one, Complex(0.0, -1.0), Complex(-1.0, 1.0), Complex(0.0, 1.0), //
    Complex.one, Complex(-1.0, -1.0), Complex(1.0, -1.0),
    Complex(1.0, 1.0)
  ]);
  final INVD1 = Array.fromList([
    -Complex.one, Complex(0.0, -1.0), Complex(-0.5, 0.5), Complex(0.0, 1.0), //
    Complex.one, Complex(-0.5, -0.5), Complex(0.5, -0.5),
    Complex(0.5, 0.5)
  ]);
  final INVD2 = Array.fromList([
    -Complex.one, Complex(0.0, 1.0), Complex(-0.5, -0.5), Complex(0.0, -1.0), //
    Complex.one, Complex(-0.5, 0.5), Complex(0.5, 0.5), Complex(0.5, -0.5)
  ]);

  C2 = PATH.substring(1, 3);

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
  M = 1;
  for (I = 2; I <= (2 * N - 1); I++) {
    TM = M;
    TI = I;
    R = (TM % TI);
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
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= N; I++) {
        A[I][J] = D1[(J % SIZE_D) + 1] *
            (M / (I + J - 1)).toComplex() *
            D1[(I % SIZE_D) + 1];
      }
    }
  } else {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= N; I++) {
        A[I][J] = D1[(J % SIZE_D) + 1] *
            (M / (I + J - 1)).toComplex() *
            D2[(I % SIZE_D) + 1];
      }
    }
  }

  // Generate matrix B as simply the first NRHS columns of M * the
  // identity.
  TMP = M.toComplex();
  zlaset('Full', N, NRHS, Complex.zero, TMP, B, LDB);

  // Generate the true solutions in X.  Because B = the first NRHS
  // columns of M*I, the true solutions are just the first NRHS columns
  // of the inverse Hilbert matrix.
  WORK[1] = N.toDouble();
  for (J = 2; J <= N; J++) {
    WORK[J] = (((WORK[J - 1] / (J - 1)) * (J - 1 - N)) / (J - 1)) * (N + J - 1);
  }

  // If we are testing SY routines,
  //       take D1_i = D2_i, else, D1_i = D2_i*
  if (lsamen(2, C2, 'SY')) {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        X[I][J] = INVD1[(J % SIZE_D) + 1] *
            ((WORK[I] * WORK[J]) / (I + J - 1)).toComplex() *
            INVD1[(I % SIZE_D) + 1];
      }
    }
  } else {
    for (J = 1; J <= NRHS; J++) {
      for (I = 1; I <= N; I++) {
        X[I][J] = INVD2[(J % SIZE_D) + 1] *
            ((WORK[I] * WORK[J]) / (I + J - 1)).toComplex() *
            INVD1[(I % SIZE_D) + 1];
      }
    }
  }
}
