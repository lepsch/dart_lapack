// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zlakf2(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> D_,
  final Matrix<Complex> E_,
  final Matrix<Complex> Z_,
  final int LDZ,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final D = D_.having(ld: LDA);
  final E = E_.having(ld: LDA);
  final Z = Z_.having(ld: LDA);
  int I, IK, J, JK, L, MN, MN2;
  MN = M * N;
  MN2 = 2 * MN;
  zlaset('Full', MN2, MN2, Complex.zero, Complex.zero, Z, LDZ);

  IK = 1;
  for (L = 1; L <= N; L++) {
    // form kron(In, A)

    for (I = 1; I <= M; I++) {
      for (J = 1; J <= M; J++) {
        Z[IK + I - 1][IK + J - 1] = A[I][J];
      }
    }

    // form kron(In, D)

    for (I = 1; I <= M; I++) {
      for (J = 1; J <= M; J++) {
        Z[IK + MN + I - 1][IK + J - 1] = D[I][J];
      }
    }

    IK += M;
  }

  IK = 1;
  for (L = 1; L <= N; L++) {
    JK = MN + 1;

    for (J = 1; J <= N; J++) {
      // form -kron(B', Im)

      for (I = 1; I <= M; I++) {
        Z[IK + I - 1][JK + I - 1] = -B[J][L];
      }

      // form -kron(E', Im)

      for (I = 1; I <= M; I++) {
        Z[IK + MN + I - 1][JK + I - 1] = -E[J][L];
      }

      JK += M;
    }

    IK += M;
  }
}
