// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlakf2(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> D_,
  final Matrix<double> E_,
  final Matrix<double> Z_,
  final int LDZ,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final D = D_.having(ld: LDA);
  final E = E_.having(ld: LDA);
  final Z = Z_.having(ld: LDZ);
  const ZERO = 0.0;
  int I, IK, J, JK, L, MN, MN2;

  // Initialize Z

  MN = M * N;
  MN2 = 2 * MN;
  dlaset('Full', MN2, MN2, ZERO, ZERO, Z, LDZ);

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
