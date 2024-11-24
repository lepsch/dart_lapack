// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/drot.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlartg.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlaqz2(
  final bool ILQ,
  final bool ILZ,
  final int K,
  final int ISTARTM,
  final int ISTOPM,
  final int IHI,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final int NQ,
  final int QSTART,
  final Matrix<double> Q_,
  final int LDQ,
  final int NZ,
  final int ZSTART,
  final Matrix<double> Z_,
  final int LDZ,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  const ZERO = 0.0;
  final C1 = Box(0.0),
      S1 = Box(0.0),
      C2 = Box(0.0),
      S2 = Box(0.0),
      TEMP = Box(0.0);
  final H = Matrix<double>(2, 3);

  if (K + 2 == IHI) {
    // Shift is located on the edge of the matrix, remove it
    //  H = B( IHI-1:IHI, IHI-2:IHI );
    for (var i = IHI - 1; i <= IHI; i++) {
      for (var j = IHI - 2; j <= IHI; j++) {
        H[i - (IHI - 1) + 1][j - (IHI - 2) + 1] = B[i][j];
      }
    }
    // Make H upper triangular
    dlartg(H[1][1], H[2][1], C1, S1, TEMP);
    H[2][1] = ZERO;
    H[1][1] = TEMP.value;
    drot(2, H(1, 2).asArray(), 2, H(2, 2).asArray(), 2, C1.value, S1.value);

    dlartg(H[2][3], H[2][2], C1, S1, TEMP);
    drot(1, H(1, 3).asArray(), 1, H(1, 2).asArray(), 1, C1.value, S1.value);
    dlartg(H[1][2], H[1][1], C2, S2, TEMP);

    drot(IHI - ISTARTM + 1, B(ISTARTM, IHI).asArray(), 1,
        B(ISTARTM, IHI - 1).asArray(), 1, C1.value, S1.value);
    drot(IHI - ISTARTM + 1, B(ISTARTM, IHI - 1).asArray(), 1,
        B(ISTARTM, IHI - 2).asArray(), 1, C2.value, S2.value);
    B[IHI - 1][IHI - 2] = ZERO;
    B[IHI][IHI - 2] = ZERO;
    drot(IHI - ISTARTM + 1, A(ISTARTM, IHI).asArray(), 1,
        A(ISTARTM, IHI - 1).asArray(), 1, C1.value, S1.value);
    drot(IHI - ISTARTM + 1, A(ISTARTM, IHI - 1).asArray(), 1,
        A(ISTARTM, IHI - 2).asArray(), 1, C2.value, S2.value);
    if (ILZ) {
      drot(NZ, Z(1, IHI - ZSTART + 1).asArray(), 1,
          Z(1, IHI - 1 - ZSTART + 1).asArray(), 1, C1.value, S1.value);
      drot(NZ, Z(1, IHI - 1 - ZSTART + 1).asArray(), 1,
          Z(1, IHI - 2 - ZSTART + 1).asArray(), 1, C2.value, S2.value);
    }

    dlartg(A[IHI - 1][IHI - 2], A[IHI][IHI - 2], C1, S1, TEMP);
    A[IHI - 1][IHI - 2] = TEMP.value;
    A[IHI][IHI - 2] = ZERO;
    drot(ISTOPM - IHI + 2, A(IHI - 1, IHI - 1).asArray(), LDA,
        A(IHI, IHI - 1).asArray(), LDA, C1.value, S1.value);
    drot(ISTOPM - IHI + 2, B(IHI - 1, IHI - 1).asArray(), LDB,
        B(IHI, IHI - 1).asArray(), LDB, C1.value, S1.value);
    if (ILQ) {
      drot(NQ, Q(1, IHI - 1 - QSTART + 1).asArray(), 1,
          Q(1, IHI - QSTART + 1).asArray(), 1, C1.value, S1.value);
    }

    dlartg(B[IHI][IHI], B[IHI][IHI - 1], C1, S1, TEMP);
    B[IHI][IHI] = TEMP.value;
    B[IHI][IHI - 1] = ZERO;
    drot(IHI - ISTARTM, B(ISTARTM, IHI).asArray(), 1,
        B(ISTARTM, IHI - 1).asArray(), 1, C1.value, S1.value);
    drot(IHI - ISTARTM + 1, A(ISTARTM, IHI).asArray(), 1,
        A(ISTARTM, IHI - 1).asArray(), 1, C1.value, S1.value);
    if (ILZ) {
      drot(NZ, Z(1, IHI - ZSTART + 1).asArray(), 1,
          Z(1, IHI - 1 - ZSTART + 1).asArray(), 1, C1.value, S1.value);
    }
  } else {
    // Normal operation, move bulge down
    //  H = B( K+1:K+2, K:K+2 );
    for (var i = K + 1; i <= K + 2; i++) {
      for (var j = K; j <= K + 2; j++) {
        H[i - (K + 1) + 1][j - K + 1] = B[i][j];
      }
    }

    // Make H upper triangular

    dlartg(H[1][1], H[2][1], C1, S1, TEMP);
    H[2][1] = ZERO;
    H[1][1] = TEMP.value;
    drot(2, H(1, 2).asArray(), 2, H(2, 2).asArray(), 2, C1.value, S1.value);

    // Calculate Z1 and Z2

    dlartg(H[2][3], H[2][2], C1, S1, TEMP);
    drot(1, H(1, 3).asArray(), 1, H(1, 2).asArray(), 1, C1.value, S1.value);
    dlartg(H[1][2], H[1][1], C2, S2, TEMP);

    // Apply transformations from the right

    drot(K + 3 - ISTARTM + 1, A(ISTARTM, K + 2).asArray(), 1,
        A(ISTARTM, K + 1).asArray(), 1, C1.value, S1.value);
    drot(K + 3 - ISTARTM + 1, A(ISTARTM, K + 1).asArray(), 1,
        A(ISTARTM, K).asArray(), 1, C2.value, S2.value);
    drot(K + 2 - ISTARTM + 1, B(ISTARTM, K + 2).asArray(), 1,
        B(ISTARTM, K + 1).asArray(), 1, C1.value, S1.value);
    drot(K + 2 - ISTARTM + 1, B(ISTARTM, K + 1).asArray(), 1,
        B(ISTARTM, K).asArray(), 1, C2.value, S2.value);
    if (ILZ) {
      drot(NZ, Z(1, K + 2 - ZSTART + 1).asArray(), 1,
          Z(1, K + 1 - ZSTART + 1).asArray(), 1, C1.value, S1.value);
      drot(NZ, Z(1, K + 1 - ZSTART + 1).asArray(), 1,
          Z(1, K - ZSTART + 1).asArray(), 1, C2.value, S2.value);
    }
    B[K + 1][K] = ZERO;
    B[K + 2][K] = ZERO;

    // Calculate Q1 and Q2

    dlartg(A[K + 2][K], A[K + 3][K], C1, S1, TEMP);
    A[K + 2][K] = TEMP.value;
    A[K + 3][K] = ZERO;
    dlartg(A[K + 1][K], A[K + 2][K], C2, S2, TEMP);
    A[K + 1][K] = TEMP.value;
    A[K + 2][K] = ZERO;

    // Apply transformations from the left

    drot(ISTOPM - K, A(K + 2, K + 1).asArray(), LDA, A(K + 3, K + 1).asArray(),
        LDA, C1.value, S1.value);
    drot(ISTOPM - K, A(K + 1, K + 1).asArray(), LDA, A(K + 2, K + 1).asArray(),
        LDA, C2.value, S2.value);

    drot(ISTOPM - K, B(K + 2, K + 1).asArray(), LDB, B(K + 3, K + 1).asArray(),
        LDB, C1.value, S1.value);
    drot(ISTOPM - K, B(K + 1, K + 1).asArray(), LDB, B(K + 2, K + 1).asArray(),
        LDB, C2.value, S2.value);
    if (ILQ) {
      drot(NQ, Q(1, K + 2 - QSTART + 1).asArray(), 1,
          Q(1, K + 3 - QSTART + 1).asArray(), 1, C1.value, S1.value);
      drot(NQ, Q(1, K + 1 - QSTART + 1).asArray(), 1,
          Q(1, K + 2 - QSTART + 1).asArray(), 1, C2.value, S2.value);
    }
  }
}
