// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlaswp(
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final int K1,
  final int K2,
  final Array<int> IPIV_,
  final int INCX,
) {
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  int I, I1, I2, INC, IP, IX, IX0, J, K, N32;
  Complex TEMP;

  // Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
  // K1 through K2.

  if (INCX > 0) {
    IX0 = K1;
    I1 = K1;
    I2 = K2;
    INC = 1;
  } else if (INCX < 0) {
    IX0 = K1 + (K1 - K2) * INCX;
    I1 = K2;
    I2 = K1;
    INC = -1;
  } else {
    return;
  }

  N32 = (N ~/ 32) * 32;
  if (N32 != 0) {
    for (J = 1; J <= N32; J += 32) {
      IX = IX0;
      for (I = I1; INC < 0 ? I >= I2 : I <= I2; I += INC) {
        IP = IPIV[IX];
        if (IP != I) {
          for (K = J; K <= J + 31; K++) {
            TEMP = A[I][K];
            A[I][K] = A[IP][K];
            A[IP][K] = TEMP;
          }
        }
        IX += INCX;
      }
    }
  }
  if (N32 != N) {
    N32++;
    IX = IX0;
    for (I = I1; INC < 0 ? I >= I2 : I <= I2; I += INC) {
      IP = IPIV[IX];
      if (IP != I) {
        for (K = N32; K <= N; K++) {
          TEMP = A[I][K];
          A[I][K] = A[IP][K];
          A[IP][K] = TEMP;
        }
      }
      IX += INCX;
    }
  }
}
