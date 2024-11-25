// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zheswapr(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final int I1,
  final int I2,
) {
  final A = A_.having(ld: LDA);
  bool UPPER;
  int I;
  Complex TMP;

  UPPER = lsame(UPLO, 'U');
  if (UPPER) {
    // UPPER
    // first swap
    //  - swap column I1 and I2 from I1 to I1-1
    zswap(I1 - 1, A(1, I1).asArray(), 1, A(1, I2).asArray(), 1);

    // second swap :
    // - swap A[I1][I1] and A[I2][I2]
    // - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
    // - swap A[I2][I1] and A[I1][I2]

    TMP = A[I1][I1];
    A[I1][I1] = A[I2][I2];
    A[I2][I2] = TMP;

    for (I = 1; I <= I2 - I1 - 1; I++) {
      TMP = A[I1][I1 + I];
      A[I1][I1 + I] = A[I1 + I][I2].conjugate();
      A[I1 + I][I2] = TMP.conjugate();
    }

    A[I1][I2] = A[I1][I2].conjugate();

    // third swap
    // - swap row I1 and I2 from I2+1 to N
    for (I = I2 + 1; I <= N; I++) {
      TMP = A[I1][I];
      A[I1][I] = A[I2][I];
      A[I2][I] = TMP;
    }
  } else {
    // LOWER
    // first swap
    //  - swap row I1 and I2 from 1 to I1-1
    zswap(I1 - 1, A(I1, 1).asArray(), LDA, A(I2, 1).asArray(), LDA);

    // second swap :
    //  - swap A[I1][I1] and A[I2][I2]
    //  - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
    //  - swap A[I2][I1] and A[I1][I2]

    TMP = A[I1][I1];
    A[I1][I1] = A[I2][I2];
    A[I2][I2] = TMP;

    for (I = 1; I <= I2 - I1 - 1; I++) {
      TMP = A[I1 + I][I1];
      A[I1 + I][I1] = A[I2][I1 + I].conjugate();
      A[I2][I1 + I] = TMP.conjugate();
    }

    A[I2][I1] = A[I2][I1].conjugate();

    // third swap
    //  - swap col I1 and I2 from I2+1 to N
    for (I = I2 + 1; I <= N; I++) {
      TMP = A[I][I1];
      A[I][I1] = A[I][I2];
      A[I][I2] = TMP;
    }
  }
}
