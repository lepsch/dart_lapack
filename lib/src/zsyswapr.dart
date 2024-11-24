// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zsyswapr(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final int I1,
  final int I2,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  bool UPPER;
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
    TMP = A[I1][I1];
    A[I1][I1] = A[I2][I2];
    A[I2][I2] = TMP;

    zswap(
        I2 - I1 - 1, A(I1, I1 + 1).asArray(), LDA, A(I1 + 1, I2).asArray(), 1);

    // third swap
    // - swap row I1 and I2 from I2+1 to N
    if (I2 < N) {
      zswap(N - I2, A(I1, I2 + 1).asArray(), LDA, A(I2, I2 + 1).asArray(), LDA);
    }
  } else {
    // LOWER
    // first swap
    //  - swap row I1 and I2 from I1 to I1-1
    zswap(I1 - 1, A(I1, 1).asArray(), LDA, A(I2, 1).asArray(), LDA);

    // second swap :
    //  - swap A[I1][I1] and A[I2][I2]
    //  - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
    TMP = A[I1][I1];
    A[I1][I1] = A[I2][I2];
    A[I2][I2] = TMP;

    zswap(
        I2 - I1 - 1, A(I1 + 1, I1).asArray(), 1, A(I2, I1 + 1).asArray(), LDA);

    // third swap
    //  - swap col I1 and I2 from I2+1 to N
    if (I2 < N) {
      zswap(N - I2, A(I2 + 1, I1).asArray(), 1, A(I2 + 1, I2).asArray(), 1);
    }
  }
}
