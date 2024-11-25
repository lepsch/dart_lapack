// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zhpr.dart';
import 'package:dart_lapack/src/blas/ztpmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/ztptri.dart';

void zpptri(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  const ONE = 1.0;
  bool UPPER;
  int J, JC, JJ, JJN;
  double AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('ZPPTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Invert the triangular Cholesky factor U or L.

  ztptri(UPLO, 'Non-unit', N, AP, INFO);
  if (INFO.value > 0) return;
  if (UPPER) {
    // Compute the product inv(U) * inv(U)**H.

    JJ = 0;
    for (J = 1; J <= N; J++) {
      JC = JJ + 1;
      JJ += J;
      if (J > 1) zhpr('Upper', J - 1, ONE, AP(JC), 1, AP);
      AJJ = AP[JJ].real;
      zdscal(J, AJJ, AP(JC), 1);
    }
  } else {
    // Compute the product inv(L)**H * inv(L).

    JJ = 1;
    for (J = 1; J <= N; J++) {
      JJN = JJ + N - J + 1;
      AP[JJ] = zdotc(N - J + 1, AP(JJ), 1, AP(JJ), 1).real.toComplex();
      if (J < N) {
        ztpmv('Lower', 'Conjugate transpose', 'Non-unit', N - J, AP(JJN),
            AP(JJ + 1), 1);
      }
      JJ = JJN;
    }
  }
}
