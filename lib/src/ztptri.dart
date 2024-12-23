// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/blas/ztpmv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void ztptri(
  final String UPLO,
  final String DIAG,
  final int N,
  final Array<Complex> AP_,
  final Box<int> INFO,
) {
  final AP = AP_.having();
  bool NOUNIT, UPPER;
  int J, JC, JCLAST = 0, JJ;
  Complex AJJ;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  NOUNIT = lsame(DIAG, 'N');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('ZTPTRI', -INFO.value);
    return;
  }

  // Check for singularity if non-unit.

  if (NOUNIT) {
    if (UPPER) {
      JJ = 0;
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        JJ += INFO.value;
        if (AP[JJ] == Complex.zero) return;
      }
    } else {
      JJ = 1;
      for (INFO.value = 1; INFO.value <= N; INFO.value++) {
        if (AP[JJ] == Complex.zero) return;
        JJ += N - INFO.value + 1;
      }
    }
    INFO.value = 0;
  }

  if (UPPER) {
    // Compute inverse of upper triangular matrix.

    JC = 1;
    for (J = 1; J <= N; J++) {
      if (NOUNIT) {
        AP[JC + J - 1] = Complex.one / AP[JC + J - 1];
        AJJ = -AP[JC + J - 1];
      } else {
        AJJ = -Complex.one;
      }

      // Compute elements 1:j-1 of j-th column.

      ztpmv('Upper', 'No transpose', DIAG, J - 1, AP, AP(JC), 1);
      zscal(J - 1, AJJ, AP(JC), 1);
      JC += J;
    }
  } else {
    // Compute inverse of lower triangular matrix.

    JC = N * (N + 1) ~/ 2;
    for (J = N; J >= 1; J--) {
      if (NOUNIT) {
        AP[JC] = Complex.one / AP[JC];
        AJJ = -AP[JC];
      } else {
        AJJ = -Complex.one;
      }
      if (J < N) {
        // Compute elements j+1:n of j-th column.

        ztpmv('Lower', 'No transpose', DIAG, N - J, AP(JCLAST), AP(JC + 1), 1);
        zscal(N - J, AJJ, AP(JC + 1), 1);
      }
      JCLAST = JC;
      JC -= N - J + 2;
    }
  }
}
