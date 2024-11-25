// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlanhp.dart';
import 'package:dart_lapack/src/zlaset.dart';

import 'zlavhp.dart';

void zhpt01(
  final String UPLO,
  final int N,
  final Array<Complex> A_,
  final Array<Complex> AFAC_,
  final Array<int> IPIV_,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having();
  final AFAC = AFAC_.having();
  final IPIV = IPIV_.having();
  final C = C_.having(ld: LDC);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  final INFO = Box(0);

  // Quick exit if N = 0.

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Determine EPS and the norm of A.

  final EPS = dlamch('Epsilon');
  final ANORM = zlanhp('1', UPLO, N, A, RWORK);

  // Check the imaginary parts of the diagonal elements and return with
  // an error code if any are nonzero.

  var JC = 1;
  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      if (AFAC[JC].imaginary != ZERO) {
        RESID.value = ONE / EPS;
        return;
      }
      JC += J + 1;
    }
  } else {
    for (var J = 1; J <= N; J++) {
      if (AFAC[JC].imaginary != ZERO) {
        RESID.value = ONE / EPS;
        return;
      }
      JC += N - J + 1;
    }
  }

  // Initialize C to the identity matrix.

  zlaset('Full', N, N, Complex.zero, Complex.one, C, LDC);

  // Call ZLAVHP to form the product D * U' (or D * L' ).

  zlavhp(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO);

  // Call ZLAVHP again to multiply by U ( or L ).

  zlavhp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO);

  // Compute the difference  C - A .

  if (lsame(UPLO, 'U')) {
    var JC = 0;
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J - 1; I++) {
        C[I][J] -= A[JC + I];
      }
      C[J][J] -= A[JC + J].real.toComplex();
      JC += J;
    }
  } else {
    var JC = 1;
    for (var J = 1; J <= N; J++) {
      C[J][J] -= A[JC].real.toComplex();
      for (var I = J + 1; I <= N; I++) {
        C[I][J] -= A[JC + I - J];
      }
      JC += N - J + 1;
    }
  }

  // Compute norm( C - A ) / ( N * norm(A) * EPS )

  RESID.value = zlanhe('1', UPLO, N, C, LDC, RWORK);

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }
}
