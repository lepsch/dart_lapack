// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlaset.dart';

import 'zlavhe.dart';

void zhet01(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> AFAC_,
  final int LDAFAC,
  final Array<int> IPIV_,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
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
  final ANORM = zlanhe('1', UPLO, N, A, LDA, RWORK);

  // Check the imaginary parts of the diagonal elements and return with
  // an error code if any are nonzero.

  for (var J = 1; J <= N; J++) {
    if (AFAC[J][J].imaginary != ZERO) {
      RESID.value = ONE / EPS;
      return;
    }
  }

  // Initialize C to the identity matrix.

  zlaset('Full', N, N, Complex.zero, Complex.one, C, LDC);

  // Call ZLAVHE to form the product D * U' (or D * L' ).

  zlavhe(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO);

  // Call ZLAVHE again to multiply by U (or L ).

  zlavhe(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO);

  // Compute the difference  C - A .

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J - 1; I++) {
        C[I][J] -= A[I][J];
      }
      C[J][J] -= A[J][J].real.toComplex();
    }
  } else {
    for (var J = 1; J <= N; J++) {
      C[J][J] -= A[J][J].real.toComplex();
      for (var I = J + 1; I <= N; I++) {
        C[I][J] -= A[I][J];
      }
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
