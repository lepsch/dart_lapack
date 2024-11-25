// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dgemv.dart';
import 'package:dart_lapack/src/blas/dger.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/iladlc.dart';
import 'package:dart_lapack/src/iladlr.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarf(
  final String SIDE,
  final int M,
  final int N,
  final Array<double> V_,
  final int INCV,
  final double TAU,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK,
) {
  final V = V_.having();
  final C = C_.having(ld: LDC);
  const ONE = 1.0, ZERO = 0.0;
  bool APPLYLEFT;
  int I, LASTV, LASTC;

  APPLYLEFT = lsame(SIDE, 'L');
  LASTV = 0;
  LASTC = 0;
  if (TAU != ZERO) {
    // Set up variables for scanning V.  LASTV begins pointing to the end
    // of V.
    if (APPLYLEFT) {
      LASTV = M;
    } else {
      LASTV = N;
    }
    if (INCV > 0) {
      I = 1 + (LASTV - 1) * INCV;
    } else {
      I = 1;
    }
    // Look for the last non-zero row in V.
    while (LASTV > 0 && V[I] == ZERO) {
      LASTV--;
      I -= INCV;
    }
    if (APPLYLEFT) {
      // Scan for the last non-zero column in C(1:lastv,:).
      LASTC = iladlc(LASTV, N, C, LDC);
    } else {
      // Scan for the last non-zero row in C(:,1:lastv).
      LASTC = iladlr(M, LASTV, C, LDC);
    }
  }
  // Note that lastc == 0 renders the BLAS operations null; no special
  // case is needed at this level.
  if (APPLYLEFT) {
    // Form  H * C

    if (LASTV > 0) {
      // w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)

      dgemv('Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, ZERO, WORK, 1);

      // C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T

      dger(LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC);
    }
  } else {
    // Form  C * H

    if (LASTV > 0) {
      // w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)

      dgemv('No transpose', LASTC, LASTV, ONE, C, LDC, V, INCV, ZERO, WORK, 1);

      // C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T

      dger(LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC);
    }
  }
}
