// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zdrscl.dart';

void zrscl(
  final int N,
  final Complex A,
  final Array<Complex> X_,
  final int INCX,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  const ZERO = 0.0, ONE = 1.0;
  double SAFMAX, SAFMIN, OV, AR, AI, ABSR, ABSI, UR, UI;

  // Quick return if possible

  if (N <= 0) return;

  // Get machine parameters

  SAFMIN = dlamch('S');
  SAFMAX = ONE / SAFMIN;
  OV = dlamch('O');

  // Initialize constants related to A.

  AR = A.real;
  AI = A.imaginary;
  ABSR = AR.abs();
  ABSI = AI.abs();

  if (AI == ZERO) {
    // If alpha is real, then we can use csrscl
    zdrscl(N, AR, X, INCX);
  } else if (AR == ZERO) {
    // If alpha has a zero real part, then we follow the same rules as if
    // alpha were real.
    if (ABSI > SAFMAX) {
      zdscal(N, SAFMIN, X, INCX);
      zscal(N, Complex(ZERO, -SAFMAX / AI), X, INCX);
    } else if (ABSI < SAFMIN) {
      zscal(N, Complex(ZERO, -SAFMIN / AI), X, INCX);
      zdscal(N, SAFMAX, X, INCX);
    } else {
      zscal(N, Complex(ZERO, -ONE / AI), X, INCX);
    }
  } else {
    // The following numbers can be computed.
    // They are the inverse of the real and imaginary parts of 1/alpha.
    // Note that a and b are always different from zero.
    // NaNs are only possible if either:
    // 1. alphaR or alphaI is NaN.
    // 2. alphaR and alphaI are both infinite, in which case it makes sense
    // to propagate a NaN.
    UR = AR + AI * (AI / AR);
    UI = AI + AR * (AR / AI);

    if ((UR.abs() < SAFMIN) || (UI.abs() < SAFMIN)) {
      // This means that both alphaR and alphaI are very small.
      zscal(N, Complex(SAFMIN / UR, -SAFMIN / UI), X, INCX);
      zdscal(N, SAFMAX, X, INCX);
    } else if ((UR.abs() > SAFMAX) || (UI.abs() > SAFMAX)) {
      if ((ABSR > OV) || (ABSI > OV)) {
        // This means that a and b are both Inf. No need for scaling.
        zscal(N, Complex(ONE / UR, -ONE / UI), X, INCX);
      } else {
        zdscal(N, SAFMIN, X, INCX);
        if ((UR.abs() > OV) || (UI.abs() > OV)) {
          // Infs were generated. We do proper scaling to avoid them.
          if (ABSR >= ABSI) {
            // ABS( UR ) <= ABS( UI )
            UR = (SAFMIN * AR) + SAFMIN * (AI * (AI / AR));
            UI = (SAFMIN * AI) + AR * ((SAFMIN * AR) / AI);
          } else {
            // ABS( UR ) > ABS( UI )
            UR = (SAFMIN * AR) + AI * ((SAFMIN * AI) / AR);
            UI = (SAFMIN * AI) + SAFMIN * (AR * (AR / AI));
          }
          zscal(N, Complex(ONE / UR, -ONE / UI), X, INCX);
        } else {
          zscal(N, Complex(SAFMAX / UR, -SAFMAX / UI), X, INCX);
        }
      }
    } else {
      zscal(N, Complex(ONE / UR, -ONE / UI), X, INCX);
    }
  }
}
