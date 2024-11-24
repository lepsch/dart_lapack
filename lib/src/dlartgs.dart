// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlartgp.dart';
import 'package:dart_lapack/src/install/dlamch.dart';

void dlartgs(
  final double X,
  final double Y,
  final double SIGMA,
  final Box<double> CS,
  final Box<double> SN,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NEGONE = -1.0, ONE = 1.0, ZERO = 0.0;
  double S, THRESH, W, Z;
  final R = Box(0.0);

  THRESH = dlamch('E');

  // Compute the first column of B**T*B - SIGMA^2*I, up to a scale
  // factor.

  if ((SIGMA == ZERO && X.abs() < THRESH) || (X.abs() == SIGMA && Y == ZERO)) {
    Z = ZERO;
    W = ZERO;
  } else if (SIGMA == ZERO) {
    if (X >= ZERO) {
      Z = X;
      W = Y;
    } else {
      Z = -X;
      W = -Y;
    }
  } else if (X.abs() < THRESH) {
    Z = -SIGMA * SIGMA;
    W = ZERO;
  } else {
    if (X >= ZERO) {
      S = ONE;
    } else {
      S = NEGONE;
    }
    Z = S * (X.abs() - SIGMA) * (S + SIGMA / X);
    W = S * Y;
  }

  // Generate the rotation.
  // CALL DLARTGP( Z, W, CS, SN, R ) might seem more natural;
  // reordering the arguments ensures that if Z = 0 then the rotation
  // is by PI/2.

  dlartgp(W, Z, SN, CS, R);
}
