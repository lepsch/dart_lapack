// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlasd5(
  final int I,
  final Array<double> D_,
  final Array<double> Z_,
  final Array<double> DELTA_,
  final double RHO,
  final Box<double> DSIGMA,
  final Array<double> WORK,
) {
  final D = D_.having();
  final Z = Z_.having();
  final DELTA = DELTA_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0;
  double B, C, DEL, DELSQ, TAU, W;

  DEL = D[2] - D[1];
  DELSQ = DEL * (D[2] + D[1]);
  if (I == 1) {
    W = ONE +
        FOUR *
            RHO *
            (Z[2] * Z[2] / (D[1] + THREE * D[2]) -
                Z[1] * Z[1] / (THREE * D[1] + D[2])) /
            DEL;
    if (W > ZERO) {
      B = DELSQ + RHO * (Z[1] * Z[1] + Z[2] * Z[2]);
      C = RHO * Z[1] * Z[1] * DELSQ;

      // B > ZERO, always

      // The following TAU is DSIGMA * DSIGMA - D[ 1 ] * D[ 1 ]

      TAU = TWO * C / (B + sqrt((B * B - FOUR * C).abs()));

      // The following TAU is DSIGMA - D[ 1 ]

      TAU /= (D[1] + sqrt(D[1] * D[1] + TAU));
      DSIGMA.value = D[1] + TAU;
      DELTA[1] = -TAU;
      DELTA[2] = DEL - TAU;
      WORK[1] = TWO * D[1] + TAU;
      WORK[2] = (D[1] + TAU) + D[2];
      // DELTA[ 1 ] = -Z[ 1 ] / TAU
      // DELTA[ 2 ] = Z[ 2 ] / ( DEL-TAU )
    } else {
      B = -DELSQ + RHO * (Z[1] * Z[1] + Z[2] * Z[2]);
      C = RHO * Z[2] * Z[2] * DELSQ;

      // The following TAU is DSIGMA * DSIGMA - D[ 2 ] * D[ 2 ]

      if (B > ZERO) {
        TAU = -TWO * C / (B + sqrt(B * B + FOUR * C));
      } else {
        TAU = (B - sqrt(B * B + FOUR * C)) / TWO;
      }

      // The following TAU is DSIGMA - D[ 2 ]

      TAU /= (D[2] + sqrt((D[2] * D[2] + TAU).abs()));
      DSIGMA.value = D[2] + TAU;
      DELTA[1] = -(DEL + TAU);
      DELTA[2] = -TAU;
      WORK[1] = D[1] + TAU + D[2];
      WORK[2] = TWO * D[2] + TAU;
      // DELTA[ 1 ] = -Z[ 1 ] / ( DEL+TAU )
      // DELTA[ 2 ] = -Z[ 2 ] / TAU
    }
    // TEMP = sqrt( DELTA[ 1 ]*DELTA[ 1 ]+DELTA[ 2 ]*DELTA[ 2 ] )
    // DELTA[ 1 ] /= TEMP
    // DELTA[ 2 ] /= TEMP
  } else {
    // Now I=2

    B = -DELSQ + RHO * (Z[1] * Z[1] + Z[2] * Z[2]);
    C = RHO * Z[2] * Z[2] * DELSQ;

    // The following TAU is DSIGMA * DSIGMA - D[ 2 ] * D[ 2 ]

    if (B > ZERO) {
      TAU = (B + sqrt(B * B + FOUR * C)) / TWO;
    } else {
      TAU = TWO * C / (-B + sqrt(B * B + FOUR * C));
    }

    // The following TAU is DSIGMA - D[ 2 ]

    TAU /= (D[2] + sqrt(D[2] * D[2] + TAU));
    DSIGMA.value = D[2] + TAU;
    DELTA[1] = -(DEL + TAU);
    DELTA[2] = -TAU;
    WORK[1] = D[1] + TAU + D[2];
    WORK[2] = TWO * D[2] + TAU;
    // DELTA[ 1 ] = -Z[ 1 ] / ( DEL+TAU )
    // DELTA[ 2 ] = -Z[ 2 ] / TAU
    // TEMP = sqrt( DELTA[ 1 ]*DELTA[ 1 ]+DELTA[ 2 ]*DELTA[ 2 ] )
    // DELTA[ 1 ] /= TEMP
    // DELTA[ 2 ] /= TEMP
  }
}
