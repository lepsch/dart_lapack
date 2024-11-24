// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';

void dlarrc(
  final String JOBT,
  final int N,
  final double VL,
  final double VU,
  final Array<double> D_,
  final Array<double> E_,
  final double PIVMIN,
  final Box<int> EIGCNT,
  final Box<int> LCNT,
  final Box<int> RCNT,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  const ZERO = 0.0;
  int I;
  bool MATT;
  double LPIVOT, RPIVOT, SL, SU, TMP, TMP2;

  INFO.value = 0;
  LCNT.value = 0;
  RCNT.value = 0;
  EIGCNT.value = 0;

  // Quick return if possible

  if (N <= 0) {
    return;
  }

  MATT = lsame(JOBT, 'T');

  if (MATT) {
    // Sturm sequence count on T
    LPIVOT = D[1] - VL;
    RPIVOT = D[1] - VU;
    if (LPIVOT <= ZERO) {
      LCNT.value++;
    }
    if (RPIVOT <= ZERO) {
      RCNT.value++;
    }
    for (I = 1; I <= N - 1; I++) {
      TMP = pow(E[I], 2).toDouble();
      LPIVOT = (D[I + 1] - VL) - TMP / LPIVOT;
      RPIVOT = (D[I + 1] - VU) - TMP / RPIVOT;
      if (LPIVOT <= ZERO) {
        LCNT.value++;
      }
      if (RPIVOT <= ZERO) {
        RCNT.value++;
      }
    }
  } else {
    // Sturm sequence count on L D L^T
    SL = -VL;
    SU = -VU;
    for (I = 1; I <= N - 1; I++) {
      LPIVOT = D[I] + SL;
      RPIVOT = D[I] + SU;
      if (LPIVOT <= ZERO) {
        LCNT.value++;
      }
      if (RPIVOT <= ZERO) {
        RCNT.value++;
      }
      TMP = E[I] * D[I] * E[I];

      TMP2 = TMP / LPIVOT;
      if (TMP2 == ZERO) {
        SL = TMP - VL;
      } else {
        SL = SL * TMP2 - VL;
      }

      TMP2 = TMP / RPIVOT;
      if (TMP2 == ZERO) {
        SU = TMP - VU;
      } else {
        SU = SU * TMP2 - VU;
      }
    }
    LPIVOT = D[N] + SL;
    RPIVOT = D[N] + SU;
    if (LPIVOT <= ZERO) {
      LCNT.value++;
    }
    if (RPIVOT <= ZERO) {
      RCNT.value++;
    }
  }
  EIGCNT.value = RCNT.value - LCNT.value;
}
