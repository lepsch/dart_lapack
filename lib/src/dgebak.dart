// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/dswap.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dgebak(
  final String JOB,
  final String SIDE,
  final int N,
  final int ILO,
  final int IHI,
  final Array<double> SCALE_,
  final int M,
  final Matrix<double> V_,
  final int LDV,
  final Box<int> INFO,
) {
  final SCALE = SCALE_.having();
  final V = V_.having(ld: LDV);
  const ONE = 1.0;
  bool LEFTV, RIGHTV;
  int I, II, K;
  double S;

  // Decode and Test the input parameters

  RIGHTV = lsame(SIDE, 'R');
  LEFTV = lsame(SIDE, 'L');

  INFO.value = 0;
  if (!lsame(JOB, 'N') &&
      !lsame(JOB, 'P') &&
      !lsame(JOB, 'S') &&
      !lsame(JOB, 'B')) {
    INFO.value = -1;
  } else if (!RIGHTV && !LEFTV) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -4;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -5;
  } else if (M < 0) {
    INFO.value = -7;
  } else if (LDV < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DGEBAK', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;
  if (M == 0) return;
  if (lsame(JOB, 'N')) return;

  if (ILO != IHI) {
    // Backward balance

    if (lsame(JOB, 'S') || lsame(JOB, 'B')) {
      if (RIGHTV) {
        for (I = ILO; I <= IHI; I++) {
          S = SCALE[I];
          dscal(M, S, V(I, 1).asArray(), LDV);
        }
      }

      if (LEFTV) {
        for (I = ILO; I <= IHI; I++) {
          S = ONE / SCALE[I];
          dscal(M, S, V(I, 1).asArray(), LDV);
        }
      }
    }
  }

  // Backward permutation
  //
  // For  I = ILO-1 step -1 until 1,
  //          IHI+1 step 1 until N do --

  if (lsame(JOB, 'P') || lsame(JOB, 'B')) {
    if (RIGHTV) {
      for (II = 1; II <= N; II++) {
        I = II;
        if (I >= ILO && I <= IHI) continue;
        if (I < ILO) I = ILO - II;
        K = SCALE[I].toInt();
        if (K == I) continue;
        dswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
      }
    }

    if (LEFTV) {
      for (II = 1; II <= N; II++) {
        I = II;
        if (I >= ILO && I <= IHI) continue;
        if (I < ILO) I = ILO - II;
        K = SCALE[I].toInt();
        if (K == I) continue;
        dswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
      }
    }
  }
}
