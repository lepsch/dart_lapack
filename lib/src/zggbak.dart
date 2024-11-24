// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zggbak(
  final String JOB,
  final String SIDE,
  final int N,
  final int ILO,
  final int IHI,
  final Array<double> LSCALE_,
  final Array<double> RSCALE_,
  final int M,
  final Matrix<Complex> V_,
  final int LDV,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having(ld: LDV);
  final LSCALE = LSCALE_.having();
  final RSCALE = RSCALE_.having();
  bool LEFTV, RIGHTV;
  int I, K;

  // Test the input parameters

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
  } else if (ILO < 1) {
    INFO.value = -4;
  } else if (N == 0 && IHI == 0 && ILO != 1) {
    INFO.value = -4;
  } else if (N > 0 && (IHI < ILO || IHI > max(1, N))) {
    INFO.value = -5;
  } else if (N == 0 && ILO == 1 && IHI != 0) {
    INFO.value = -5;
  } else if (M < 0) {
    INFO.value = -8;
  } else if (LDV < max(1, N)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZGGBAK', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;
  if (M == 0) return;
  if (lsame(JOB, 'N')) return;

  if (ILO != IHI) {
    // Backward balance

    if (lsame(JOB, 'S') || lsame(JOB, 'B')) {
      // Backward transformation on right eigenvectors

      if (RIGHTV) {
        for (I = ILO; I <= IHI; I++) {
          zdscal(M, RSCALE[I], V(I, 1).asArray(), LDV);
        }
      }

      // Backward transformation on left eigenvectors

      if (LEFTV) {
        for (I = ILO; I <= IHI; I++) {
          zdscal(M, LSCALE[I], V(I, 1).asArray(), LDV);
        }
      }
    }

    // Backward permutation
  }
  if (lsame(JOB, 'P') || lsame(JOB, 'B')) {
    // Backward permutation on right eigenvectors

    if (RIGHTV) {
      if (ILO != 1) {
        for (I = ILO - 1; I >= 1; I--) {
          K = RSCALE[I].toInt();
          if (K == I) continue;
          zswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
        }
      }
      if (IHI != N) {
        for (I = IHI + 1; I <= N; I++) {
          K = RSCALE[I].toInt();
          if (K == I) continue;
          zswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
        }
      }
    }

    // Backward permutation on left eigenvectors

    if (LEFTV) {
      if (ILO != 1) {
        for (I = ILO - 1; I >= 1; I--) {
          K = LSCALE[I].toInt();
          if (K == I) continue;
          zswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
        }
      }
      if (IHI == N) {
        for (I = IHI + 1; I <= N; I++) {
          K = LSCALE[I].toInt();
          if (K == I) continue;
          zswap(M, V(I, 1).asArray(), LDV, V(K, 1).asArray(), LDV);
        }
      }
    }
  }
}
