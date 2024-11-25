// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';

void drotm(
  final int N,
  final Array<double> DX_,
  final int INCX,
  final Array<double> DY_,
  final int INCY,
  final Array<double> DPARAM,
) {
  final DX = DX_.having();
  final DY = DY_.having();
  final (ZERO, TWO) = (0.0, 2.0);

  final DFLAG = DPARAM[1];
  if (N <= 0 || (DFLAG + TWO == ZERO)) return;
  if (INCX == INCY && INCX > 0) {
    final NSTEPS = N * INCX;
    if (DFLAG < ZERO) {
      final DH11 = DPARAM[2],
          DH12 = DPARAM[4],
          DH21 = DPARAM[3],
          DH22 = DPARAM[5];
      for (var I = 1; I <= NSTEPS; I += INCX) {
        final W = DX[I], Z = DY[I];
        DX[I] = W * DH11 + Z * DH12;
        DY[I] = W * DH21 + Z * DH22;
      }
    } else if (DFLAG == ZERO) {
      final DH12 = DPARAM[4], DH21 = DPARAM[3];
      for (var I = 1; I <= NSTEPS; I += INCX) {
        final W = DX[I], Z = DY[I];
        DX[I] = W + Z * DH12;
        DY[I] = W * DH21 + Z;
      }
    } else {
      final DH11 = DPARAM[2], DH22 = DPARAM[5];
      for (var I = 1; I <= NSTEPS; I += INCX) {
        final W = DX[I], Z = DY[I];
        DX[I] = W * DH11 + Z;
        DY[I] = -W + DH22 * Z;
      }
    }
  } else {
    var KX = 1, KY = 1;
    if (INCX < 0) KX = 1 + (1 - N) * INCX;
    if (INCY < 0) KY = 1 + (1 - N) * INCY;

    if (DFLAG < ZERO) {
      final DH11 = DPARAM[2],
          DH12 = DPARAM[4],
          DH21 = DPARAM[3],
          DH22 = DPARAM[5];
      for (var I = 1; I <= N; I++) {
        final W = DX[KX], Z = DY[KY];
        DX[KX] = W * DH11 + Z * DH12;
        DY[KY] = W * DH21 + Z * DH22;
        KX += INCX;
        KY += INCY;
      }
    } else if (DFLAG == ZERO) {
      final DH12 = DPARAM[4], DH21 = DPARAM[3];
      for (var I = 1; I <= N; I++) {
        final W = DX[KX], Z = DY[KY];
        DX[KX] = W + Z * DH12;
        DY[KY] = W * DH21 + Z;
        KX += INCX;
        KY += INCY;
      }
    } else {
      final DH11 = DPARAM[2], DH22 = DPARAM[5];
      for (var I = 1; I <= N; I++) {
        final W = DX[KX], Z = DY[KY];
        DX[KX] = W * DH11 + Z;
        DY[KY] = -W + DH22 * Z;
        KX += INCX;
        KY += INCY;
      }
    }
  }
}
