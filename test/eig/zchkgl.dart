// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zggbal.dart';
import 'package:dart_lapack/src/zlange.dart';

Future<void> zchkgl(
  final Nin NIN,
  final Nout NOUT,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const LDA = 20, LDB = 20, LWORK = 6 * LDA;
  const ZERO = 0.0;
  int I, IHIIN, ILOIN, J, KNT, N, NINFO;
  double ANORM, BNORM, EPS, RMAX, VMAX;
  final LMAX = Array<int>(3);
  final LSCALE = Array<double>(LDA),
      LSCLIN = Array<double>(LDA),
      RSCALE = Array<double>(LDA),
      RSCLIN = Array<double>(LDA),
      WORK = Array<double>(LWORK);
  final A = Matrix<Complex>(LDA, LDA),
      AIN = Matrix<Complex>(LDA, LDA),
      B = Matrix<Complex>(LDB, LDB),
      BIN = Matrix<Complex>(LDB, LDB);
  final INFO = Box(0), ILO = Box(0), IHI = Box(0);

  LMAX[1] = 0;
  LMAX[2] = 0;
  LMAX[3] = 0;
  NINFO = 0;
  KNT = 0;
  RMAX = ZERO;

  EPS = dlamch('Precision');
  try {
    while (true) {
      N = await NIN.readInt();
      if (N == 0) break;

      await NIN.readMatrix(A, N, N);
      await NIN.readMatrix(B, N, N);
      (ILOIN, IHIIN) = await NIN.readInt2();

      await NIN.readMatrix(AIN, N, N);
      await NIN.readMatrix(BIN, N, N);

      await NIN.readArray(LSCLIN, N);
      await NIN.readArray(RSCLIN, N);

      ANORM = zlange('M', N, N, A, LDA, WORK);
      BNORM = zlange('M', N, N, B, LDB, WORK);

      KNT++;

      zggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO);

      if (INFO.value != 0) {
        NINFO++;
        LMAX[1] = KNT;
      }

      if (ILO.value != ILOIN || IHI.value != IHIIN) {
        NINFO++;
        LMAX[2] = KNT;
      }

      VMAX = ZERO;
      for (I = 1; I <= N; I++) {
        for (J = 1; J <= N; J++) {
          VMAX = max(VMAX, (A[I][J] - AIN[I][J]).abs());
          VMAX = max(VMAX, (B[I][J] - BIN[I][J]).abs());
        }
      }

      for (I = 1; I <= N; I++) {
        VMAX = max(VMAX, (LSCALE[I] - LSCLIN[I]).abs());
        VMAX = max(VMAX, (RSCALE[I] - RSCLIN[I]).abs());
      }

      VMAX /= EPS * max(ANORM, BNORM);

      if (VMAX > RMAX) {
        LMAX[3] = KNT;
        RMAX = VMAX;
      }
    }
  } catch (_) {}

  NOUT.println(' .. test output of ZGGBAL .. ');

  NOUT.println(' ratio of largest test error              = ${RMAX.d12_3}');
  NOUT.println(' example number where info is not zero    = ${LMAX[1].i4}');
  NOUT.println(' example number where ILO or IHI is wrong = ${LMAX[2].i4}');
  NOUT.println(' example number having largest error      = ${LMAX[3].i4}');
  NOUT.println(' number of examples where info is not 0   = ${NINFO.i4}');
  NOUT.println(' total number of examples tested          = ${KNT.i4}');
}
