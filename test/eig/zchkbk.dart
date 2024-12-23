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
import 'package:dart_lapack/src/zgebak.dart';

Future<void> zchkbk(final Nin NIN, final Nout NOUT) async {
  const LDE = 20;
  const ZERO = 0.0;
  int I, IHI, ILO, J, KNT, N, NINFO;
  double EPS, RMAX, SAFMIN, VMAX, X;
  final LMAX = Array<int>(2);
  final SCALE = Array<double>(LDE);
  final E = Matrix<Complex>(LDE, LDE), EIN = Matrix<Complex>(LDE, LDE);
  final INFO = Box(0);

  LMAX[1] = 0;
  LMAX[2] = 0;
  NINFO = 0;
  KNT = 0;
  RMAX = ZERO;
  EPS = dlamch('E');
  SAFMIN = dlamch('S');

  while (true) {
    (N, ILO, IHI) = await NIN.readInt3();
    if (N == 0) break;

    await NIN.readArray(SCALE, N);
    await NIN.readMatrix(E, N, N);
    await NIN.readMatrix(EIN, N, N);

    KNT++;
    zgebak('B', 'R', N, ILO, IHI, SCALE, N, E, LDE, INFO);

    if (INFO.value != 0) {
      NINFO++;
      LMAX[1] = KNT;
    }

    VMAX = ZERO;
    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        X = (E[I][J] - EIN[I][J]).cabs1() / EPS;
        if (E[I][J].cabs1() > SAFMIN) X /= E[I][J].cabs1();
        VMAX = max(VMAX, X);
      }
    }

    if (VMAX > RMAX) {
      LMAX[2] = KNT;
      RMAX = VMAX;
    }
  }

  NOUT.println(' .. test output of ZGEBAK .. ');
  NOUT.println(' value of largest test error             = ${RMAX.d12_3}');
  NOUT.println(' example number where info is not zero   = ${LMAX[1].i4}');
  NOUT.println(' example number having largest error     = ${LMAX[2].i4}');
  NOUT.println(' number of examples where info is not 0  = ${NINFO.i4}');
  NOUT.println(' total number of examples tested         = ${KNT.i4}');
}
