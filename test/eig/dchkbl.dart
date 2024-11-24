// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgebal.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';

Future<void> dchkbl(
  final Nin NIN,
  final Nout NOUT,
  final TestDriver test,
  final String group,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const LDA = 20;
  const ZERO = 0.0, ONE = 1.0;
  final A = Matrix<double>(LDA, LDA), AIN = Matrix<double>(LDA, LDA);
  final SCALE = Array<double>(LDA), SCALIN = Array<double>(LDA);

  final LMAX = Array.fromList([0, 0, 0]);
  var NINFO = 0;
  var KNT = 0;
  var RMAX = ZERO;
  var VMAX = ZERO;
  final SFMIN = dlamch('S');

  while (true) {
    final N = await NIN.readInt();
    if (N == 0) break;
    await NIN.readMatrix(A, N, N);
    final (ILOIN, IHIIN) = await NIN.readInt2();
    await NIN.readMatrix(AIN, N, N);
    await NIN.readArray(SCALIN, N);

    final ctx = (A: A.copy(), AIN: AIN.copy(), SCALIN: SCALIN.copy());
    test.group(group, () {
      final (:A, :AIN, :SCALIN) = ctx;
      test('DCHKBL (N=$N, ILOIN=$ILOIN, IHIIN=$IHIIN)', () {
        final INFO = Box(0), IHI = Box(0), ILO = Box(0);

        // ANORM = dlange('M', N, N, A, LDA, DUMMY);
        KNT++;

        dgebal('B', N, A, LDA, ILO, IHI, SCALE, INFO);

        if (INFO.value != 0) {
          NINFO++;
          LMAX[1] = KNT;
        }

        if (ILO.value != ILOIN || IHI.value != IHIIN) {
          NINFO++;
          LMAX[2] = KNT;
        }

        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            final TEMP = max(max(A[I][J], AIN[I][J]), SFMIN);
            VMAX = max(VMAX, (A[I][J] - AIN[I][J]).abs() / TEMP);
          }
        }

        for (var I = 1; I <= N; I++) {
          final TEMP = max(max(SCALE[I], SCALIN[I]), SFMIN);
          VMAX = max(VMAX, (SCALE[I] - SCALIN[I]).abs() / TEMP);
        }

        test.expect(VMAX, lessThanOrEqualTo(ONE));
        if (VMAX > RMAX) {
          LMAX[3] = KNT;
          RMAX = VMAX;
        }
      });
    });
  }

  NOUT.println(' .. test output of DGEBAL .. ');
  NOUT.println(' value of largest test error            = ${RMAX.d12_3}');
  NOUT.println(' example number where info is not zero  = ${LMAX[1].i4}');
  NOUT.println(' example number where ILO or IHI wrong  = ${LMAX[2].i4}');
  NOUT.println(' example number having largest error    = ${LMAX[3].i4}');
  NOUT.println(' number of examples where info is not 0 = ${NINFO.i4}');
  NOUT.println(' total number of examples tested        = ${KNT.i4}');
}
