// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtgexc.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'dget51.dart';

Future<void> dget40(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Array<int> NINFO_,
  final Box<int> KNT,
  final Nin NIN,
  final TestDriver test,
  final String group,
  final double THRESH,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NINFO = NINFO_.having();
  const ZERO = 0.0, ONE = 1.0;
  const LDT = 10, LWORK = 100 + 4 * LDT + 16;
  final Q = Matrix<double>(LDT, LDT),
      Z = Matrix<double>(LDT, LDT),
      T = Matrix<double>(LDT, LDT),
      T1 = Matrix<double>(LDT, LDT),
      T2 = Matrix<double>(LDT, LDT),
      S = Matrix<double>(LDT, LDT),
      S1 = Matrix<double>(LDT, LDT),
      S2 = Matrix<double>(LDT, LDT),
      TMP = Matrix<double>(LDT, LDT);
  final RESULT = Array<double>(4), WORK = Array<double>(LWORK);

  final EPS = dlamch('P');
  RMAX.value = ZERO;
  LMAX.value = 0;
  KNT.value = 0;
  NINFO[1] = 0;
  NINFO[2] = 0;

  // Read input data until N=0

  while (true) {
    final (N, IFST, ILST) = await NIN.readInt3();
    if (N == 0) return;
    KNT.value++;
    await NIN.readMatrix(TMP, N, N);
    dlacpy('F', N, N, TMP, LDT, T, LDT);
    dlacpy('F', N, N, TMP, LDT, T1, LDT);
    dlacpy('F', N, N, TMP, LDT, T2, LDT);
    await NIN.readMatrix(TMP, N, N);
    dlacpy('F', N, N, TMP, LDT, S, LDT);
    dlacpy('F', N, N, TMP, LDT, S1, LDT);
    dlacpy('F', N, N, TMP, LDT, S2, LDT);

    final ctx = (
      T: T.copy(),
      T1: T1.copy(),
      T2: T2.copy(),
      S: S.copy(),
      S1: S1.copy(),
      S2: S2.copy(),
    );

    test.group(group, () {
      final (:T, :T1, :T2, :S, :S1, :S2) = ctx;

      test('DTGEXC DGET40 N=$N IFST=$IFST ILST=$ILST', () {
        // IFSTSV = IFST;
        // ILSTSV = ILST;
        final IFST1 = Box(IFST),
            ILST1 = Box(ILST),
            IFST2 = Box(IFST),
            ILST2 = Box(ILST);
        var RES = ZERO;

        // Test without accumulating Q and Z

        dlaset('Full', N, N, ZERO, ONE, Q, LDT);
        dlaset('Full', N, N, ZERO, ONE, Z, LDT);
        dtgexc(false, false, N, T1, LDT, S1, LDT, Q, LDT, Z, LDT, IFST1, ILST1,
            WORK, LWORK, NINFO(1));
        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            if (I == J && Q[I][J] != ONE) RES += ONE / EPS;
            if (I != J && Q[I][J] != ZERO) RES += ONE / EPS;
            if (I == J && Z[I][J] != ONE) RES += ONE / EPS;
            if (I != J && Z[I][J] != ZERO) RES += ONE / EPS;
          }
        }

        // Test with accumulating Q

        dlaset('Full', N, N, ZERO, ONE, Q, LDT);
        dlaset('Full', N, N, ZERO, ONE, Z, LDT);
        dtgexc(true, true, N, T2, LDT, S2, LDT, Q, LDT, Z, LDT, IFST2, ILST2,
            WORK, LWORK, NINFO(2));

        // Compare T1 with T2 and S1 with S2

        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            if (T1[I][J] != T2[I][J]) RES += ONE / EPS;
            if (S1[I][J] != S2[I][J]) RES += ONE / EPS;
          }
        }
        if (IFST1.value != IFST2.value) RES += ONE / EPS;
        if (ILST1.value != ILST2.value) RES += ONE / EPS;
        if (NINFO[1] != NINFO[2]) RES += ONE / EPS;

        // Test orthogonality of Q and Z and backward error on T2 and S2

        dget51(1, N, T, LDT, T2, LDT, Q, LDT, Z, LDT, WORK, RESULT(1));
        dget51(1, N, S, LDT, S2, LDT, Q, LDT, Z, LDT, WORK, RESULT(2));
        dget51(3, N, T, LDT, T2, LDT, Q, LDT, Q, LDT, WORK, RESULT(3));
        dget51(3, N, T, LDT, T2, LDT, Z, LDT, Z, LDT, WORK, RESULT(4));
        RES += RESULT[1] + RESULT[2] + RESULT[3] + RESULT[4];

        if (RES > RMAX.value) {
          RMAX.value = RES;
          LMAX.value = KNT.value;
        }

        test.expect(RES, lessThan(THRESH));

        // Read next matrix pair
      });
    });
  }
}
