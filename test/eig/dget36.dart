// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtrexc.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'dhst01.dart';

Future<void> dget36(
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
  const LDT = 10, LWORK = 2 * LDT * LDT;
  final Q = Matrix<double>(LDT, LDT),
      T1 = Matrix<double>(LDT, LDT),
      T2 = Matrix<double>(LDT, LDT),
      TMP = Matrix<double>(LDT, LDT);
  final RESULT = Array<double>(2), WORK = Array<double>(LWORK);

  final EPS = dlamch('P');

  RMAX.value = ZERO;
  LMAX.value = 0;
  KNT.value = 0;
  NINFO[1] = 0;
  NINFO[2] = 0;
  NINFO[3] = 0;

  // Read input data until N=0

  while (true) {
    final (N, IFST, ILST) = await NIN.readInt3();
    if (N == 0) return;
    KNT.value++;
    await NIN.readMatrix(TMP, N, N);
    dlacpy('F', N, N, TMP, LDT, T1, LDT);
    dlacpy('F', N, N, TMP, LDT, T2, LDT);

    final ctx = (
      TMP: TMP.copy(),
      T1: T1.copy(),
      T2: T2.copy(),
    );

    test.group(group, () {
      test('DGET36', () {
        final (:TMP, :T1, :T2) = ctx;
        final IFSTSV = IFST;
        final ILSTSV = ILST;
        final IFST1 = Box(IFST),
            IFST2 = Box(IFST),
            ILST1 = Box(ILST),
            ILST2 = Box(ILST);
        var RES = ZERO;
        final INFO1 = Box(0), INFO2 = Box(0);

        // Test without accumulating Q

        dlaset('Full', N, N, ZERO, ONE, Q, LDT);
        dtrexc('N', N, T1, LDT, Q, LDT, IFST1, ILST1, WORK, INFO1);
        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            if (I == J && Q[I][J] != ONE) RES += ONE / EPS;
            if (I != J && Q[I][J] != ZERO) RES += ONE / EPS;
          }
        }

        // Test with accumulating Q

        dlaset('Full', N, N, ZERO, ONE, Q, LDT);
        dtrexc('V', N, T2, LDT, Q, LDT, IFST2, ILST2, WORK, INFO2);

        // Compare T1 with T2

        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            if (T1[I][J] != T2[I][J]) RES += ONE / EPS;
          }
        }
        if (IFST1.value != IFST2.value) RES += ONE / EPS;
        if (ILST1.value != ILST2.value) RES += ONE / EPS;
        if (INFO1.value != INFO2.value) RES += ONE / EPS;

        // Test for successful reordering of T2

        if (INFO2.value != 0) {
          NINFO[INFO2.value]++;
        } else {
          if ((IFST2.value - IFSTSV).abs() > 1) RES += ONE / EPS;
          if ((ILST2.value - ILSTSV).abs() > 1) RES += ONE / EPS;
        }

        // Test for small residual, and orthogonality of Q

        dhst01(N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RESULT);
        RES += RESULT[1] + RESULT[2];

        // Test for T2 being in Schur form

        var LOC = 1;
        do {
          if (T2[LOC + 1][LOC] != ZERO) {
            // 2 by 2 block

            if (T2[LOC][LOC + 1] == ZERO ||
                T2[LOC][LOC] != T2[LOC + 1][LOC + 1] ||
                sign(ONE, T2[LOC][LOC + 1]) == sign(ONE, T2[LOC + 1][LOC])) {
              RES += ONE / EPS;
            }
            for (var I = LOC + 2; I <= N; I++) {
              if (T2[I][LOC] != ZERO) RES += ONE / RES;
              if (T2[I][LOC + 1] != ZERO) RES += ONE / RES;
            }
            LOC += 2;
          } else {
            // 1 by 1 block

            for (var I = LOC + 1; I <= N; I++) {
              if (T2[I][LOC] != ZERO) RES += ONE / RES;
            }
            LOC++;
          }
        } while (LOC < N);

        if (RES > RMAX.value) {
          RMAX.value = RES;
          LMAX.value = KNT.value;
        }

        test.expect(RES, lessThan(THRESH));
        test.expect(NINFO[3], 0);
      });
    });
  }
}
