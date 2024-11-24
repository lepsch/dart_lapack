// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'zerrec.dart';
import 'zget35.dart';
import 'zget36.dart';
import 'zget37.dart';
import 'zget38.dart';
import 'zsyl01.dart';

Future<void> zchkec(
  final double THRESH,
  final bool TSTERR,
  final Nin NIN,
  final Nout NOUT,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool OK;
  String PATH;
  int NTESTS = 0;
  double EPS, SFMIN;
  final FTRSYL = Array<int>(3),
      ITRSYL = Array<int>(2),
      LTRSEN = Array<int>(3),
      LTRSNA = Array<int>(3),
      NTRSEN = Array<int>(3),
      NTRSNA = Array<int>(3);
  final RTRSEN = Array<double>(3),
      RTRSNA = Array<double>(3),
      RTRSYL = Array<double>(2);
  final LTRSYL = Box(0),
      NTRSYL = Box(0),
      KTRSYL = Box(0),
      KTRSYL3 = Box(0),
      KTRSNA = Box(0),
      KTRSEN = Box(0),
      LTREXC = Box(0),
      NTREXC = Box(0),
      KTREXC = Box(0);
  final RTREXC = Box(0.0);

  PATH = '${'Zomplex precision'[0]}EC';
  EPS = dlamch('P');
  SFMIN = dlamch('S');
  NOUT.println(
      ' Tests of the Nonsymmetric eigenproblem condition estimation routines\n ZTRSYL, ZTREXC, ZTRSNA, ZTRSEN\n');
  NOUT.println(
      ' Relative machine precision (EPS) = ${EPS.d16_6}\n Safe minimum (SFMIN)             = ${SFMIN.d16_6}\n');
  NOUT.println(
      ' Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n\n');

  // Test error exits if TSTERR is true;

  if (TSTERR) zerrec(PATH, NOUT);

  OK = true;
  await zget35(RTRSYL(1), LTRSYL, NTRSYL, KTRSYL, NIN);
  if (RTRSYL[1] > THRESH) {
    OK = false;
    NOUT.println(
        ' Error in ZTRSYL: RMAX =${RTRSYL[1].d12_3}\n LMAX = ${LTRSYL.value.i8} NINFO=${NTRSYL.value.i8} KNT=${KTRSYL.value.i8}');
  }

  zsyl01(THRESH, FTRSYL, RTRSYL, ITRSYL, KTRSYL3);
  if (FTRSYL[1] > 0) {
    OK = false;
    NOUT.println(
        'Error in ZTRSYL: ${FTRSYL[1].i8} tests fail the threshold.\nMaximum test ratio =${RTRSYL[1].d12_3} threshold =${THRESH.d12_3}');
  }
  if (FTRSYL[2] > 0) {
    OK = false;
    NOUT.println(
        'Error in ZTRSYL3: ${FTRSYL[2].i8} tests fail the threshold.\nMaximum test ratio =${RTRSYL[2].d12_3} threshold =${THRESH.d12_3}');
  }
  if (FTRSYL[3] > 0) {
    OK = false;
    NOUT.println(
        'ZTRSYL and ZTRSYL3 compute an inconsistent scale factor in ${FTRSYL[3].i8} tests.');
  }

  await zget36(RTREXC, LTREXC, NTREXC, KTREXC, NIN);
  if (RTREXC.value > THRESH || NTREXC.value > 0) {
    OK = false;
    NOUT.println(
        ' Error in ZTREXC: RMAX =${RTREXC.value.d12_3}\n LMAX = ${LTREXC.value.i8} NINFO=${NTREXC.value.i8} KNT=${KTREXC.value.i8}');
  }

  await zget37(RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN);
  if (RTRSNA[1] > THRESH ||
      RTRSNA[2] > THRESH ||
      NTRSNA[1] != 0 ||
      NTRSNA[2] != 0 ||
      NTRSNA[3] != 0) {
    OK = false;
    NOUT.println(
        ' Error in ZTRSNA: RMAX =${RTRSNA.d12_3(3)}\n LMAX = ${LTRSNA.i8(3)} NINFO=${NTRSNA.i8(3)} KNT=${KTRSNA.value.i8}');
  }

  await zget38(RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN);
  if (RTRSEN[1] > THRESH ||
      RTRSEN[2] > THRESH ||
      NTRSEN[1] != 0 ||
      NTRSEN[2] != 0 ||
      NTRSEN[3] != 0) {
    OK = false;
    NOUT.println(
        ' Error in ZTRSEN: RMAX =${RTRSEN.d12_3(3)}\n LMAX = ${LTRSEN.i8(3)} NINFO=${NTRSEN.i8(3)} KNT=${KTRSEN.value.i8}');
  }

  NTESTS =
      KTRSYL.value + KTRSYL3.value + KTREXC.value + KTRSNA.value + KTRSEN.value;
  if (OK) {
    NOUT.println(
        '\n All tests for ${PATH.a3} routines passed the threshold ( ${NTESTS.i6} tests run)');
  }
}
