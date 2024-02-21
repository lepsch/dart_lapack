import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'zerrec.dart';
import 'zget35.dart';
import 'zget36.dart';
import 'zget37.dart';
import 'zget38.dart';
import 'zsyl01.dart';

void zchkec(
  final double THRESH,
  final bool TSTERR,
  final Nin NIN,
  final Nout NOUT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool OK;
  String PATH;
  int KTREXC = 0,
      KTRSEN = 0,
      KTRSNA = 0,
      KTRSYL = 0,
      KTRSYL3 = 0,
      LTREXC = 0,
      LTRSYL = 0,
      NTESTS = 0,
      NTREXC = 0,
      NTRSYL = 0;
  double EPS, RTREXC = 0, SFMIN;
  final FTRSYL = Array<int>(3),
      ITRSYL = Array<int>(2),
      LTRSEN = Array<int>(3),
      LTRSNA = Array<int>(3),
      NTRSEN = Array<int>(3),
      NTRSNA = Array<int>(3);
  final RTRSEN = Array<double>(3),
      RTRSNA = Array<double>(3),
      RTRSYL = Array<double>(2);
  // ..
  // .. External Subroutines ..
  // EXTERNAL ZERREC, ZGET35, ZGET36, ZGET37, ZGET38, ZSYL01
  // ..
  // .. External Functions ..
  //- double             DLAMCH;
  // EXTERNAL DLAMCH

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
  zget35(RTRSYL[1], LTRSYL, NTRSYL, KTRSYL, NIN);
  if (RTRSYL[1] > THRESH) {
    OK = false;
    NOUT.println(
        ' Error in ZTRSYL: RMAX =${RTRSYL[1].d12_3}\n LMAX = ${LTRSYL.i8} NINFO=${NTRSYL.i8} KNT=${KTRSYL.i8}');
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

  zget36(RTREXC, LTREXC, NTREXC, KTREXC, NIN);
  if (RTREXC > THRESH || NTREXC > 0) {
    OK = false;
    NOUT.println(
        ' Error in ZTREXC: RMAX =${RTREXC.d12_3}\n LMAX = ${LTREXC.i8} NINFO=${NTREXC.i8} KNT=${KTREXC.i8}');
  }

  zget37(RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN);
  if (RTRSNA[1] > THRESH ||
      RTRSNA[2] > THRESH ||
      NTRSNA[1] != 0 ||
      NTRSNA[2] != 0 ||
      NTRSNA[3] != 0) {
    OK = false;
    NOUT.println(
        ' Error in ZTRSNA: RMAX =${RTRSNA.d12_3(3)}\n LMAX = ${LTRSNA.i8(3)} NINFO=${NTRSNA.i8(3)} KNT=${KTRSNA.i8}');
  }

  zget38(RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN);
  if (RTRSEN[1] > THRESH ||
      RTRSEN[2] > THRESH ||
      NTRSEN[1] != 0 ||
      NTRSEN[2] != 0 ||
      NTRSEN[3] != 0) {
    OK = false;
    NOUT.println(
        ' Error in ZTRSEN: RMAX =${RTRSEN.d12_3(3)}\n LMAX = ${LTRSEN.i8(3)} NINFO=${NTRSEN.i8(3)} KNT=${KTRSEN.i8}');
  }

  NTESTS = KTRSYL + KTRSYL3 + KTREXC + KTRSNA + KTRSEN;
  if (OK)
    NOUT.println(
        '\n All tests for ${PATH.a3} routines passed the threshold ( ${NTESTS.i6} tests run)');
}
