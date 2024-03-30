import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';
import 'derrec.dart';
import 'dget31.dart';
import 'dget32.dart';
import 'dget33.dart';
import 'dget34.dart';
import 'dget35.dart';
import 'dget36.dart';
import 'dget37.dart';
import 'dget38.dart';
import 'dget39.dart';
import 'dget40.dart';
import 'dsyl01.dart';

Future<void> dchkec(
  final double THRESH,
  final bool TSTERR,
  final Nin NIN,
  final Nout NOUT,
  final TestDriver test,
  final String group,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool OK;
  int NTESTS = 0;

  final PATH = '${'Double precision'[0]}EC';
  final EPS = dlamch('P');
  final SFMIN = dlamch('S');

  // print header information

  NOUT.println(
      ' Tests of the Nonsymmetric eigenproblem condition estimation routines\n DLALN2, DLASY2, DLANV2, DLAEXC, DTRSYL, DTREXC, DTRSNA, DTRSEN, DLAQTR, DTGEXC\n');
  NOUT.println(
      ' Relative machine precision (EPS) = ${EPS.d16_6}\n Safe minimum (SFMIN)             = ${SFMIN.d16_6}\n');
  NOUT.println(
      ' Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n\n');

  OK = true;

  // Test error exits if TSTERR is true;
  test.group(group, () {
    test.group('error exits', () {
      if (TSTERR) derrec(PATH, NOUT, test);
    });

    test('DLALN2', () {
      final RLALN2 = Box(0.0),
          LLALN2 = Box(0),
          NLALN2 = Array<int>(2),
          KLALN2 = Box(0);
      dget31(RLALN2, LLALN2, NLALN2, KLALN2);
      NTESTS += KLALN2.value;

      test.expect(RLALN2.value, lessThan(THRESH));
      test.expect(NLALN2[1], 0);
      if (RLALN2.value > THRESH || NLALN2[1] != 0) {
        OK = false;
        NOUT.println(
            ' Error in DLALN2: RMAX =${RLALN2.value.d12_3}\n LMAX = ${LLALN2.value.i8} NINFO=${NLALN2.i8(2)} KNT=${KLALN2.value.i8}');
      }
    });

    test('DLASY2', () {
      final RLASY2 = Box(0.0),
          LLASY2 = Box(0),
          NLASY2 = Box(0),
          KLASY2 = Box(0);
      dget32(RLASY2, LLASY2, NLASY2, KLASY2);
      NTESTS += KLASY2.value;

      test.expect(RLASY2.value, lessThan(THRESH));
      if (RLASY2.value > THRESH) {
        OK = false;
        NOUT.println(
            ' Error in DLASY2: RMAX =${RLASY2.value.d12_3}\n LMAX = ${LLASY2.value.i8} NINFO=${NLASY2.value.i8} KNT=${KLASY2.value.i8}');
      }
    });

    test('DLANV2', () {
      final RLANV2 = Box(0.0),
          LLANV2 = Box(0),
          NLANV2 = Box(0),
          KLANV2 = Box(0);
      dget33(RLANV2, LLANV2, NLANV2, KLANV2);
      NTESTS += KLANV2.value;

      test.expect(RLANV2.value, lessThan(THRESH));
      test.expect(NLANV2.value, 0);
      if (RLANV2.value > THRESH || NLANV2.value != 0) {
        OK = false;
        NOUT.println(
            ' Error in DLANV2: RMAX =${RLANV2.value.d12_3}\n LMAX = ${LLANV2.value.i8} NINFO=${NLANV2.value.i8} KNT=${KLANV2.value.i8}');
      }
    });

    test('DLAEXC', () {
      final RLAEXC = Box(0.0),
          LLAEXC = Box(0),
          NLAEXC = Array<int>(2),
          KLAEXC = Box(0);
      dget34(RLAEXC, LLAEXC, NLAEXC, KLAEXC);
      NTESTS += KLAEXC.value;

      test.expect(RLAEXC.value, lessThan(THRESH));
      test.expect(NLAEXC[2], 0);
      if (RLAEXC.value > THRESH || NLAEXC[2] != 0) {
        OK = false;
        NOUT.println(
            ' Error in DLAEXC: RMAX =${RLAEXC.value.d12_3}\n LMAX = ${LLAEXC.value.i8} NINFO=${NLAEXC.i8(2)} KNT=${KLAEXC.value.i8}');
      }
    });

    test.group('DTRSYL', () {
      final RTRSYL = Array<double>(2),
          LTRSYL = Box(0),
          NTRSYL = Box(0),
          KTRSYL = Box(0);
      dget35(RTRSYL(1), LTRSYL, NTRSYL, KTRSYL, test, THRESH);
      NTESTS += KTRSYL.value;

      if (RTRSYL[1] > THRESH) {
        OK = false;
        NOUT.println(
            ' Error in DTRSYL: RMAX =${RTRSYL[1].d12_3}\n LMAX = ${LTRSYL.value.i8} NINFO=${NTRSYL.value.i8} KNT=${KTRSYL.value.i8}');
      }
    });

    test.group('DSYL01', () {
      final FTRSYL = Array<int>(3),
          RTRSYL = Array<double>(2),
          ITRSYL = Array<int>(2),
          KTRSYL3 = Box(0);
      dsyl01(THRESH, FTRSYL, RTRSYL, ITRSYL, KTRSYL3, test);
      NTESTS += KTRSYL3.value;

      if (FTRSYL[1] > 0) {
        OK = false;
        NOUT.println(
            'Error in DTRSYL: ${FTRSYL[1].i8} tests fail the threshold.\nMaximum test ratio =${RTRSYL[1].d12_3} threshold =${THRESH.d12_3}');
      }

      if (FTRSYL[2] > 0) {
        OK = false;
        NOUT.println(
            'Error in DTRSYL3: ${FTRSYL[2].i8} tests fail the threshold.\nMaximum test ratio =${RTRSYL[2].d12_3} threshold =${THRESH.d12_3}');
      }

      if (FTRSYL[3] > 0) {
        OK = false;
        NOUT.println(
            'DTRSYL and DTRSYL3 compute an inconsistent result factor in ${FTRSYL[3].i8} tests.');
      }
    });
  });

  {
    final RTREXC = Box(0.0),
        LTREXC = Box(0),
        NTREXC = Array<int>(3),
        KTREXC = Box(0);
    await dget36(RTREXC, LTREXC, NTREXC, KTREXC, NIN, test, group, THRESH);
    NTESTS += KTREXC.value;

    if (RTREXC.value > THRESH || NTREXC[3] > 0) {
      OK = false;
      NOUT.println(
          ' Error in DTREXC: RMAX =${RTREXC.value.d12_3}\n LMAX = ${LTREXC.value.i8} NINFO=${NTREXC.i8(3)} KNT=${KTREXC.value.i8}');
    }
  }

  {
    final RTRSNA = Array<double>(3),
        LTRSNA = Array<int>(3),
        NTRSNA = Array<int>(3),
        KTRSNA = Box(0);
    await dget37(RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN, test, group, THRESH);
    NTESTS += KTRSNA.value;

    if (RTRSNA[1] > THRESH ||
        RTRSNA[2] > THRESH ||
        NTRSNA[1] != 0 ||
        NTRSNA[2] != 0 ||
        NTRSNA[3] != 0) {
      OK = false;
      NOUT.println(
          ' Error in DTRSNA: RMAX =${RTRSNA.d12_3(3)}\n LMAX = ${LTRSNA.i8(3)} NINFO=${NTRSNA.i8(3)} KNT=${KTRSNA.value.i8}');
    }
  }

  {
    final RTRSEN = Array<double>(3),
        LTRSEN = Array<int>(3),
        NTRSEN = Array<int>(3),
        KTRSEN = Box(0);
    await dget38(RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN, test, group, THRESH);
    NTESTS += KTRSEN.value;

    if (RTRSEN[1] > THRESH ||
        RTRSEN[2] > THRESH ||
        NTRSEN[1] != 0 ||
        NTRSEN[2] != 0 ||
        NTRSEN[3] != 0) {
      OK = false;
      NOUT.println(
          ' Error in DTRSEN: RMAX =${RTRSEN.d12_3(3)}\n LMAX = ${LTRSEN.i8(3)} NINFO=${NTRSEN.i8(3)} KNT=${KTRSEN.value.i8}');
    }
  }

  test.group(group, () {
    final RLAQTR = Box(0.0), LLAQTR = Box(0), NLAQTR = Box(0), KLAQTR = Box(0);
    dget39(RLAQTR, LLAQTR, NLAQTR, KLAQTR, test, THRESH);
    NTESTS += KLAQTR.value;

    if (RLAQTR.value > THRESH) {
      OK = false;
      NOUT.println(
          ' Error in DLAQTR: RMAX =${RLAQTR.value.d12_3}\n LMAX = ${LLAQTR.value.i8} NINFO=${NLAQTR.value.i8} KNT=${KLAQTR.value.i8}');
    }
  });

  {
    final RTGEXC = Box(0.0),
        LTGEXC = Box(0),
        NTGEXC = Array<int>(2),
        KTGEXC = Box(0);
    await dget40(RTGEXC, LTGEXC, NTGEXC, KTGEXC, NIN, test, group, THRESH);
    NTESTS += KTGEXC.value;

    if (RTGEXC.value > THRESH) {
      OK = false;
      NOUT.println(
          ' Error in DTGEXC: RMAX =${RTGEXC.value.d12_3}\n LMAX = ${LTGEXC.value.i8} NINFO=${NTGEXC.i8(2)} KNT=${KTGEXC.value.i8}');
    }
  }

  if (OK) {
    NOUT.println(
        ' All tests for ${PATH.a3} routines passed the threshold ( ${NTESTS.i6} tests run');
  }
}
