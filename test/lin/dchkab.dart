// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';
import 'dart:math';

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/install/dsecnd.dart';
import 'package:dart_lapack/src/install/ilaver.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/install/slamch.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import 'alareq.dart';
import 'ddrvab.dart';
import 'ddrvac.dart';
import 'derrab.dart';
import 'derrac.dart';

void main() async {
  final NIN = Nin(stdin), NOUT = Nout(stdout);
  const NMAX = 132;
  const MAXIN = 12;
  const MAXRHS = 16;
  const MATMAX = 30;
  const LDAMAX = NMAX;
  final DOTYPE = Array<bool>(MATMAX);
  final IWORK = Array<int>(NMAX),
      MVAL = Array<int>(MAXIN),
      NSVAL = Array<int>(MAXIN);
  final A = Matrix<double>(LDAMAX * NMAX, 2),
      B = Matrix<double>(NMAX * MAXRHS, 2),
      RWORK = Array<double>(NMAX),
      WORK = Array<double>(NMAX * MAXRHS * 2),
      SWORK = Array<double>(NMAX * (NMAX + MAXRHS));
  const INTSTR = '0123456789';

  final S1 = dsecnd();
  final LDA = NMAX;
  var FATAL = false;

  // Read a dummy line.

  await NIN.readLine();

  // Report values of parameters.

  final VERS_MAJOR = Box(0), VERS_MINOR = Box(0), VERS_PATCH = Box(0);
  ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
  NOUT.println(
      ' Tests of the DOUBLE PRECISION LAPACK DSGESV/DSPOSV routines \n LAPACK VERSION ${VERS_MAJOR.value.i1}.${VERS_MINOR.value.i1}.${VERS_PATCH.value.i1}\n\n The following parameter values will be used:');

  // Read the values of M

  var NM = await NIN.readInt();
  if (NM < 1) {
    NOUT.print9996(' NM ', NM, 1);
    NM = 0;
    FATAL = true;
  } else if (NM > MAXIN) {
    NOUT.print9995(' NM ', NM, MAXIN);
    NM = 0;
    FATAL = true;
  }
  await NIN.readArray(MVAL, NM);
  for (var I = 1; I <= NM; I++) {
    if (MVAL[I] < 0) {
      NOUT.print9996(' M  ', MVAL[I], 0);
      FATAL = true;
    } else if (MVAL[I] > NMAX) {
      NOUT.print9995(' M  ', MVAL[I], NMAX);
      FATAL = true;
    }
  }
  if (NM > 0) NOUT.print9993('M   ', MVAL, NM);

  // Read the values of NRHS

  var NNS = await NIN.readInt();
  if (NNS < 1) {
    NOUT.print9996(' NNS', NNS, 1);
    NNS = 0;
    FATAL = true;
  } else if (NNS > MAXIN) {
    NOUT.print9995(' NNS', NNS, MAXIN);
    NNS = 0;
    FATAL = true;
  }
  await NIN.readArray(NSVAL, NNS);
  for (var I = 1; I <= NNS; I++) {
    if (NSVAL[I] < 0) {
      NOUT.print9996('NRHS', NSVAL[I], 0);
      FATAL = true;
    } else if (NSVAL[I] > MAXRHS) {
      NOUT.print9995('NRHS', NSVAL[I], MAXRHS);
      FATAL = true;
    }
  }
  if (NNS > 0) NOUT.print9993('NRHS', NSVAL, NNS);

  // Read the threshold value for the test ratios.

  final THRESH = await NIN.readDouble();
  NOUT.println(
      '\n Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n');

  // Read the flag that indicates whether to test the driver routine.

  final TSTDRV = await NIN.readBool();

  // Read the flag that indicates whether to test the error exits.

  final TSTERR = await NIN.readBool();

  if (FATAL) {
    NOUT.println('\n Execution not attempted due to input errors');
    return;
  }

  // Calculate and print the machine dependent constants.

  var SEPS = slamch('Underflow threshold');
  NOUT.print9991('(single precision) underflow', SEPS);
  SEPS = slamch('Overflow threshold');
  NOUT.print9991('(single precision) overflow ', SEPS);
  SEPS = slamch('Epsilon');
  NOUT.print9991('(single precision) precision', SEPS);
  NOUT.println();

  var EPS = dlamch('Underflow threshold');
  NOUT.print9991('(double precision) underflow', EPS);
  EPS = dlamch('Overflow threshold');
  NOUT.print9991('(double precision) overflow ', EPS);
  EPS = dlamch('Epsilon');
  NOUT.print9991('(double precision) precision', EPS);
  NOUT.println();

  while (true) {
    // Read a test path and the number of matrix types to use.
    final String ALINE;
    try {
      ALINE = await NIN.readLine();
    } on EOF catch (_) {
      break;
    }
    final PATH = ALINE.substring(0, 3);
    var NMATS = MATMAX;
    var I = 3;
    do {
      I++;
      if (I > 72) {
        NMATS = MATMAX;
        break;
      }
    } while (ALINE[I - 1] == ' ');

    if (I <= 72) {
      NMATS = 0;
      while (true) {
        final C1 = ALINE[I - 1];
        var isDigit = false;
        var IC = 0;
        for (var K = 1; K <= 10; K++) {
          if (C1 == INTSTR[K - 1]) {
            IC = K - 1;
            isDigit = true;
            break;
          }
        }
        if (!isDigit) break;

        NMATS = NMATS * 10 + IC;
        I++;
        if (I > 72) break;
      }
    }
    final C1 = PATH[0];
    final C2 = PATH.substring(1, 3);
    // final NRHS = NSVAL[1];

    // Check first character for correct precision.

    if (!lsame(C1, 'Double precision')) {
      NOUT.println('\n ${PATH.a6} routines were not tested');
    } else if (NMATS <= 0) {
      // Check for a positive number of tests requested.

      NOUT.print9989(PATH);
      break;
    } else if (lsamen(2, C2, 'GE')) {
      // GE:  general matrices

      final NTYPES = 11;
      await alareq('DGE', NMATS, DOTYPE, NTYPES, NIN, NOUT);

      // Test the error exits

      if (TSTERR) derrab(NOUT);

      if (TSTDRV) {
        ddrvab(
            DOTYPE,
            NM,
            MVAL,
            NNS,
            NSVAL,
            THRESH,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            WORK,
            RWORK,
            SWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989('DSGESV');
      }
    } else if (lsamen(2, C2, 'PO')) {
      // PO:  positive definite matrices

      final NTYPES = 9;
      await alareq('DPO', NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTERR) derrac(NOUT);

      if (TSTDRV) {
        ddrvac(
            DOTYPE,
            NM,
            MVAL,
            NNS,
            NSVAL,
            THRESH,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            WORK,
            RWORK,
            SWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else {}

    // Go back to get another input line.
  }

  // Branch to this line when the last record is read.

  await NIN.close();
  final S2 = dsecnd();
  NOUT.println('\n End of tests');
  NOUT.println(' Total time used = ${(S2 - S1).f12_2} seconds\n');
}

extension on Nout {
  void print9996(String s, int actual, int expected) {
    println(
        ' Invalid input value: ${s.a4}=${actual.i6}; must be >=${expected.i6}');
  }

  void print9995(String s, int actual, int expected) {
    println(
        ' Invalid input value: ${s.a4}=${actual.i6}; must be <=${expected.i6}');
  }

  void print9993(final String s, final Array<int> a, int n) {
    var prefix = '    $s:  ';
    var i = 1;
    while (n > 0) {
      println('$prefix${a(i).i6(min(n, 10))}');
      prefix = ' ' * 11;
      n -= 10;
      i += 10;
    }
  }

  void print9991(String s, double EPS) {
    println(' Relative machine $s is taken to be${EPS.d16_6}');
  }

  void print9989(String s) {
    println('\n ${s.a6} driver routines were not tested');
  }
}
