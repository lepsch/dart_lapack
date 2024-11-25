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
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import 'ddrvrf1.dart';
import 'ddrvrf2.dart';
import 'ddrvrf3.dart';
import 'ddrvrf4.dart';
import 'ddrvrfp.dart';
import 'derrrfp.dart';

void main() async {
  final NIN = Nin(stdin), NOUT = Nout(stdout);
  const MAXIN = 12;
  const NMAX = 50;
  const MAXRHS = 16;
  const NTYPES = 9;

  final NVAL = Array<int>(MAXIN),
      NSVAL = Array<int>(MAXIN),
      NTVAL = Array<int>(NTYPES);
  final WORKA = Matrix<double>(NMAX, NMAX),
      WORKASAV = Matrix<double>(NMAX, NMAX),
      WORKB = Matrix<double>(NMAX, MAXRHS),
      WORKXACT = Matrix<double>(NMAX, MAXRHS),
      WORKBSAV = Matrix<double>(NMAX, MAXRHS),
      WORKX = Matrix<double>(NMAX, MAXRHS),
      WORKAFAC = Matrix<double>(NMAX, NMAX),
      WORKAINV = Matrix<double>(NMAX, NMAX),
      WORKARF = Array<double>((NMAX * (NMAX + 1)) ~/ 2),
      WORKAP = Array<double>((NMAX * (NMAX + 1)) ~/ 2),
      WORKARFINV = Array<double>((NMAX * (NMAX + 1)) ~/ 2),
      D_WORK_DLATMS = Array<double>(3 * NMAX),
      D_WORK_DPOT01 = Array<double>(NMAX),
      D_TEMP_DPOT02 = Matrix<double>(NMAX, MAXRHS),
      D_TEMP_DPOT03 = Matrix<double>(NMAX, NMAX),
      D_WORK_DLANSY = Array<double>(NMAX),
      D_WORK_DPOT02 = Array<double>(NMAX),
      D_WORK_DPOT03 = Array<double>(NMAX);

  final S1 = dsecnd();
  var FATAL = false;

  // Read a dummy line.

  await NIN.readLine();

  // Report LAPACK version tag (e.g. LAPACK-3.2.0)

  final VERS_MAJOR = Box(0), VERS_MINOR = Box(0), VERS_PATCH = Box(0);
  ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
  NOUT.println(
      '\n Tests of the DOUBLE PRECISION LAPACK RFP routines  LAPACK VERSION ${VERS_MAJOR.value.i1}.${VERS_MINOR.value.i1}.${VERS_PATCH.value.i1}\n\n The following parameter values will be used:');

  // Read the values of N

  var NN = await NIN.readInt();
  if (NN < 1) {
    NOUT.print9996(' NN ', NN, 1);
    NN = 0;
    FATAL = true;
  } else if (NN > MAXIN) {
    NOUT.print9995(' NN ', NN, MAXIN);
    NN = 0;
    FATAL = true;
  }
  await NIN.readArray(NVAL, NN);
  for (var I = 1; I <= NN; I++) {
    if (NVAL[I] < 0) {
      NOUT.print9996(' M  ', NVAL[I], 0);
      FATAL = true;
    } else if (NVAL[I] > NMAX) {
      NOUT.print9995(' M  ', NVAL[I], NMAX);
      FATAL = true;
    }
  }
  if (NN > 0) NOUT.print9993('N   ', NVAL, NN);

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

  // Read the matrix types

  var NNT = await NIN.readInt();
  if (NNT < 1) {
    NOUT.print9996(' NMA', NNT, 1);
    NNT = 0;
    FATAL = true;
  } else if (NNT > NTYPES) {
    NOUT.print9995(' NMA', NNT, NTYPES);
    NNT = 0;
    FATAL = true;
  }
  await NIN.readArray(NTVAL, NNT);
  for (var I = 1; I <= NNT; I++) {
    if (NTVAL[I] < 0) {
      NOUT.print9996('TYPE', NTVAL[I], 0);
      FATAL = true;
    } else if (NTVAL[I] > NTYPES) {
      NOUT.print9995('TYPE', NTVAL[I], NTYPES);
      FATAL = true;
    }
  }
  if (NNT > 0) NOUT.print9993('TYPE', NTVAL, NNT);

  // Read the threshold value for the test ratios.

  final THRESH = await NIN.readDouble();
  NOUT.println(
      '\n Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n');

  // Read the flag that indicates whether to test the error exits.

  final TSTERR = await NIN.readBool();

  if (FATAL) {
    NOUT.println('\n Execution not attempted due to input errors');
    return;
  }

  // Calculate and print the machine dependent constants.

  var EPS = dlamch('Underflow threshold');
  NOUT.print9991('underflow', EPS);
  EPS = dlamch('Overflow threshold');
  NOUT.print9991('overflow ', EPS);
  EPS = dlamch('Epsilon');
  NOUT.print9991('precision', EPS);
  NOUT.println();

  // Test the error exit of:

  if (TSTERR) derrrfp(NOUT);

  // Test the routines: dpftrf, dpftri, dpftrs (as in DDRVPO).
  // This also tests the routines: dtfsm, dtftri, dtfttr, dtrttf.

  ddrvrfp(
      NOUT,
      NN,
      NVAL,
      NNS,
      NSVAL,
      NNT,
      NTVAL,
      THRESH,
      WORKA.asArray(),
      WORKASAV.asArray(),
      WORKAFAC.asArray(),
      WORKAINV.asArray(),
      WORKB.asArray(),
      WORKBSAV.asArray(),
      WORKXACT.asArray(),
      WORKX.asArray(),
      WORKARF,
      WORKARFINV,
      D_WORK_DLATMS,
      D_WORK_DPOT01,
      D_TEMP_DPOT02.asArray(),
      D_TEMP_DPOT03.asArray(),
      D_WORK_DLANSY,
      D_WORK_DPOT02,
      D_WORK_DPOT03);

  // Test the routine: dlansf

  ddrvrf1(NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, D_WORK_DLANSY);

  // Test the conversion routines:
  //   dtfttp, dtpttf, dtfttr, dtrttf, dtrttp and dtpttr.

  ddrvrf2(NOUT, NN, NVAL, WORKA, NMAX, WORKARF, WORKAP, WORKASAV);

  // Test the routine: dtfsm

  ddrvrf3(NOUT, NN, NVAL, THRESH, WORKA, NMAX, WORKARF, WORKAINV, WORKAFAC,
      D_WORK_DLANSY, D_WORK_DPOT03, D_WORK_DPOT01);

  // Test the routine: dsfrk

  ddrvrf4(NOUT, NN, NVAL, THRESH, WORKA, WORKAFAC, NMAX, WORKARF, WORKAINV,
      NMAX, D_WORK_DLANSY);

  await NIN.close();
  final S2 = dsecnd();
  NOUT.println('\n End of tests');
  NOUT.println(' Total time used = ${(S2 - S1).f12_2} seconds\n');
}

extension on Nout {
  void print9996(String s, int actual, int expected) {
    println(
        ' !! Invalid input value: ${s.a4}=${actual.i6}; must be >=${expected.i6}');
  }

  void print9995(String s, int actual, int expected) {
    println(
        ' !! Invalid input value: ${s.a4}=${actual.i6}; must be <=${expected.i6}');
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

  void print9991(String s, double d) {
    println(' Relative machine $s is taken to be${d.d16_6}');
  }
}
