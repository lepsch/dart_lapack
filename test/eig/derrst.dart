import 'package:lapack/src/box.dart';
import 'package:lapack/src/dopgtr.dart';
import 'package:lapack/src/dopmtr.dart';
import 'package:lapack/src/dorgtr.dart';
import 'package:lapack/src/dormtr.dart';
import 'package:lapack/src/dpteqr.dart';
import 'package:lapack/src/dsbev.dart';
import 'package:lapack/src/dsbev_2stage.dart';
import 'package:lapack/src/dsbevd.dart';
import 'package:lapack/src/dsbevd_2stage.dart';
import 'package:lapack/src/dsbevx.dart';
import 'package:lapack/src/dsbevx_2stage.dart';
import 'package:lapack/src/dsbtrd.dart';
import 'package:lapack/src/dspev.dart';
import 'package:lapack/src/dspevd.dart';
import 'package:lapack/src/dspevx.dart';
import 'package:lapack/src/dsptrd.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dstedc.dart';
import 'package:lapack/src/dstein.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/dstev.dart';
import 'package:lapack/src/dstevd.dart';
import 'package:lapack/src/dstevr.dart';
import 'package:lapack/src/dstevx.dart';
import 'package:lapack/src/dsyev.dart';
import 'package:lapack/src/dsyev_2stage.dart';
import 'package:lapack/src/dsyevd.dart';
import 'package:lapack/src/dsyevd_2stage.dart';
import 'package:lapack/src/dsyevr.dart';
import 'package:lapack/src/dsyevr_2stage.dart';
import 'package:lapack/src/dsyevx.dart';
import 'package:lapack/src/dsyevx_2stage.dart';
import 'package:lapack/src/dsytd2.dart';
import 'package:lapack/src/dsytrd.dart';
import 'package:lapack/src/dsytrd_2stage.dart';
import 'package:lapack/src/dsytrd_sb2st.dart';
import 'package:lapack/src/dsytrd_sy2sb.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../test_driver.dart';
import 'chkxer.dart';
import 'common.dart';

void derrst(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // NMAX has to be at least 3 or LIW may be too small
  const NMAX = 3, LIW = 12 * NMAX, LW = 20 * NMAX;
  int N, NT;
  final I1 = Array<int>(NMAX),
      I2 = Array<int>(NMAX),
      I3 = Array<int>(NMAX),
      IW = Array<int>(LIW);
  final A = Matrix<double>(NMAX, NMAX),
      C = Matrix<double>(NMAX, NMAX),
      Q = Matrix<double>(NMAX, NMAX),
      Z = Matrix<double>(NMAX, NMAX);
  final D = Array<double>(NMAX),
      E = Array<double>(NMAX),
      R = Array<double>(NMAX),
      TAU = Array<double>(NMAX),
      W = Array<double>(LW),
      X = Array<double>(NMAX);
  final INFO = Box(0), M = Box(0), NSPLIT = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();
  final C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
    }
  }
  for (var J = 1; J <= NMAX; J++) {
    D[J] = J.toDouble();
    E[J] = 0.0;
    I1[J] = J;
    I2[J] = J;
    TAU[J] = 1.0;
  }
  OK.value = true;
  NT = 0;

  // Test error exits for the ST path.

  if (lsamen(2, C2, 'ST')) {
    test('DSYTRD', () {
      srnamc.SRNAMT = 'DSYTRD';
      infoc.INFOT = 1;
      dsytrd('/', 0, A, 1, D, E, TAU, W, 1, INFO);
      chkxer('DSYTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrd('U', -1, A, 1, D, E, TAU, W, 1, INFO);
      chkxer('DSYTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrd('U', 2, A, 1, D, E, TAU, W, 1, INFO);
      chkxer('DSYTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsytrd('U', 0, A, 1, D, E, TAU, W, 0, INFO);
      chkxer('DSYTRD', infoc.INFOT, NOUT, LERR, OK);
      NT += 4;
    });

    test('DSYTD2', () {
      srnamc.SRNAMT = 'DSYTD2';
      infoc.INFOT = 1;
      dsytd2('/', 0, A, 1, D, E, TAU, INFO);
      chkxer('DSYTD2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytd2('U', -1, A, 1, D, E, TAU, INFO);
      chkxer('DSYTD2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytd2('U', 2, A, 1, D, E, TAU, INFO);
      chkxer('DSYTD2', infoc.INFOT, NOUT, LERR, OK);
      NT += 3;
    });

    test('DSYTRD_2STAGE', () {
      srnamc.SRNAMT = 'DSYTRD_2STAGE';
      infoc.INFOT = 1;
      dsytrd_2stage('/', 'U', 0, A, 1, D, E, TAU, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsytrd_2stage('H', 'U', 0, A, 1, D, E, TAU, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrd_2stage('N', '/', 0, A, 1, D, E, TAU, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrd_2stage('N', 'U', -1, A, 1, D, E, TAU, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrd_2stage('N', 'U', 2, A, 1, D, E, TAU, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsytrd_2stage('N', 'U', 0, A, 1, D, E, TAU, C.asArray(), 0, W, 1, INFO);
      chkxer('DSYTRD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dsytrd_2stage('N', 'U', 0, A, 1, D, E, TAU, C.asArray(), 1, W, 0, INFO);
      chkxer('DSYTRD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      NT += 7;
    });

    test('DSYTRD_SY2SB', () {
      srnamc.SRNAMT = 'DSYTRD_SY2SB';
      infoc.INFOT = 1;
      dsytrd_sy2sb('/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO);
      chkxer('DSYTRD_SY2SB', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrd_sy2sb('U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO);
      chkxer('DSYTRD_SY2SB', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrd_sy2sb('U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO);
      chkxer('DSYTRD_SY2SB', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrd_sy2sb('U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO);
      chkxer('DSYTRD_SY2SB', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrd_sy2sb('U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO);
      chkxer('DSYTRD_SY2SB', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsytrd_sy2sb('U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO);
      chkxer('DSYTRD_SY2SB', infoc.INFOT, NOUT, LERR, OK);
      NT += 6;
    });

    test('DSYTRD_SB2ST', () {
      srnamc.SRNAMT = 'DSYTRD_SB2ST';
      infoc.INFOT = 1;
      dsytrd_sb2st('/', 'N', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrd_sb2st('Y', '/', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrd_sb2st('Y', 'H', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrd_sb2st('Y', 'N', '/', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrd_sb2st(
          'Y', 'N', 'U', -1, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrd_sb2st(
          'Y', 'N', 'U', 0, -1, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrd_sb2st('Y', 'N', 'U', 0, 1, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsytrd_sb2st('Y', 'N', 'U', 0, 0, A, 1, D, E, C.asArray(), 0, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsytrd_sb2st('Y', 'N', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 0, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DORGTR', () {
      srnamc.SRNAMT = 'DORGTR';
      infoc.INFOT = 1;
      dorgtr('/', 0, A, 1, TAU, W, 1, INFO);
      chkxer('DORGTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dorgtr('U', -1, A, 1, TAU, W, 1, INFO);
      chkxer('DORGTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dorgtr('U', 2, A, 1, TAU, W, 1, INFO);
      chkxer('DORGTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dorgtr('U', 3, A, 3, TAU, W, 1, INFO);
      chkxer('DORGTR', infoc.INFOT, NOUT, LERR, OK);
      NT += 4;
    });

    test('DORMTR', () {
      srnamc.SRNAMT = 'DORMTR';
      infoc.INFOT = 1;
      dormtr('/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dormtr('L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dormtr('L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dormtr('L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dormtr('L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dormtr('L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dormtr('R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dormtr('L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dormtr('L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dormtr('R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO);
      chkxer('DORMTR', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    test('DSPTRD', () {
      srnamc.SRNAMT = 'DSPTRD';
      infoc.INFOT = 1;
      dsptrd('/', 0, A.asArray(), D, E, TAU, INFO);
      chkxer('DSPTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsptrd('U', -1, A.asArray(), D, E, TAU, INFO);
      chkxer('DSPTRD', infoc.INFOT, NOUT, LERR, OK);
      NT += 2;
    });

    test('DOPGTR', () {
      srnamc.SRNAMT = 'DOPGTR';
      infoc.INFOT = 1;
      dopgtr('/', 0, A.asArray(), TAU, Z, 1, W, INFO);
      chkxer('DOPGTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dopgtr('U', -1, A.asArray(), TAU, Z, 1, W, INFO);
      chkxer('DOPGTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dopgtr('U', 2, A.asArray(), TAU, Z, 1, W, INFO);
      chkxer('DOPGTR', infoc.INFOT, NOUT, LERR, OK);
      NT += 3;
    });

    test('DOPMTR', () {
      srnamc.SRNAMT = 'DOPMTR';
      infoc.INFOT = 1;
      dopmtr('/', 'U', 'N', 0, 0, A.asArray(), TAU, C, 1, W, INFO);
      chkxer('DOPMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dopmtr('L', '/', 'N', 0, 0, A.asArray(), TAU, C, 1, W, INFO);
      chkxer('DOPMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dopmtr('L', 'U', '/', 0, 0, A.asArray(), TAU, C, 1, W, INFO);
      chkxer('DOPMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dopmtr('L', 'U', 'N', -1, 0, A.asArray(), TAU, C, 1, W, INFO);
      chkxer('DOPMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dopmtr('L', 'U', 'N', 0, -1, A.asArray(), TAU, C, 1, W, INFO);
      chkxer('DOPMTR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dopmtr('L', 'U', 'N', 2, 0, A.asArray(), TAU, C, 1, W, INFO);
      chkxer('DOPMTR', infoc.INFOT, NOUT, LERR, OK);
      NT += 6;
    });

    test('DPTEQR', () {
      srnamc.SRNAMT = 'DPTEQR';
      infoc.INFOT = 1;
      dpteqr('/', 0, D, E, Z, 1, W, INFO);
      chkxer('DPTEQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dpteqr('N', -1, D, E, Z, 1, W, INFO);
      chkxer('DPTEQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dpteqr('V', 2, D, E, Z, 1, W, INFO);
      chkxer('DPTEQR', infoc.INFOT, NOUT, LERR, OK);
      NT += 3;
    });

    test('DSTEBZ', () {
      srnamc.SRNAMT = 'DSTEBZ';
      infoc.INFOT = 1;
      dstebz('/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dstebz('A', '/', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dstebz('A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dstebz('V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dstebz('I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dstebz('I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dstebz('I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dstebz('I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, D, E, M, NSPLIT, X, I1, I2, W,
          IW, INFO);
      chkxer('DSTEBZ', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    test('DSTEIN', () {
      srnamc.SRNAMT = 'DSTEIN';
      infoc.INFOT = 1;
      dstein(-1, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEIN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dstein(0, D, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEIN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dstein(0, D, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEIN', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dstein(2, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEIN', infoc.INFOT, NOUT, LERR, OK);
      NT += 4;
    });

    test('DSTEQR', () {
      srnamc.SRNAMT = 'DSTEQR';
      infoc.INFOT = 1;
      dsteqr('/', 0, D, E, Z, 1, W, INFO);
      chkxer('DSTEQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsteqr('N', -1, D, E, Z, 1, W, INFO);
      chkxer('DSTEQR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsteqr('V', 2, D, E, Z, 1, W, INFO);
      chkxer('DSTEQR', infoc.INFOT, NOUT, LERR, OK);
      NT += 3;
    });

    test('DSTERF', () {
      srnamc.SRNAMT = 'DSTERF';
      infoc.INFOT = 1;
      dsterf(-1, D, E, INFO);
      chkxer('DSTERF', infoc.INFOT, NOUT, LERR, OK);
      NT++;
    });

    test('DSTEDC', () {
      srnamc.SRNAMT = 'DSTEDC';
      infoc.INFOT = 1;
      dstedc('/', 0, D, E, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dstedc('N', -1, D, E, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dstedc('V', 2, D, E, Z, 1, W, 23, IW, 28, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstedc('N', 1, D, E, Z, 1, W, 0, IW, 1, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstedc('I', 2, D, E, Z, 2, W, 0, IW, 12, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstedc('V', 2, D, E, Z, 2, W, 0, IW, 28, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dstedc('N', 1, D, E, Z, 1, W, 1, IW, 0, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dstedc('I', 2, D, E, Z, 2, W, 19, IW, 0, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dstedc('V', 2, D, E, Z, 2, W, 23, IW, 0, INFO);
      chkxer('DSTEDC', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DSTEVD', () {
      srnamc.SRNAMT = 'DSTEVD';
      infoc.INFOT = 1;
      dstevd('/', 0, D, E, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSTEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dstevd('N', -1, D, E, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSTEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dstevd('V', 2, D, E, Z, 1, W, 19, IW, 12, INFO);
      chkxer('DSTEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstevd('N', 1, D, E, Z, 1, W, 0, IW, 1, INFO);
      chkxer('DSTEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstevd('V', 2, D, E, Z, 2, W, 12, IW, 12, INFO);
      chkxer('DSTEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dstevd('N', 0, D, E, Z, 1, W, 1, IW, 0, INFO);
      chkxer('DSTEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dstevd('V', 2, D, E, Z, 2, W, 19, IW, 11, INFO);
      chkxer('DSTEVD', infoc.INFOT, NOUT, LERR, OK);
      NT += 7;
    });

    test('DSTEV', () {
      srnamc.SRNAMT = 'DSTEV';
      infoc.INFOT = 1;
      dstev('/', 0, D, E, Z, 1, W, INFO);
      chkxer('DSTEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dstev('N', -1, D, E, Z, 1, W, INFO);
      chkxer('DSTEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dstev('V', 2, D, E, Z, 1, W, INFO);
      chkxer('DSTEV', infoc.INFOT, NOUT, LERR, OK);
      NT += 3;
    });

    test('DSTEVX', () {
      srnamc.SRNAMT = 'DSTEVX';
      infoc.INFOT = 1;
      dstevx(
          '/', 'A', 0, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dstevx(
          'N', '/', 0, D, E, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dstevx(
          'N', 'A', -1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dstevx(
          'N', 'V', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstevx(
          'N', 'I', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstevx(
          'N', 'I', 1, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dstevx(
          'N', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dstevx(
          'N', 'I', 1, D, E, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dstevx(
          'V', 'A', 2, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO);
      chkxer('DSTEVX', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DSTEVR', () {
      N = 1;
      srnamc.SRNAMT = 'DSTEVR';
      infoc.INFOT = 1;
      dstevr('/', 'A', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dstevr('V', '/', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dstevr('V', 'A', -1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dstevr('V', 'V', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dstevr('V', 'I', 1, D, E, 0.0, 0.0, 0, 1, 0.0, M, W, Z, 1, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      N = 2;
      dstevr('V', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, W, Z, 1, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      N = 1;
      dstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 0, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 1, IW, X,
          20 * N - 1, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 19;
      dstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 1, IW, X, 20 * N,
          IW(2 * N + 1), 10 * N - 1, INFO);
      chkxer('DSTEVR', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DSYEVD', () {
      srnamc.SRNAMT = 'DSYEVD';
      infoc.INFOT = 1;
      dsyevd('/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyevd('N', '/', 0, A, 1, X, W, 1, IW, 1, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyevd('N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsyevd('N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevd('N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevd('N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevd('V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevd('N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevd('N', 'U', 2, A, 2, X, W, 5, IW, 0, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevd('V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO);
      chkxer('DSYEVD', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    test('DSYEVD_2STAGE', () {
      srnamc.SRNAMT = 'DSYEVD_2STAGE';
      infoc.INFOT = 1;
      dsyevd_2stage('/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsyevd_2stage('V', 'U', 0, A, 1, X, W, 1, IW, 1, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyevd_2stage('N', '/', 0, A, 1, X, W, 1, IW, 1, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyevd_2stage('N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsyevd_2stage('N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevd_2stage('N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevd_2stage('N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // infoc.INFOT = 8
      // CALL DSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
      // CALL CHKXER( 'DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK )
      infoc.INFOT = 10;
      dsyevd_2stage('N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevd_2stage('N', 'U', 2, A, 2, X, W, 25, IW, 0, INFO);
      chkxer('DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // infoc.INFOT = 10
      // CALL DSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
      // CALL CHKXER( 'DSYEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK )
      NT += 9;
    });

    test('DSYEVR', () {
      srnamc.SRNAMT = 'DSYEVR';
      N = 1;
      infoc.INFOT = 1;
      dsyevr('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyevr('V', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyevr('V', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsyevr('V', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsyevr('V', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevr('V', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;

      dsyevr('V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 0, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 0, INFO);
      chkxer('DSYEVR', infoc.INFOT, NOUT, LERR, OK);
      NT += 11;
    });

    test('DSYEVR_2STAGE', () {
      srnamc.SRNAMT = 'DSYEVR_2STAGE';
      N = 1;
      infoc.INFOT = 1;
      dsyevr_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsyevr_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyevr_2stage('N', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyevr_2stage('N', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1,
          IW, Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsyevr_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1,
          IW, Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsyevr_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevr_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevr_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 0, IW(2 * N + 1), 10 * N, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 20;
      dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW,
          Q.asArray(), 26 * N, IW(2 * N + 1), 0, INFO);
      chkxer('DSYEVR_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      NT += 12;
    });

    test('DSYEV', () {
      srnamc.SRNAMT = 'DSYEV';
      infoc.INFOT = 1;
      dsyev('/', 'U', 0, A, 1, X, W, 1, INFO);
      chkxer('DSYEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyev('N', '/', 0, A, 1, X, W, 1, INFO);
      chkxer('DSYEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyev('N', 'U', -1, A, 1, X, W, 1, INFO);
      chkxer('DSYEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsyev('N', 'U', 2, A, 1, X, W, 3, INFO);
      chkxer('DSYEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyev('N', 'U', 1, A, 1, X, W, 1, INFO);
      chkxer('DSYEV', infoc.INFOT, NOUT, LERR, OK);
      NT += 5;
    });

    test('DSYEV_2STAGE', () {
      srnamc.SRNAMT = 'DSYEV_2STAGE';
      infoc.INFOT = 1;
      dsyev_2stage('/', 'U', 0, A, 1, X, W, 1, INFO);
      chkxer('DSYEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsyev_2stage('V', 'U', 0, A, 1, X, W, 1, INFO);
      chkxer('DSYEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyev_2stage('N', '/', 0, A, 1, X, W, 1, INFO);
      chkxer('DSYEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyev_2stage('N', 'U', -1, A, 1, X, W, 1, INFO);
      chkxer('DSYEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsyev_2stage('N', 'U', 2, A, 1, X, W, 3, INFO);
      chkxer('DSYEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyev_2stage('N', 'U', 1, A, 1, X, W, 1, INFO);
      chkxer('DSYEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      NT += 6;
    });

    test('DSYEVX', () {
      srnamc.SRNAMT = 'DSYEVX';
      infoc.INFOT = 1;
      dsyevx('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyevx('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyevx('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW,
          I3, INFO);
      infoc.INFOT = 4;
      dsyevx('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsyevx('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevx('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 8, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevx('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 16, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 8, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dsyevx('V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dsyevx('V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW,
          I3, INFO);
      chkxer('DSYEVX', infoc.INFOT, NOUT, LERR, OK);
      NT += 12;
    });

    test('DSYEVX_2STAGE', () {
      srnamc.SRNAMT = 'DSYEVX_2STAGE';
      infoc.INFOT = 1;
      dsyevx_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          1, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsyevx_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          1, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyevx_2stage('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W,
          1, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyevx_2stage('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          1, IW, I3, INFO);
      infoc.INFOT = 4;
      dsyevx_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          1, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsyevx_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          16, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsyevx_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          8, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          8, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W,
          8, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevx_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W,
          16, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W,
          8, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 15;
      dsyevx_2stage('N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 0, W,
          16, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 17;
      dsyevx_2stage('N', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          0, IW, I3, INFO);
      chkxer('DSYEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      NT += 13;
    });

    test('DSPEVD', () {
      srnamc.SRNAMT = 'DSPEVD';
      infoc.INFOT = 1;
      dspevd('/', 'U', 0, A.asArray(), X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dspevd('N', '/', 0, A.asArray(), X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dspevd('N', 'U', -1, A.asArray(), X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dspevd('V', 'U', 2, A.asArray(), X, Z, 1, W, 23, IW, 12, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dspevd('N', 'U', 1, A.asArray(), X, Z, 1, W, 0, IW, 1, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dspevd('N', 'U', 2, A.asArray(), X, Z, 1, W, 3, IW, 1, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dspevd('V', 'U', 2, A.asArray(), X, Z, 2, W, 16, IW, 12, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dspevd('N', 'U', 1, A.asArray(), X, Z, 1, W, 1, IW, 0, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dspevd('N', 'U', 2, A.asArray(), X, Z, 1, W, 4, IW, 0, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dspevd('V', 'U', 2, A.asArray(), X, Z, 2, W, 23, IW, 11, INFO);
      chkxer('DSPEVD', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    test('DSPEV', () {
      srnamc.SRNAMT = 'DSPEV';
      infoc.INFOT = 1;
      dspev('/', 'U', 0, A.asArray(), W, Z, 1, X, INFO);
      chkxer('DSPEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dspev('N', '/', 0, A.asArray(), W, Z, 1, X, INFO);
      chkxer('DSPEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dspev('N', 'U', -1, A.asArray(), W, Z, 1, X, INFO);
      chkxer('DSPEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dspev('V', 'U', 2, A.asArray(), W, Z, 1, X, INFO);
      chkxer('DSPEV', infoc.INFOT, NOUT, LERR, OK);
      NT += 4;
    });

    test('DSPEVX', () {
      srnamc.SRNAMT = 'DSPEVX';
      infoc.INFOT = 1;
      dspevx('/', 'A', 'U', 0, A.asArray(), 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dspevx('N', '/', 'U', 0, A.asArray(), 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dspevx('N', 'A', '/', 0, A.asArray(), 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      infoc.INFOT = 4;
      dspevx('N', 'A', 'U', -1, A.asArray(), 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dspevx('N', 'V', 'U', 1, A.asArray(), 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dspevx('N', 'I', 'U', 1, A.asArray(), 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dspevx('N', 'I', 'U', 1, A.asArray(), 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dspevx('N', 'I', 'U', 2, A.asArray(), 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dspevx('N', 'I', 'U', 1, A.asArray(), 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 14;
      dspevx('V', 'A', 'U', 2, A.asArray(), 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W,
          IW, I3, INFO);
      chkxer('DSPEVX', infoc.INFOT, NOUT, LERR, OK);
      NT += 10;
    });

    // Test error exits for the SB path.
  } else if (lsamen(2, C2, 'SB')) {
    test('DSBTRD', () {
      srnamc.SRNAMT = 'DSBTRD';
      infoc.INFOT = 1;
      dsbtrd('/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO);
      chkxer('DSBTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsbtrd('N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO);
      chkxer('DSBTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsbtrd('N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO);
      chkxer('DSBTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsbtrd('N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO);
      chkxer('DSBTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsbtrd('N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO);
      chkxer('DSBTRD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsbtrd('V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO);
      chkxer('DSBTRD', infoc.INFOT, NOUT, LERR, OK);
      NT += 6;
    });

    test('DSYTRD_SB2ST', () {
      srnamc.SRNAMT = 'DSYTRD_SB2ST';
      infoc.INFOT = 1;
      dsytrd_sb2st('/', 'N', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrd_sb2st('N', '/', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrd_sb2st('N', 'H', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrd_sb2st('N', 'N', '/', 0, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrd_sb2st(
          'N', 'N', 'U', -1, 0, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrd_sb2st(
          'N', 'N', 'U', 0, -1, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrd_sb2st('N', 'N', 'U', 0, 1, A, 1, D, E, C.asArray(), 1, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsytrd_sb2st('N', 'N', 'U', 0, 0, A, 1, D, E, C.asArray(), 0, W, 1, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsytrd_sb2st('N', 'N', 'U', 0, 0, A, 1, D, E, C.asArray(), 1, W, 0, INFO);
      chkxer('DSYTRD_SB2ST', infoc.INFOT, NOUT, LERR, OK);
      NT += 9;
    });

    test('DSBEVD', () {
      srnamc.SRNAMT = 'DSBEVD';
      infoc.INFOT = 1;
      dsbevd('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsbevd('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsbevd('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsbevd('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsbevd('N', 'U', 2, 1, A, 1, X, Z, 1, W, 4, IW, 1, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsbevd('V', 'U', 2, 1, A, 2, X, Z, 1, W, 25, IW, 12, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, IW, 1, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsbevd('N', 'U', 2, 0, A, 1, X, Z, 1, W, 3, IW, 1, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsbevd('V', 'U', 2, 0, A, 1, X, Z, 2, W, 18, IW, 12, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, IW, 0, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsbevd('V', 'U', 2, 0, A, 1, X, Z, 2, W, 25, IW, 11, INFO);
      chkxer('DSBEVD', infoc.INFOT, NOUT, LERR, OK);
      NT += 11;
    });

    test('DSBEVD_2STAGE', () {
      srnamc.SRNAMT = 'DSBEVD_2STAGE';
      infoc.INFOT = 1;
      dsbevd_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsbevd_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsbevd_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsbevd_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsbevd_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsbevd_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 4, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // infoc.INFOT = 9
      // CALL DSBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
      // $                                      25, IW, 12, INFO )
      // CALL CHKXER( 'DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK )
      infoc.INFOT = 11;
      dsbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsbevd_2stage('N', 'U', 2, 0, A, 1, X, Z, 1, W, 3, IW, 1, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // infoc.INFOT = 11
      // CALL DSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
      // $                                      18, IW, 12, INFO )
      // CALL CHKXER( 'DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK )
      infoc.INFOT = 13;
      dsbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, IW, 0, INFO);
      chkxer('DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // infoc.INFOT = 13
      // CALL DSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
      // $                                      25, IW, 11, INFO )
      // CALL CHKXER( 'DSBEVD_2STAGE', infoc.INFOT, NOUT, LERR, OK )
      // NT += 12
      NT += 9;
    });

    test('DSBEV', () {
      srnamc.SRNAMT = 'DSBEV';
      infoc.INFOT = 1;
      dsbev('/', 'U', 0, 0, A, 1, X, Z, 1, W, INFO);
      chkxer('DSBEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsbev('N', '/', 0, 0, A, 1, X, Z, 1, W, INFO);
      chkxer('DSBEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsbev('N', 'U', -1, 0, A, 1, X, Z, 1, W, INFO);
      chkxer('DSBEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsbev('N', 'U', 0, -1, A, 1, X, Z, 1, W, INFO);
      chkxer('DSBEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsbev('N', 'U', 2, 1, A, 1, X, Z, 1, W, INFO);
      chkxer('DSBEV', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsbev('V', 'U', 2, 0, A, 1, X, Z, 1, W, INFO);
      chkxer('DSBEV', infoc.INFOT, NOUT, LERR, OK);
      NT += 6;
    });

    test('DSBEV_2STAGE', () {
      srnamc.SRNAMT = 'DSBEV_2STAGE';
      infoc.INFOT = 1;
      dsbev_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsbev_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsbev_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsbev_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsbev_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsbev_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsbev_2stage('N', 'U', 2, 0, A, 1, X, Z, 0, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsbev_2stage('N', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO);
      chkxer('DSBEV_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      NT += 8;
    });

    test('DSBEVX', () {
      srnamc.SRNAMT = 'DSBEVX';
      infoc.INFOT = 1;
      dsbevx('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsbevx('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsbevx('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsbevx('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsbevx('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsbevx('N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsbevx('V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsbevx('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dsbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dsbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsbevx('N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 18;
      dsbevx('V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1,
          W, IW, I3, INFO);
      chkxer('DSBEVX', infoc.INFOT, NOUT, LERR, OK);
      NT += 13;
    });

    test('DSBEVX_2STAGE', () {
      srnamc.SRNAMT = 'DSBEVX_2STAGE';
      infoc.INFOT = 1;
      dsbevx_2stage('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 1;
      dsbevx_2stage('V', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsbevx_2stage('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsbevx_2stage('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsbevx_2stage('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsbevx_2stage('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsbevx_2stage('N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // infoc.INFOT = 9
      // CALL DSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0,
      // $          0.0, 0, 0, 0.0, M, X, Z, 2, W, 0, IW, I3, INFO )
      // CALL CHKXER( 'DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK )
      infoc.INFOT = 11;
      dsbevx_2stage('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dsbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dsbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsbevx_2stage('N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 13;
      dsbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // infoc.INFOT = 18
      // CALL DSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0,
      // $          0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO )
      // CALL CHKXER( 'DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK )
      infoc.INFOT = 20;
      dsbevx_2stage('N', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X,
          Z, 1, W, 0, IW, I3, INFO);
      chkxer('DSBEVX_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      // NT += 15
      NT += 13;
    });
  }

  // Print a summary line.

  if (OK.value) {
    infoc.NOUT.println(
        ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)');
  } else {
    infoc.NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
