import 'package:lapack/src/box.dart';
import 'package:lapack/src/dtrexc.dart';
import 'package:lapack/src/dtrsen.dart';
import 'package:lapack/src/dtrsna.dart';
import 'package:lapack/src/dtrsyl.dart';
import 'package:lapack/src/dtrsyl3.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'chkxer.dart';
import 'common.dart';

void derrec(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4, ONE = 1.0, ZERO = 0.0;
  int I, J, NT;
  final SCALE = Box(0.0);
  final SEL = Array<bool>(NMAX);
  final IWORK = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      B = Matrix<double>(NMAX, NMAX),
      C = Matrix<double>(NMAX, NMAX);
  final S = Array<double>(NMAX),
      SEP = Array<double>(NMAX),
      WI = Array<double>(NMAX),
      WORK = Array<double>(NMAX),
      WR = Array<double>(NMAX);
  final INFO = Box(0), M = Box(0), IFST = Box(0), ILST = Box(0);

  infoc.NOUT = NUNIT;
  infoc.OK.value = true;
  NT = 0;

  // Initialize A, B and SEL

  for (J = 1; J <= NMAX; J++) {
    for (I = 1; I <= NMAX; I++) {
      A[I][J] = ZERO;
      B[I][J] = ZERO;
    }
  }
  for (I = 1; I <= NMAX; I++) {
    A[I][I] = ONE;
    SEL[I] = true;
  }

  // Test DTRSYL

  srnamc.SRNAMT = 'DTRSYL';
  infoc.INFOT = 1;
  dtrsyl('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtrsyl('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtrsyl('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtrsyl('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dtrsyl('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dtrsyl('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  dtrsyl('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  dtrsyl('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, INFO);
  chkxer('DTRSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT = NT + 8;

  // Test DTRSYL3

  srnamc.SRNAMT = 'DTRSYL3';
  infoc.INFOT = 1;
  dtrsyl3('X', 'N', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtrsyl3('N', 'X', 1, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtrsyl3('N', 'N', 0, 0, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtrsyl3('N', 'N', 1, -1, 0, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dtrsyl3('N', 'N', 1, 0, -1, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dtrsyl3('N', 'N', 1, 2, 0, A, 1, B, 1, C, 2, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  dtrsyl3('N', 'N', 1, 0, 2, A, 1, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  dtrsyl3('N', 'N', 1, 2, 0, A, 2, B, 1, C, 1, SCALE, IWORK, NMAX, WORK, NMAX,
      INFO);
  chkxer('DTRSYL3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT = NT + 8;

  // Test DTREXC

  srnamc.SRNAMT = 'DTREXC';
  IFST.value = 1;
  ILST.value = 1;
  infoc.INFOT = 1;
  dtrexc('X', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtrexc('N', -1, A, 1, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ILST.value = 2;
  dtrexc('N', 2, A, 1, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dtrexc('V', 2, A, 2, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  IFST.value = 0;
  ILST.value = 1;
  dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  IFST.value = 2;
  dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  IFST.value = 1;
  ILST.value = 0;
  dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  ILST.value = 2;
  dtrexc('V', 1, A, 1, B, 1, IFST, ILST, WORK, INFO);
  chkxer('DTREXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT = NT + 8;

  // Test DTRSNA

  srnamc.SRNAMT = 'DTRSNA';
  infoc.INFOT = 1;
  dtrsna('X', 'A', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(1), 1,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtrsna('B', 'X', SEL, 0, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(1), 1,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtrsna('B', 'A', SEL, -1, A, 1, B, 1, C, 1, S, SEP, 1, M, WORK.asMatrix(1), 1,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dtrsna('V', 'A', SEL, 2, A, 1, B, 1, C, 1, S, SEP, 2, M, WORK.asMatrix(2), 2,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  dtrsna('B', 'A', SEL, 2, A, 2, B, 1, C, 2, S, SEP, 2, M, WORK.asMatrix(2), 2,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  dtrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 1, S, SEP, 2, M, WORK.asMatrix(2), 2,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 13;
  dtrsna('B', 'A', SEL, 1, A, 1, B, 1, C, 1, S, SEP, 0, M, WORK.asMatrix(1), 1,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 13;
  dtrsna('B', 'S', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 1, M, WORK.asMatrix(2), 2,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 16;
  dtrsna('B', 'A', SEL, 2, A, 2, B, 2, C, 2, S, SEP, 2, M, WORK.asMatrix(1), 1,
      IWORK, INFO);
  chkxer('DTRSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT = NT + 9;

  // Test DTRSEN

  SEL[1] = false;
  srnamc.SRNAMT = 'DTRSEN';
  infoc.INFOT = 1;
  dtrsen('X', 'N', SEL, 0, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK, 1,
      IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtrsen('N', 'X', SEL, 0, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK, 1,
      IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtrsen('N', 'N', SEL, -1, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK,
      1, IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dtrsen('N', 'N', SEL, 2, A, 1, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK, 2,
      IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  dtrsen('N', 'V', SEL, 2, A, 2, B, 1, WR, WI, M, S.box(1), SEP.box(1), WORK, 1,
      IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 15;
  dtrsen('N', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S.box(1), SEP.box(1), WORK, 0,
      IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 15;
  dtrsen('E', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S.box(1), SEP.box(1), WORK, 1,
      IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 15;
  dtrsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S.box(1), SEP.box(1), WORK, 3,
      IWORK, 2, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 17;
  dtrsen('E', 'V', SEL, 2, A, 2, B, 2, WR, WI, M, S.box(1), SEP.box(1), WORK, 1,
      IWORK, 0, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 17;
  dtrsen('V', 'V', SEL, 3, A, 3, B, 3, WR, WI, M, S.box(1), SEP.box(1), WORK, 4,
      IWORK, 1, INFO);
  chkxer('DTRSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  NT = NT + 10;

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT.println(
      ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)',
    );
  } else {
    infoc.NOUT.println(
      ' *** ${PATH.a3} routines failed the tests of the error exits ***',
    );
  }
}
