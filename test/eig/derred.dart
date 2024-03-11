import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgees.dart';
import 'package:lapack/src/dgeesx.dart';
import 'package:lapack/src/dgeev.dart';
import 'package:lapack/src/dgeevx.dart';
import 'package:lapack/src/dgejsv.dart';
import 'package:lapack/src/dgesdd.dart';
import 'package:lapack/src/dgesvd.dart';
import 'package:lapack/src/dgesvdq.dart';
import 'package:lapack/src/dgesvdx.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'chkxer.dart';
import 'common.dart';
import 'dslect.dart';

void derred(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4, ONE = 1.0, ZERO = 0.0;
  String C2;
  int I, J, NT;
  final ABNRM = Box(0.0);
  final B = Array<bool>(NMAX);
  final IW = Array<int>(2 * NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      U = Matrix<double>(NMAX, NMAX),
      VL = Matrix<double>(NMAX, NMAX),
      VR = Matrix<double>(NMAX, NMAX),
      VT = Matrix<double>(NMAX, NMAX);
  final R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      S = Array<double>(NMAX),
      W = Array<double>(10 * NMAX),
      WI = Array<double>(NMAX),
      WR = Array<double>(NMAX);
  final INFO = Box(0), NS = Box(0), IHI = Box(0), ILO = Box(0), SDIM = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println('');
  C2 = PATH.substring(1, 3);

  // Initialize A

  for (J = 1; J <= NMAX; J++) {
    for (I = 1; I <= NMAX; I++) {
      A[I][J] = ZERO;
    }
  }
  for (I = 1; I <= NMAX; I++) {
    A[I][I] = ONE;
  }
  infoc.OK.value = true;
  NT = 0;

  if (lsamen(2, C2, 'EV')) {
    // Test DGEEV

    srnamc.SRNAMT = 'DGEEV ';
    infoc.INFOT = 1;
    dgeev('X', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO);
    chkxer('DGEEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgeev('N', 'X', 0, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO);
    chkxer('DGEEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgeev('N', 'N', -1, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO);
    chkxer('DGEEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgeev('N', 'N', 2, A, 1, WR, WI, VL, 1, VR, 1, W, 6, INFO);
    chkxer('DGEEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    dgeev('V', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, W, 8, INFO);
    chkxer('DGEEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    dgeev('N', 'V', 2, A, 2, WR, WI, VL, 1, VR, 1, W, 8, INFO);
    chkxer('DGEEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    dgeev('V', 'V', 1, A, 1, WR, WI, VL, 1, VR, 1, W, 3, INFO);
    chkxer('DGEEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;
  } else if (lsamen(2, C2, 'ES')) {
    // Test DGEES

    srnamc.SRNAMT = 'DGEES ';
    infoc.INFOT = 1;
    dgees('X', 'N', dslect, 0, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO);
    chkxer('DGEES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgees('N', 'X', dslect, 0, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO);
    chkxer('DGEES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgees('N', 'S', dslect, -1, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO);
    chkxer('DGEES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgees('N', 'S', dslect, 2, A, 1, SDIM, WR, WI, VL, 1, W, 6, B, INFO);
    chkxer('DGEES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    dgees('V', 'S', dslect, 2, A, 2, SDIM, WR, WI, VL, 1, W, 6, B, INFO);
    chkxer('DGEES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    dgees('N', 'S', dslect, 1, A, 1, SDIM, WR, WI, VL, 1, W, 2, B, INFO);
    chkxer('DGEES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 6;
  } else if (lsamen(2, C2, 'VX')) {
    // Test DGEEVX

    srnamc.SRNAMT = 'DGEEVX';
    infoc.INFOT = 1;
    dgeevx('X', 'N', 'N', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 1, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgeevx('N', 'X', 'N', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 1, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgeevx('N', 'N', 'X', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 1, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgeevx('N', 'N', 'N', 'X', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 1, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgeevx('N', 'N', 'N', 'N', -1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 1, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgeevx('N', 'N', 'N', 'N', 2, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 1, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    dgeevx('N', 'V', 'N', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 6, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    dgeevx('N', 'N', 'V', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 6, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 21;
    dgeevx('N', 'N', 'N', 'N', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 1, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 21;
    dgeevx('N', 'V', 'N', 'N', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 2, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 21;
    dgeevx('N', 'N', 'V', 'V', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S,
        ABNRM, R1, R2, W, 3, IW, INFO);
    chkxer('DGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 11;
  } else if (lsamen(2, C2, 'SX')) {
    // Test DGEESX

    srnamc.SRNAMT = 'DGEESX';
    infoc.INFOT = 1;
    dgeesx('X', 'N', dslect, 'N', 0, A, 1, SDIM, WR, WI, VL, 1, R1.box(1),
        R2.box(1), W, 1, IW, 1, B, INFO);
    chkxer('DGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgeesx('N', 'X', dslect, 'N', 0, A, 1, SDIM, WR, WI, VL, 1, R1.box(1),
        R2.box(1), W, 1, IW, 1, B, INFO);
    chkxer('DGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgeesx('N', 'N', dslect, 'X', 0, A, 1, SDIM, WR, WI, VL, 1, R1.box(1),
        R2.box(1), W, 1, IW, 1, B, INFO);
    chkxer('DGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgeesx('N', 'N', dslect, 'N', -1, A, 1, SDIM, WR, WI, VL, 1, R1.box(1),
        R2.box(1), W, 1, IW, 1, B, INFO);
    chkxer('DGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgeesx('N', 'N', dslect, 'N', 2, A, 1, SDIM, WR, WI, VL, 1, R1.box(1),
        R2.box(1), W, 6, IW, 1, B, INFO);
    chkxer('DGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    dgeesx('V', 'N', dslect, 'N', 2, A, 2, SDIM, WR, WI, VL, 1, R1.box(1),
        R2.box(1), W, 6, IW, 1, B, INFO);
    chkxer('DGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    dgeesx('N', 'N', dslect, 'N', 1, A, 1, SDIM, WR, WI, VL, 1, R1.box(1),
        R2.box(1), W, 2, IW, 1, B, INFO);
    chkxer('DGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;
  } else if (lsamen(2, C2, 'BD')) {
    // Test DGESVD

    srnamc.SRNAMT = 'DGESVD';
    infoc.INFOT = 1;
    dgesvd('X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgesvd('N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgesvd('O', 'O', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgesvd('N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgesvd('N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgesvd('N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    dgesvd('A', 'N', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    dgesvd('N', 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, INFO);
    chkxer('DGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 8;
    if (infoc.OK.value) {
      print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      print9998(infoc.NOUT);
    }

    // Test DGESDD

    srnamc.SRNAMT = 'DGESDD';
    infoc.INFOT = 1;
    dgesdd('X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
    chkxer('DGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgesdd('N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
    chkxer('DGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgesdd('N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
    chkxer('DGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgesdd('N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, IW, INFO);
    chkxer('DGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    dgesdd('A', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, IW, INFO);
    chkxer('DGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    dgesdd('A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, IW, INFO);
    chkxer('DGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 6;
    if (infoc.OK.value) {
      print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      print9998(infoc.NOUT);
    }

    // Test DGEJSV

    srnamc.SRNAMT = 'DGEJSV';
    infoc.INFOT = 1;
    dgejsv('X', 'U', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgejsv('G', 'X', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgejsv('G', 'U', 'X', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgejsv('G', 'U', 'V', 'X', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgejsv('G', 'U', 'V', 'R', 'X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgejsv('G', 'U', 'V', 'R', 'N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgejsv('G', 'U', 'V', 'R', 'N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    dgejsv('G', 'U', 'V', 'R', 'N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    dgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    dgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 1, VT, 2, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    dgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 2, VT, 1, W, 1, IW,
        INFO);
    chkxer('DGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 11;
    if (infoc.OK.value) {
      print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      print9998(infoc.NOUT);
    }

    // Test DGESVDX

    srnamc.SRNAMT = 'DGESVDX';
    infoc.INFOT = 1;
    dgesvdx('X', 'N', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgesvdx('N', 'X', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgesvdx('N', 'N', 'X', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgesvdx('N', 'N', 'A', -1, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgesvdx('N', 'N', 'A', 0, -1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgesvdx('N', 'N', 'A', 2, 1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    dgesvdx('N', 'N', 'V', 2, 1, A, 2, -ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    dgesvdx('N', 'N', 'V', 2, 1, A, 2, ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    dgesvdx('N', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 0, 1, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    dgesvdx('V', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 1, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    dgesvdx('V', 'N', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 17;
    dgesvdx('N', 'V', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, IW, INFO);
    chkxer('DGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 12;
    if (infoc.OK.value) {
      print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      print9998(infoc.NOUT);
    }

    // Test DGESVDQ

    srnamc.SRNAMT = 'DGESVDQ';
    infoc.INFOT = 1;
    dgesvdq('X', 'P', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgesvdq('A', 'X', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgesvdq('A', 'P', 'X', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgesvdq('A', 'P', 'T', 'X', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgesvdq('A', 'P', 'T', 'A', 'X', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgesvdq('A', 'P', 'T', 'A', 'A', -1, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgesvdq('A', 'P', 'T', 'A', 'A', 0, 1, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 0, S, U, 0, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, -1, VT, 0, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, -1, NS, IW, 1, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 17;
    dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, 1, NS, IW, -5, W,
        1, W, 1, INFO);
    chkxer('DGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 11;
    if (infoc.OK.value) {
      print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      print9998(infoc.NOUT);
    }
  }

  // Print a summary line.

  if (!lsamen(2, C2, 'BD')) {
    if (infoc.OK.value) {
      print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      print9998(infoc.NOUT);
    }
  }
}

void print9999(final Nout nout, final String srnamt, final int nt) {
  nout.println(
      ' ${srnamt.trimRight()} passed the tests of the error exits (${nt.i3} tests done)');
}

void print9998(final Nout nout, [final String s = '']) {
  nout.println(' *** $s failed the tests of the error exits ***');
}
