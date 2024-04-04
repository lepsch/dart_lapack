import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgees.dart';
import 'package:lapack/src/zgeesx.dart';
import 'package:lapack/src/zgeev.dart';
import 'package:lapack/src/zgeevx.dart';
import 'package:lapack/src/zgejsv.dart';
import 'package:lapack/src/zgesdd.dart';
import 'package:lapack/src/zgesvd.dart';
import 'package:lapack/src/zgesvdq.dart';
import 'package:lapack/src/zgesvdx.dart';

import 'chkxer.dart';
import 'common.dart';
import 'zslect.dart';

void zerred(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4, LW = 5 * NMAX;
  const ONE = 1.0, ZERO = 0.0;
  String C2;
  int I, J, NT;
  final B = Array<bool>(NMAX);
  final IW = Array<int>(4 * NMAX);
  final R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      RW = Array<double>(LW),
      S = Array<double>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      U = Matrix<Complex>(NMAX, NMAX),
      VL = Matrix<Complex>(NMAX, NMAX),
      VR = Matrix<Complex>(NMAX, NMAX),
      VT = Matrix<Complex>(NMAX, NMAX);
  final W = Array<Complex>(10 * NMAX), X = Array<Complex>(NMAX);
  final INFO = Box(0), SDIM = Box(0), IHI = Box(0), ILO = Box(0), NS = Box(0);
  final ABNRM = Box(0.0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();
  C2 = PATH.substring(1, 3);

  // Initialize A

  for (J = 1; J <= NMAX; J++) {
    for (I = 1; I <= NMAX; I++) {
      A[I][J] = Complex.zero;
    }
  }
  for (I = 1; I <= NMAX; I++) {
    A[I][I] = Complex.one;
  }
  infoc.OK.value = true;
  NT = 0;

  if (lsamen(2, C2, 'EV')) {
    // Test ZGEEV

    srnamc.SRNAMT = 'ZGEEV';
    infoc.INFOT = 1;
    zgeev('X', 'N', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO);
    chkxer('ZGEEV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgeev('N', 'X', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO);
    chkxer('ZGEEV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgeev('N', 'N', -1, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO);
    chkxer('ZGEEV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgeev('N', 'N', 2, A, 1, X, VL, 1, VR, 1, W, 4, RW, INFO);
    chkxer('ZGEEV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgeev('V', 'N', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW, INFO);
    chkxer('ZGEEV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgeev('N', 'V', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW, INFO);
    chkxer('ZGEEV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgeev('V', 'V', 1, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO);
    chkxer('ZGEEV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;
  } else if (lsamen(2, C2, 'ES')) {
    // Test ZGEES

    srnamc.SRNAMT = 'ZGEES';
    infoc.INFOT = 1;
    zgees('X', 'N', zslect, 0, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO);
    chkxer('ZGEES', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgees('N', 'X', zslect, 0, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO);
    chkxer('ZGEES', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgees('N', 'S', zslect, -1, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO);
    chkxer('ZGEES', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgees('N', 'S', zslect, 2, A, 1, SDIM, X, VL, 1, W, 4, RW, B, INFO);
    chkxer('ZGEES', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgees('V', 'S', zslect, 2, A, 2, SDIM, X, VL, 1, W, 4, RW, B, INFO);
    chkxer('ZGEES', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgees('N', 'S', zslect, 1, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO);
    chkxer('ZGEES', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 6;
  } else if (lsamen(2, C2, 'VX')) {
    // Test ZGEEVX

    srnamc.SRNAMT = 'ZGEEVX';
    infoc.INFOT = 1;
    zgeevx('X', 'N', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 1, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgeevx('N', 'X', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 1, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgeevx('N', 'N', 'X', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 1, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgeevx('N', 'N', 'N', 'X', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 1, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgeevx('N', 'N', 'N', 'N', -1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM,
        R1, R2, W, 1, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgeevx('N', 'N', 'N', 'N', 2, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 4, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgeevx('N', 'V', 'N', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 4, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgeevx('N', 'N', 'V', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 4, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    zgeevx('N', 'N', 'N', 'N', 1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 1, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    zgeevx('N', 'N', 'V', 'V', 1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1,
        R2, W, 2, RW, INFO);
    chkxer('ZGEEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 10;
  } else if (lsamen(2, C2, 'SX')) {
    // Test ZGEESX

    srnamc.SRNAMT = 'ZGEESX';
    infoc.INFOT = 1;
    zgeesx('X', 'N', zslect, 'N', 0, A, 1, SDIM, X, VL, 1, R1(1), R2(1), W, 1,
        RW, B, INFO);
    chkxer('ZGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgeesx('N', 'X', zslect, 'N', 0, A, 1, SDIM, X, VL, 1, R1(1), R2(1), W, 1,
        RW, B, INFO);
    chkxer('ZGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgeesx('N', 'N', zslect, 'X', 0, A, 1, SDIM, X, VL, 1, R1(1), R2(1), W, 1,
        RW, B, INFO);
    chkxer('ZGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgeesx('N', 'N', zslect, 'N', -1, A, 1, SDIM, X, VL, 1, R1(1), R2(1), W, 1,
        RW, B, INFO);
    chkxer('ZGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgeesx('N', 'N', zslect, 'N', 2, A, 1, SDIM, X, VL, 1, R1(1), R2(1), W, 4,
        RW, B, INFO);
    chkxer('ZGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zgeesx('V', 'N', zslect, 'N', 2, A, 2, SDIM, X, VL, 1, R1(1), R2(1), W, 4,
        RW, B, INFO);
    chkxer('ZGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zgeesx('N', 'N', zslect, 'N', 1, A, 1, SDIM, X, VL, 1, R1(1), R2(1), W, 1,
        RW, B, INFO);
    chkxer('ZGEESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;
  } else if (lsamen(2, C2, 'BD')) {
    // Test ZGESVD

    srnamc.SRNAMT = 'ZGESVD';
    infoc.INFOT = 1;
    zgesvd('X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgesvd('N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgesvd('O', 'O', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgesvd('N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgesvd('N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgesvd('N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgesvd('A', 'N', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zgesvd('N', 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW, INFO);
    chkxer('ZGESVD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 8;
    if (infoc.OK.value) {
      _print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      _print9998(infoc.NOUT, srnamc.SRNAMT);
    }

    // Test ZGESDD

    srnamc.SRNAMT = 'ZGESDD';
    infoc.INFOT = 1;
    zgesdd('X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO);
    chkxer('ZGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgesdd('N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO);
    chkxer('ZGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgesdd('N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO);
    chkxer('ZGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgesdd('N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW, IW, INFO);
    chkxer('ZGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgesdd('A', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW, IW, INFO);
    chkxer('ZGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgesdd('A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW, IW, INFO);
    chkxer('ZGESDD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT -= 2;
    if (infoc.OK.value) {
      _print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      _print9998(infoc.NOUT, srnamc.SRNAMT);
    }

    // Test ZGEJSV

    srnamc.SRNAMT = 'ZGEJSV';
    infoc.INFOT = 1;
    zgejsv('X', 'U', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgejsv('G', 'X', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgejsv('G', 'U', 'X', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgejsv('G', 'U', 'V', 'X', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgejsv('G', 'U', 'V', 'R', 'X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgejsv('G', 'U', 'V', 'R', 'N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgejsv('G', 'U', 'V', 'R', 'N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgejsv('G', 'U', 'V', 'R', 'N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 1, VT, 2, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 2, VT, 1, W, 1, RW,
        1, IW, INFO);
    chkxer('ZGEJSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 11;
    if (infoc.OK.value) {
      _print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      _print9998(infoc.NOUT, srnamc.SRNAMT);
    }

    // Test ZGESVDX

    srnamc.SRNAMT = 'ZGESVDX';
    infoc.INFOT = 1;
    zgesvdx('X', 'N', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgesvdx('N', 'X', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgesvdx('N', 'N', 'X', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgesvdx('N', 'N', 'A', -1, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgesvdx('N', 'N', 'A', 0, -1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgesvdx('N', 'N', 'A', 2, 1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgesvdx('N', 'N', 'V', 2, 1, A, 2, -ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgesvdx('N', 'N', 'V', 2, 1, A, 2, ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgesvdx('N', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 0, 1, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zgesvdx('V', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 1, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zgesvdx('V', 'N', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 17;
    zgesvdx('N', 'V', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W,
        1, RW, IW, INFO);
    chkxer('ZGESVDX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 12;
    if (infoc.OK.value) {
      _print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      _print9998(infoc.NOUT, srnamc.SRNAMT);
    }

    // Test ZGESVDQ

    srnamc.SRNAMT = 'ZGESVDQ';
    infoc.INFOT = 1;
    zgesvdq('X', 'P', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgesvdq('A', 'X', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgesvdq('A', 'P', 'X', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgesvdq('A', 'P', 'T', 'X', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgesvdq('A', 'P', 'T', 'A', 'X', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zgesvdq('A', 'P', 'T', 'A', 'A', -1, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgesvdq('A', 'P', 'T', 'A', 'A', 0, 1, A, 1, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 0, S, U, 0, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, -1, VT, 0, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, -1, NS, IW, 1, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 17;
    zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, 1, NS, IW, -5, W,
        1, RW, 1, INFO);
    chkxer('ZGESVDQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = 11;
    if (infoc.OK.value) {
      _print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      _print9998(infoc.NOUT, srnamc.SRNAMT);
    }
  }

  // Print a summary line.

  if (!lsamen(2, C2, 'BD')) {
    if (infoc.OK.value) {
      _print9999(infoc.NOUT, srnamc.SRNAMT, NT);
    } else {
      _print9998(infoc.NOUT, srnamc.SRNAMT);
    }
  }
}

void _print9999(Nout nout, String s, int t) {
  nout.println(' $s passed the tests of the error exits (${t.i3} tests done)');
}

void _print9998(Nout nout, String s) {
  nout.println(' *** $s failed the tests of the error exits ***');
}
