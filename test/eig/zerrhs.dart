import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgebak.dart';
import 'package:lapack/src/zgebal.dart';
import 'package:lapack/src/zgehd2.dart';
import 'package:lapack/src/zgehrd.dart';
import 'package:lapack/src/zhsein.dart';
import 'package:lapack/src/zhseqr.dart';
import 'package:lapack/src/ztrevc.dart';
import 'package:lapack/src/ztrevc3.dart';
import 'package:lapack/src/zunghr.dart';
import 'package:lapack/src/zunmhr.dart';

import 'chkxer.dart';
import 'common.dart';

void zerrhs(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 3, LW = NMAX * NMAX;
  String C2;
  int I, J, NT;
  final SEL = Array<bool>(NMAX);
  final IFAILL = Array<int>(NMAX), IFAILR = Array<int>(NMAX);
  final RW = Array<double>(NMAX), S = Array<double>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      C = Matrix<Complex>(NMAX, NMAX),
      VL = Matrix<Complex>(NMAX, NMAX),
      VR = Matrix<Complex>(NMAX, NMAX);
  final TAU = Array<Complex>(NMAX),
      W = Array<Complex>(LW),
      X = Array<Complex>(NMAX);
  final IHI = Box(0), ILO = Box(0), INFO = Box(0), M = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();
  C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (J = 1; J <= NMAX; J++) {
    for (I = 1; I <= NMAX; I++) {
      A[I][J] = (1.0 / (I + J)).toComplex();
    }
    SEL[J] = true;
  }
  infoc.OK.value = true;
  NT = 0;

  // Test error exits of the nonsymmetric eigenvalue routines.

  if (lsamen(2, C2, 'HS')) {
    // ZGEBAL

    srnamc.SRNAMT = 'ZGEBAL';
    infoc.INFOT = 1;
    zgebal('/', 0, A, 1, ILO, IHI, S, INFO);
    chkxer('ZGEBAL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgebal('N', -1, A, 1, ILO, IHI, S, INFO);
    chkxer('ZGEBAL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgebal('N', 2, A, 1, ILO, IHI, S, INFO);
    chkxer('ZGEBAL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 3;

    // ZGEBAK

    srnamc.SRNAMT = 'ZGEBAK';
    infoc.INFOT = 1;
    zgebak('/', 'R', 0, 1, 0, S, 0, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgebak('N', '/', 0, 1, 0, S, 0, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgebak('N', 'R', -1, 1, 0, S, 0, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgebak('N', 'R', 0, 0, 0, S, 0, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgebak('N', 'R', 0, 2, 0, S, 0, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgebak('N', 'R', 2, 2, 1, S, 0, A, 2, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgebak('N', 'R', 0, 1, 1, S, 0, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgebak('N', 'R', 0, 1, 0, S, -1, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgebak('N', 'R', 2, 1, 2, S, 0, A, 1, INFO);
    chkxer('ZGEBAK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 9;

    // ZGEHRD

    srnamc.SRNAMT = 'ZGEHRD';
    infoc.INFOT = 1;
    zgehrd(-1, 1, 1, A, 1, TAU, W, 1, INFO);
    chkxer('ZGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgehrd(0, 0, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgehrd(0, 2, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgehrd(1, 1, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgehrd(0, 1, 1, A, 1, TAU, W, 1, INFO);
    chkxer('ZGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgehrd(2, 1, 1, A, 1, TAU, W, 2, INFO);
    chkxer('ZGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zgehrd(2, 1, 2, A, 2, TAU, W, 1, INFO);
    chkxer('ZGEHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;

    // ZGEHD2

    srnamc.SRNAMT = 'ZGEHD2';
    infoc.INFOT = 1;
    zgehd2(-1, 1, 1, A, 1, TAU, W, INFO);
    chkxer('ZGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgehd2(0, 0, 0, A, 1, TAU, W, INFO);
    chkxer('ZGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgehd2(0, 2, 0, A, 1, TAU, W, INFO);
    chkxer('ZGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgehd2(1, 1, 0, A, 1, TAU, W, INFO);
    chkxer('ZGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgehd2(0, 1, 1, A, 1, TAU, W, INFO);
    chkxer('ZGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgehd2(2, 1, 1, A, 1, TAU, W, INFO);
    chkxer('ZGEHD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 6;

    // ZUNGHR

    srnamc.SRNAMT = 'ZUNGHR';
    infoc.INFOT = 1;
    zunghr(-1, 1, 1, A, 1, TAU, W, 1, INFO);
    chkxer('ZUNGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zunghr(0, 0, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZUNGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zunghr(0, 2, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZUNGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zunghr(1, 1, 0, A, 1, TAU, W, 1, INFO);
    chkxer('ZUNGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zunghr(0, 1, 1, A, 1, TAU, W, 1, INFO);
    chkxer('ZUNGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zunghr(2, 1, 1, A, 1, TAU, W, 1, INFO);
    chkxer('ZUNGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zunghr(3, 1, 3, A, 3, TAU, W, 1, INFO);
    chkxer('ZUNGHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;

    // ZUNMHR

    srnamc.SRNAMT = 'ZUNMHR';
    infoc.INFOT = 1;
    zunmhr('/', 'N', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zunmhr('L', '/', 0, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zunmhr('L', 'N', -1, 0, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zunmhr('L', 'N', 0, -1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zunmhr('L', 'N', 0, 0, 0, 0, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zunmhr('L', 'N', 0, 0, 2, 0, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zunmhr('L', 'N', 1, 2, 2, 1, A, 1, TAU, C, 1, W, 2, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zunmhr('R', 'N', 2, 1, 2, 1, A, 1, TAU, C, 2, W, 2, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zunmhr('L', 'N', 1, 1, 1, 0, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zunmhr('L', 'N', 0, 1, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zunmhr('R', 'N', 1, 0, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zunmhr('L', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zunmhr('R', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zunmhr('L', 'N', 2, 1, 1, 1, A, 2, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zunmhr('L', 'N', 1, 2, 1, 1, A, 1, TAU, C, 1, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zunmhr('R', 'N', 2, 1, 1, 1, A, 1, TAU, C, 2, W, 1, INFO);
    chkxer('ZUNMHR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 16;

    // ZHSEQR

    srnamc.SRNAMT = 'ZHSEQR';
    infoc.INFOT = 1;
    zhseqr('/', 'N', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhseqr('E', '/', 0, 1, 0, A, 1, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhseqr('E', 'N', -1, 1, 0, A, 1, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zhseqr('E', 'N', 0, 0, 0, A, 1, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zhseqr('E', 'N', 0, 2, 0, A, 1, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhseqr('E', 'N', 1, 1, 0, A, 1, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhseqr('E', 'N', 1, 1, 2, A, 1, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zhseqr('E', 'N', 2, 1, 2, A, 1, X, C, 2, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhseqr('E', 'V', 2, 1, 2, A, 2, X, C, 1, W, 1, INFO);
    chkxer('ZHSEQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 9;

    // ZHSEIN

    srnamc.SRNAMT = 'ZHSEIN';
    infoc.INFOT = 1;
    zhsein('/', 'N', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhsein('R', '/', 'N', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhsein('R', 'N', '/', SEL, 0, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhsein('R', 'N', 'N', SEL, -1, A, 1, X, VL, 1, VR, 1, 0, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zhsein('R', 'N', 'N', SEL, 2, A, 1, X, VL, 1, VR, 2, 4, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhsein('L', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zhsein('R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 1, 4, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zhsein('R', 'N', 'N', SEL, 2, A, 2, X, VL, 1, VR, 2, 1, M, W, RW, IFAILL,
        IFAILR, INFO);
    chkxer('ZHSEIN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 8;

    // ZTREVC

    srnamc.SRNAMT = 'ZTREVC';
    infoc.INFOT = 1;
    ztrevc('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO);
    chkxer('ZTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztrevc('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO);
    chkxer('ZTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztrevc('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, RW, INFO);
    chkxer('ZTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztrevc('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, RW, INFO);
    chkxer('ZTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztrevc('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO);
    chkxer('ZTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztrevc('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, RW, INFO);
    chkxer('ZTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    ztrevc('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, RW, INFO);
    chkxer('ZTREVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;

    // ZTREVC3

    srnamc.SRNAMT = 'ZTREVC3';
    infoc.INFOT = 1;
    ztrevc3('/', 'A', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztrevc3('L', '/', SEL, 0, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztrevc3('L', 'A', SEL, -1, A, 1, VL, 1, VR, 1, 0, M, W, LW, RW, 1, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztrevc3('L', 'A', SEL, 2, A, 1, VL, 2, VR, 1, 4, M, W, LW, RW, 2, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztrevc3('L', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztrevc3('R', 'A', SEL, 2, A, 2, VL, 1, VR, 1, 4, M, W, LW, RW, 2, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    ztrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 1, M, W, LW, RW, 2, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    ztrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, 2, RW, 2, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    ztrevc3('L', 'A', SEL, 2, A, 2, VL, 2, VR, 1, 2, M, W, LW, RW, 1, INFO);
    chkxer('ZTREVC3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 9;
  }

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT.println(
        ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)');
  } else {
    infoc.NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
