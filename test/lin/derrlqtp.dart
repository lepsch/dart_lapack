import 'package:lapack/src/box.dart';
import 'package:lapack/src/dtplqt2.dart';
import 'package:lapack/src/dtpmlqt.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';
import 'dtplqt.dart';

void derrlqtp(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      T = Matrix<double>(NMAX, NMAX),
      W = Array<double>(NMAX),
      B = Matrix<double>(NMAX, NMAX),
      C = Matrix<double>(NMAX, NMAX);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      C[I][J] = 1.0 / (I + J);
      T[I][J] = 1.0 / (I + J);
    }
    W[J] = 0.0;
  }
  infoc.OK.value = true;

  // Error exits for TPLQT factorization

  // DTPLQT

  srnamc.SRNAMT = 'DTPLQT';
  infoc.INFOT = 1;
  dtplqt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtplqt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtplqt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtplqt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtplqt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtplqt(1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dtplqt(2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  dtplqt(2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  dtplqt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO);
  chkxer('DTPLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DTPLQT2

  srnamc.SRNAMT = 'DTPLQT2';
  infoc.INFOT = 1;
  dtplqt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO);
  chkxer('DTPLQT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtplqt2(0, -1, 0, A, 1, B, 1, T, 1, INFO);
  chkxer('DTPLQT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtplqt2(0, 0, -1, A, 1, B, 1, T, 1, INFO);
  chkxer('DTPLQT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dtplqt2(2, 2, 0, A, 1, B, 2, T, 2, INFO);
  chkxer('DTPLQT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dtplqt2(2, 2, 0, A, 2, B, 1, T, 2, INFO);
  chkxer('DTPLQT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  dtplqt2(2, 2, 0, A, 2, B, 2, T, 1, INFO);
  chkxer('DTPLQT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DTPMLQT

  srnamc.SRNAMT = 'DTPMLQT';
  infoc.INFOT = 1;
  dtpmlqt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtpmlqt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtpmlqt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtpmlqt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dtpmlqt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  infoc.INFOT = 6;
  dtpmlqt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dtpmlqt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  dtpmlqt('R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  dtpmlqt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 13;
  dtpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 15;
  dtpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO);
  chkxer('DTPMLQT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
