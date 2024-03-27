import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgelqt.dart';
import 'package:lapack/src/zgelqt3.dart';
import 'package:lapack/src/zgemlqt.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrlqt(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<Complex>(NMAX, NMAX),
      T = Matrix<Complex>(NMAX, NMAX),
      W = Array<Complex>(NMAX),
      C = Matrix<Complex>(NMAX, NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = Complex.one / (I + J).toComplex();
      C[I][J] = Complex.one / (I + J).toComplex();
      T[I][J] = Complex.one / (I + J).toComplex();
    }
    W[J] = Complex.zero;
  }
  infoc.OK.value = true;

  // Error exits for LQT factorization

  // ZGELQT

  srnamc.SRNAMT = 'ZGELQT';
  infoc.INFOT = 1;
  zgelqt(-1, 0, 1, A, 1, T, 1, W, INFO);
  chkxer('ZGELQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgelqt(0, -1, 1, A, 1, T, 1, W, INFO);
  chkxer('ZGELQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zgelqt(0, 0, 0, A, 1, T, 1, W, INFO);
  chkxer('ZGELQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zgelqt(2, 1, 1, A, 1, T, 1, W, INFO);
  chkxer('ZGELQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zgelqt(2, 2, 2, A, 2, T, 1, W, INFO);
  chkxer('ZGELQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGELQT3

  srnamc.SRNAMT = 'ZGELQT3';
  infoc.INFOT = 1;
  zgelqt3(-1, 0, A, 1, T, 1, INFO);
  chkxer('ZGELQT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgelqt3(0, -1, A, 1, T, 1, INFO);
  chkxer('ZGELQT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgelqt3(2, 2, A, 1, T, 1, INFO);
  chkxer('ZGELQT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  zgelqt3(2, 2, A, 2, T, 1, INFO);
  chkxer('ZGELQT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZGEMLQT

  srnamc.SRNAMT = 'ZGEMLQT';
  infoc.INFOT = 1;
  zgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  zgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  zgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  zgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO);
  chkxer('ZGEMLQT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
