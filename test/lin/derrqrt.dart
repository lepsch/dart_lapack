import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgemqrt.dart';
import 'package:lapack/src/dgeqrt.dart';
import 'package:lapack/src/dgeqrt2.dart';
import 'package:lapack/src/dgeqrt3.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrqrt(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  int I, J;
  final A = Matrix<double>(NMAX, NMAX),
      T = Matrix<double>(NMAX, NMAX),
      W = Array<double>(NMAX),
      C = Matrix<double>(NMAX, NMAX);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();

  // Set the variables to innocuous values.

  for (J = 1; J <= NMAX; J++) {
    for (I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      C[I][J] = 1.0 / (I + J);
      T[I][J] = 1.0 / (I + J);
    }
    W[J] = 0.0;
  }
  infoc.OK.value = true;

  // Error exits for QRT factorization

  // DGEQRT

  srnamc.SRNAMT = 'DGEQRT';
  infoc.INFOT = 1;
  dgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO);
  chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgeqrt(0, -1, 1, A, 1, T, 1, W, INFO);
  chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dgeqrt(0, 0, 0, A, 1, T, 1, W, INFO);
  chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dgeqrt(2, 1, 1, A, 1, T, 1, W, INFO);
  chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dgeqrt(2, 2, 2, A, 2, T, 1, W, INFO);
  chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DGEQRT2

  srnamc.SRNAMT = 'DGEQRT2';
  infoc.INFOT = 1;
  dgeqrt2(-1, 0, A, 1, T, 1, INFO);
  chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgeqrt2(0, -1, A, 1, T, 1, INFO);
  chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dgeqrt2(2, 1, A, 1, T, 1, INFO);
  chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dgeqrt2(2, 2, A, 2, T, 1, INFO);
  chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DGEQRT3

  srnamc.SRNAMT = 'DGEQRT3';
  infoc.INFOT = 1;
  dgeqrt3(-1, 0, A, 1, T, 1, INFO);
  chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgeqrt3(0, -1, A, 1, T, 1, INFO);
  chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dgeqrt3(2, 1, A, 1, T, 1, INFO);
  chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dgeqrt3(2, 2, A, 2, T, 1, INFO);
  chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DGEMQRT

  srnamc.SRNAMT = 'DGEMQRT';
  infoc.INFOT = 1;
  dgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  dgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  dgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  dgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  dgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO);
  chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
