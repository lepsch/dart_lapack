import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqr2.dart';
import 'package:lapack/src/dgeqr2p.dart';
import 'package:lapack/src/dgeqrf.dart';
import 'package:lapack/src/dgeqrfp.dart';
import 'package:lapack/src/dorg2r.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/dorm2r.dart';
import 'package:lapack/src/dormqr.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrqr(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      W = Array<double>(NMAX),
      X = Array<double>(NMAX);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      AF[I][J] = 1.0 / (I + J);
    }
    B[J] = 0.0;
    W[J] = 0.0;
    X[J] = 0.0;
  }
  infoc.OK.value = true;

  // Error exits for QR factorization

  // DGEQRF

  srnamc.SRNAMT = 'DGEQRF';
  infoc.INFOT = 1;
  dgeqrf(-1, 0, A, 1, B, W, 1, INFO);
  chkxer('DGEQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgeqrf(0, -1, A, 1, B, W, 1, INFO);
  chkxer('DGEQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dgeqrf(2, 1, A, 1, B, W, 1, INFO);
  chkxer('DGEQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dgeqrf(1, 2, A, 1, B, W, 1, INFO);
  chkxer('DGEQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DGEQRFP

  srnamc.SRNAMT = 'DGEQRFP';
  infoc.INFOT = 1;
  dgeqrfp(-1, 0, A, 1, B, W, 1, INFO);
  chkxer('DGEQRFP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgeqrfp(0, -1, A, 1, B, W, 1, INFO);
  chkxer('DGEQRFP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dgeqrfp(2, 1, A, 1, B, W, 1, INFO);
  chkxer('DGEQRFP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dgeqrfp(1, 2, A, 1, B, W, 1, INFO);
  chkxer('DGEQRFP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DGEQR2

  srnamc.SRNAMT = 'DGEQR2';
  infoc.INFOT = 1;
  dgeqr2(-1, 0, A, 1, B, W, INFO);
  chkxer('DGEQR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgeqr2(0, -1, A, 1, B, W, INFO);
  chkxer('DGEQR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dgeqr2(2, 1, A, 1, B, W, INFO);
  chkxer('DGEQR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DGEQR2P

  srnamc.SRNAMT = 'DGEQR2P';
  infoc.INFOT = 1;
  dgeqr2p(-1, 0, A, 1, B, W, INFO);
  chkxer('DGEQR2P', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dgeqr2p(0, -1, A, 1, B, W, INFO);
  chkxer('DGEQR2P', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dgeqr2p(2, 1, A, 1, B, W, INFO);
  chkxer('DGEQR2P', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DORGQR

  srnamc.SRNAMT = 'DORGQR';
  infoc.INFOT = 1;
  dorgqr(-1, 0, 0, A, 1, X, W, 1, INFO);
  chkxer('DORGQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dorgqr(0, -1, 0, A, 1, X, W, 1, INFO);
  chkxer('DORGQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dorgqr(1, 2, 0, A, 1, X, W, 2, INFO);
  chkxer('DORGQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dorgqr(0, 0, -1, A, 1, X, W, 1, INFO);
  chkxer('DORGQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dorgqr(1, 1, 2, A, 1, X, W, 1, INFO);
  chkxer('DORGQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dorgqr(2, 2, 0, A, 1, X, W, 2, INFO);
  chkxer('DORGQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  dorgqr(2, 2, 0, A, 2, X, W, 1, INFO);
  chkxer('DORGQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DORG2R

  srnamc.SRNAMT = 'DORG2R';
  infoc.INFOT = 1;
  dorg2r(-1, 0, 0, A, 1, X, W, INFO);
  chkxer('DORG2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dorg2r(0, -1, 0, A, 1, X, W, INFO);
  chkxer('DORG2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dorg2r(1, 2, 0, A, 1, X, W, INFO);
  chkxer('DORG2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dorg2r(0, 0, -1, A, 1, X, W, INFO);
  chkxer('DORG2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dorg2r(2, 1, 2, A, 2, X, W, INFO);
  chkxer('DORG2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dorg2r(2, 1, 0, A, 1, X, W, INFO);
  chkxer('DORG2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DORMQR

  srnamc.SRNAMT = 'DORMQR';
  infoc.INFOT = 1;
  dormqr('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dormqr('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dormqr('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dormqr('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dormqr('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dormqr('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dormqr('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dormqr('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dormqr('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  dormqr('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  dormqr('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 12;
  dormqr('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO);
  chkxer('DORMQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // DORM2R

  srnamc.SRNAMT = 'DORM2R';
  infoc.INFOT = 1;
  dorm2r('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dorm2r('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dorm2r('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dorm2r('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dorm2r('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dorm2r('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dorm2r('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dorm2r('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dorm2r('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 10;
  dorm2r('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO);
  chkxer('DORM2R', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
