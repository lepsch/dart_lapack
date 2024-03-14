import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zcgesv.dart';

import 'chkxer.dart';
import 'common.dart';

void zerrab(Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4;
  final IP = Array<int>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      AF = Matrix<Complex>(NMAX, NMAX),
      B = Array<Complex>(NMAX),
      C = Array<Complex>(NMAX),
      R = Array<Complex>(NMAX),
      R1 = Array<Complex>(NMAX),
      R2 = Array<Complex>(NMAX),
      W = Array<Complex>(2 * NMAX),
      X = Array<Complex>(NMAX);
  final WORK = Array<Complex>(1);
  final SWORK = Array<Complex>(1);
  final RWORK = Array<double>(1);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = (1.0 / (I + J)).toComplex();
      AF[I][J] = (1.0 / (I + J)).toComplex();
    }
    B[J] = Complex.zero;
    R1[J] = Complex.zero;
    R2[J] = Complex.zero;
    W[J] = Complex.zero;
    X[J] = Complex.zero;
    C[J] = Complex.zero;
    R[J] = Complex.zero;
    IP[J] = J;
  }
  infoc.OK.value = true;

  final ITER = Box(0);
  srnamc.SRNAMT = 'ZCGESV';
  infoc.INFOT = 1;
  zcgesv(-1, 0, A, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zcgesv(0, -1, A, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zcgesv(2, 1, A, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zcgesv(2, 1, A, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 9;
  zcgesv(2, 1, A, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1, WORK.asMatrix(),
      SWORK, RWORK, ITER, INFO);
  chkxer('ZCGESV', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  if (infoc.OK.value) {
    NOUT.println(' ZCGESV drivers passed the tests of the error exits');
  } else {
    NOUT.println(' *** ZCGESV drivers failed the tests of the error exits ***');
  }
}
