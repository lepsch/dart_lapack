import 'package:lapack/src/box.dart';
import 'package:lapack/src/dpstf2.dart';
import 'package:lapack/src/dpstrf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrps(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4;
  final A = Matrix<double>(NMAX, NMAX), WORK = Array<double>(2 * NMAX);
  final PIV = Array<int>(NMAX);
  final INFO = Box(0), RANK = Box(0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
    }
    PIV[J] = J;
    WORK[J] = 0.0;
    WORK[NMAX + J] = 0.0;
  }
  OK.value = true;

  // Test error exits of the routines that use the Cholesky
  // decomposition of a symmetric positive semidefinite matrix.

  test('DPSTRF', () {
    srnamc.SRNAMT = 'DPSTRF';
    infoc.INFOT = 1;
    dpstrf('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO);
    chkxer('DPSTRF', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dpstrf('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO);
    chkxer('DPSTRF', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dpstrf('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO);
    chkxer('DPSTRF', infoc.INFOT, NOUT, LERR, OK, test);
  });

  test('DPSTF2', () {
    srnamc.SRNAMT = 'DPSTF2';
    infoc.INFOT = 1;
    dpstf2('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO);
    chkxer('DPSTF2', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 2;
    dpstf2('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO);
    chkxer('DPSTF2', infoc.INFOT, NOUT, LERR, OK, test);
    infoc.INFOT = 4;
    dpstf2('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO);
    chkxer('DPSTF2', infoc.INFOT, NOUT, LERR, OK, test);
  });

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
