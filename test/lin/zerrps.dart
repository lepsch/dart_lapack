import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zpstf2.dart';
import 'package:lapack/src/zpstrf.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrps(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4;
  final A = Matrix<Complex>(NMAX, NMAX);
  final RWORK = Array<double>(2 * NMAX);
  final PIV = Array<int>(NMAX);
  final INFO = Box(0);

  final NOUT = NUNIT;
  NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = Complex.one / (I + J).toComplex();
    }
    PIV[J] = J;
    RWORK[J] = 0.0;
    RWORK[NMAX + J] = 0.0;
  }
  infoc.OK.value = true;

  // Test error exits of the routines that use the Cholesky
  // decomposition of an Hermitian positive semidefinite matrix.

  // ZPSTRF

  srnamc.SRNAMT = 'ZPSTRF';
  infoc.INFOT = 1;
  final RANK = Box(0);
  zpstrf('/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO);
  chkxer('ZPSTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zpstrf('U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO);
  chkxer('ZPSTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zpstrf('U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO);
  chkxer('ZPSTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // ZPSTF2

  srnamc.SRNAMT = 'ZPSTF2';
  infoc.INFOT = 1;
  zpstf2('/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO);
  chkxer('ZPSTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zpstf2('U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO);
  chkxer('ZPSTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zpstf2('U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO);
  chkxer('ZPSTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, NOUT);
}
