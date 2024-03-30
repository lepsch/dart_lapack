import 'package:lapack/lapack.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrpo(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4;
  final IW = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      W = Array<double>(3 * NMAX),
      X = Array<double>(NMAX);
  final INFO = Box(0);
  final ANRM = Box(0.0), RCOND = Box(0.0);

  final NOUT = infoc.NOUT = NUNIT;
  final OK = infoc.OK;
  final LERR = infoc.LERR;
  NOUT.println();
  final C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      AF[I][J] = 1.0 / (I + J);
    }
    B[J] = 0.0;
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = 0.0;
    X[J] = 0.0;
    IW[J] = J;
  }
  infoc.OK.value = true;

  if (lsamen(2, C2, 'PO')) {
    // Test error exits of the routines that use the Cholesky
    // decomposition of a symmetric positive definite matrix.

    test('DPOTRF', () {
      srnamc.SRNAMT = 'DPOTRF';
      infoc.INFOT = 1;
      dpotrf('/', 0, A, 1, INFO);
      chkxer('DPOTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpotrf('U', -1, A, 1, INFO);
      chkxer('DPOTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpotrf('U', 2, A, 1, INFO);
      chkxer('DPOTRF', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPOTF2', () {
      srnamc.SRNAMT = 'DPOTF2';
      infoc.INFOT = 1;
      dpotf2('/', 0, A, 1, INFO);
      chkxer('DPOTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpotf2('U', -1, A, 1, INFO);
      chkxer('DPOTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpotf2('U', 2, A, 1, INFO);
      chkxer('DPOTF2', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPOTRI', () {
      srnamc.SRNAMT = 'DPOTRI';
      infoc.INFOT = 1;
      dpotri('/', 0, A, 1, INFO);
      chkxer('DPOTRI', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpotri('U', -1, A, 1, INFO);
      chkxer('DPOTRI', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpotri('U', 2, A, 1, INFO);
      chkxer('DPOTRI', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPOTRS', () {
      srnamc.SRNAMT = 'DPOTRS';
      infoc.INFOT = 1;
      dpotrs('/', 0, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPOTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpotrs('U', -1, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPOTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpotrs('U', 0, -1, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPOTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dpotrs('U', 2, 1, A, 1, B.asMatrix(), 2, INFO);
      chkxer('DPOTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dpotrs('U', 2, 1, A, 2, B.asMatrix(), 1, INFO);
      chkxer('DPOTRS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPORFS', () {
      srnamc.SRNAMT = 'DPORFS';
      infoc.INFOT = 1;
      dporfs('/', 0, 0, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
          W, IW, INFO);
      chkxer('DPORFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dporfs('U', -1, 0, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
          W, IW, INFO);
      chkxer('DPORFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dporfs('U', 0, -1, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
          W, IW, INFO);
      chkxer('DPORFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dporfs('U', 2, 1, A, 1, AF, 2, B.asMatrix(), 2, X.asMatrix(), 2, R1, R2,
          W, IW, INFO);
      chkxer('DPORFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dporfs('U', 2, 1, A, 2, AF, 1, B.asMatrix(), 2, X.asMatrix(), 2, R1, R2,
          W, IW, INFO);
      chkxer('DPORFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dporfs('U', 2, 1, A, 2, AF, 2, B.asMatrix(), 1, X.asMatrix(), 2, R1, R2,
          W, IW, INFO);
      chkxer('DPORFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 11;
      dporfs('U', 2, 1, A, 2, AF, 2, B.asMatrix(), 2, X.asMatrix(), 1, R1, R2,
          W, IW, INFO);
      chkxer('DPORFS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPOCON', () {
      srnamc.SRNAMT = 'DPOCON';
      infoc.INFOT = 1;
      dpocon('/', 0, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPOCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpocon('U', -1, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPOCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpocon('U', 2, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPOCON', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPOEQU', () {
      srnamc.SRNAMT = 'DPOEQU';
      infoc.INFOT = 1;
      dpoequ(-1, A, 1, R1, RCOND, ANRM, INFO);
      chkxer('DPOEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpoequ(2, A, 1, R1, RCOND, ANRM, INFO);
      chkxer('DPOEQU', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'PP')) {
    // Test error exits of the routines that use the Cholesky
    // decomposition of a symmetric positive definite packed matrix.

    test('DPPTRF', () {
      srnamc.SRNAMT = 'DPPTRF';
      infoc.INFOT = 1;
      dpptrf('/', 0, A.asArray(), INFO);
      chkxer('DPPTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpptrf('U', -1, A.asArray(), INFO);
      chkxer('DPPTRF', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPPTRI', () {
      srnamc.SRNAMT = 'DPPTRI';
      infoc.INFOT = 1;
      dpptri('/', 0, A.asArray(), INFO);
      chkxer('DPPTRI', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpptri('U', -1, A.asArray(), INFO);
      chkxer('DPPTRI', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPPTRS', () {
      srnamc.SRNAMT = 'DPPTRS';
      infoc.INFOT = 1;
      dpptrs('/', 0, 0, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpptrs('U', -1, 0, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpptrs('U', 0, -1, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dpptrs('U', 2, 1, A.asArray(), B.asMatrix(), 1, INFO);
      chkxer('DPPTRS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPPRFS', () {
      srnamc.SRNAMT = 'DPPRFS';
      infoc.INFOT = 1;
      dpprfs('/', 0, 0, A.asArray(), AF.asArray(), B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DPPRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpprfs('U', -1, 0, A.asArray(), AF.asArray(), B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DPPRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpprfs('U', 0, -1, A.asArray(), AF.asArray(), B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DPPRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 7;
      dpprfs('U', 2, 1, A.asArray(), AF.asArray(), B.asMatrix(), 1,
          X.asMatrix(), 2, R1, R2, W, IW, INFO);
      chkxer('DPPRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 9;
      dpprfs('U', 2, 1, A.asArray(), AF.asArray(), B.asMatrix(), 2,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DPPRFS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPPCON', () {
      srnamc.SRNAMT = 'DPPCON';
      infoc.INFOT = 1;
      dppcon('/', 0, A.asArray(), ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPPCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dppcon('U', -1, A.asArray(), ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPPCON', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPPEQU', () {
      srnamc.SRNAMT = 'DPPEQU';
      infoc.INFOT = 1;
      dppequ('/', 0, A.asArray(), R1, RCOND, ANRM, INFO);
      chkxer('DPPEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dppequ('U', -1, A.asArray(), R1, RCOND, ANRM, INFO);
      chkxer('DPPEQU', infoc.INFOT, NOUT, LERR, OK, test);
    });
  } else if (lsamen(2, C2, 'PB')) {
    // Test error exits of the routines that use the Cholesky
    // decomposition of a symmetric positive definite band matrix.

    test('DPBTRF', () {
      srnamc.SRNAMT = 'DPBTRF';
      infoc.INFOT = 1;
      dpbtrf('/', 0, 0, A, 1, INFO);
      chkxer('DPBTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbtrf('U', -1, 0, A, 1, INFO);
      chkxer('DPBTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbtrf('U', 1, -1, A, 1, INFO);
      chkxer('DPBTRF', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dpbtrf('U', 2, 1, A, 1, INFO);
      chkxer('DPBTRF', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPBTF2', () {
      srnamc.SRNAMT = 'DPBTF2';
      infoc.INFOT = 1;
      dpbtf2('/', 0, 0, A, 1, INFO);
      chkxer('DPBTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbtf2('U', -1, 0, A, 1, INFO);
      chkxer('DPBTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbtf2('U', 1, -1, A, 1, INFO);
      chkxer('DPBTF2', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dpbtf2('U', 2, 1, A, 1, INFO);
      chkxer('DPBTF2', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPBTRS', () {
      srnamc.SRNAMT = 'DPBTRS';
      infoc.INFOT = 1;
      dpbtrs('/', 0, 0, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbtrs('U', -1, 0, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbtrs('U', 1, -1, 0, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpbtrs('U', 0, 0, -1, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dpbtrs('U', 2, 1, 1, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBTRS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dpbtrs('U', 2, 0, 1, A, 1, B.asMatrix(), 1, INFO);
      chkxer('DPBTRS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPBRFS', () {
      srnamc.SRNAMT = 'DPBRFS';
      infoc.INFOT = 1;
      dpbrfs('/', 0, 0, 0, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbrfs('U', -1, 0, 0, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbrfs('U', 1, -1, 0, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 4;
      dpbrfs('U', 0, 0, -1, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 6;
      dpbrfs('U', 2, 1, 1, A, 1, AF, 2, B.asMatrix(), 2, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 8;
      dpbrfs('U', 2, 1, 1, A, 2, AF, 1, B.asMatrix(), 2, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 10;
      dpbrfs('U', 2, 0, 1, A, 1, AF, 1, B.asMatrix(), 1, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 12;
      dpbrfs('U', 2, 0, 1, A, 1, AF, 1, B.asMatrix(), 2, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DPBRFS', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPBCON', () {
      srnamc.SRNAMT = 'DPBCON';
      infoc.INFOT = 1;
      dpbcon('/', 0, 0, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPBCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbcon('U', -1, 0, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPBCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbcon('U', 1, -1, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPBCON', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dpbcon('U', 2, 1, A, 1, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DPBCON', infoc.INFOT, NOUT, LERR, OK, test);
    });

    test('DPBEQU', () {
      srnamc.SRNAMT = 'DPBEQU';
      infoc.INFOT = 1;
      dpbequ('/', 0, 0, A, 1, R1, RCOND, ANRM, INFO);
      chkxer('DPBEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 2;
      dpbequ('U', -1, 0, A, 1, R1, RCOND, ANRM, INFO);
      chkxer('DPBEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 3;
      dpbequ('U', 1, -1, A, 1, R1, RCOND, ANRM, INFO);
      chkxer('DPBEQU', infoc.INFOT, NOUT, LERR, OK, test);
      infoc.INFOT = 5;
      dpbequ('U', 2, 1, A, 1, R1, RCOND, ANRM, INFO);
      chkxer('DPBEQU', infoc.INFOT, NOUT, LERR, OK, test);
    });
  }

  // Print a summary line.
  alaesm(PATH, OK.value, NOUT);
}
