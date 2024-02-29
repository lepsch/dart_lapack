import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zbdsqr.dart';
import 'package:lapack/src/zgebd2.dart';
import 'package:lapack/src/zgebrd.dart';
import 'package:lapack/src/zungbr.dart';
import 'package:lapack/src/zunmbr.dart';

import 'chkxer.dart';
import 'common.dart';

void zerrbd(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4, LW = NMAX;
  String C2;
  int I, J, NT;
  final D = Array<double>(NMAX),
      E = Array<double>(NMAX),
      RW = Array<double>(4 * NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      U = Matrix<Complex>(NMAX, NMAX),
      V = Matrix<Complex>(NMAX, NMAX);
  final TP = Array<Complex>(NMAX),
      TQ = Array<Complex>(NMAX),
      W = Array<Complex>(LW);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();
  C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (J = 1; J <= NMAX; J++) {
    // 20
    for (I = 1; I <= NMAX; I++) {
      // 10
      A[I][J] = (1.0 / (I + J)).toComplex();
    } // 10
  } // 20
  infoc.OK.value = true;
  NT = 0;

  // Test error exits of the SVD routines.

  if (lsamen(2, C2, 'BD')) {
    // ZGEBRD

    srnamc.SRNAMT = 'ZGEBRD';
    infoc.INFOT = 1;
    zgebrd(-1, 0, A, 1, D, E, TQ, TP, W, 1, INFO);
    chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgebrd(0, -1, A, 1, D, E, TQ, TP, W, 1, INFO);
    chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgebrd(2, 1, A, 1, D, E, TQ, TP, W, 2, INFO);
    chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zgebrd(2, 1, A, 2, D, E, TQ, TP, W, 1, INFO);
    chkxer('ZGEBRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = NT + 4;

    // ZGEBD2

    srnamc.SRNAMT = 'ZGEBD2';
    infoc.INFOT = 1;
    zgebd2(-1, 0, A, 1, D, E, TQ, TP, W, INFO);
    chkxer('ZGEBD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgebd2(0, -1, A, 1, D, E, TQ, TP, W, INFO);
    chkxer('ZGEBD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgebd2(2, 1, A, 1, D, E, TQ, TP, W, INFO);
    chkxer('ZGEBD2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = NT + 3;

    // ZUNGBR

    srnamc.SRNAMT = 'ZUNGBR';
    infoc.INFOT = 1;
    zungbr('/', 0, 0, 0, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zungbr('Q', -1, 0, 0, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zungbr('Q', 0, -1, 0, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zungbr('Q', 0, 1, 0, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zungbr('Q', 1, 0, 1, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zungbr('P', 1, 0, 0, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zungbr('P', 0, 1, 1, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zungbr('Q', 0, 0, -1, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zungbr('Q', 2, 1, 1, A, 1, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zungbr('Q', 2, 2, 1, A, 2, TQ, W, 1, INFO);
    chkxer('ZUNGBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = NT + 10;

    // ZUNMBR

    srnamc.SRNAMT = 'ZUNMBR';
    infoc.INFOT = 1;
    zunmbr('/', 'L', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zunmbr('Q', '/', 'T', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zunmbr('Q', 'L', '/', 0, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zunmbr('Q', 'L', 'C', -1, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zunmbr('Q', 'L', 'C', 0, -1, 0, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zunmbr('Q', 'L', 'C', 0, 0, -1, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zunmbr('Q', 'L', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zunmbr('Q', 'R', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zunmbr('P', 'L', 'C', 2, 0, 2, A, 1, TQ, U, 2, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zunmbr('P', 'R', 'C', 0, 2, 2, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zunmbr('Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 1, W, 1, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zunmbr('Q', 'L', 'C', 0, 2, 0, A, 1, TQ, U, 1, W, 0, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zunmbr('Q', 'R', 'C', 2, 0, 0, A, 1, TQ, U, 2, W, 0, INFO);
    chkxer('ZUNMBR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = NT + 13;

    // ZBDSQR

    srnamc.SRNAMT = 'ZBDSQR';
    infoc.INFOT = 1;
    zbdsqr('/', 0, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zbdsqr('U', -1, 0, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zbdsqr('U', 0, -1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zbdsqr('U', 0, 0, -1, 0, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zbdsqr('U', 0, 0, 0, -1, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zbdsqr('U', 2, 1, 0, 0, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zbdsqr('U', 0, 0, 2, 0, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zbdsqr('U', 2, 0, 0, 1, D, E, V, 1, U, 1, A, 1, RW, INFO);
    chkxer('ZBDSQR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT = NT + 8;
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
