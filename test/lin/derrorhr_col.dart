import 'package:lapack/src/box.dart';
import 'package:lapack/src/dorhr_col.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrorhr_col(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 2;
  final A = Matrix<double>(NMAX, NMAX),
      T = Matrix<double>(NMAX, NMAX),
      D = Array<double>(NMAX);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      T[I][J] = 1.0 / (I + J);
    }
    D[J] = 0.0;
  }
  infoc.OK.value = true;

  // Error exits for Householder reconstruction

  // DORHR_COL

  srnamc.SRNAMT = 'DORHR_COL';

  infoc.INFOT = 1;
  dorhr_col(-1, 0, 1, A, 1, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 2;
  dorhr_col(0, -1, 1, A, 1, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  dorhr_col(1, 2, 1, A, 1, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 3;
  dorhr_col(0, 0, -1, A, 1, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  dorhr_col(0, 0, 0, A, 1, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 5;
  dorhr_col(0, 0, 1, A, -1, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  dorhr_col(0, 0, 1, A, 0, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  dorhr_col(2, 0, 1, A, 1, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  infoc.INFOT = 7;
  dorhr_col(0, 0, 1, A, 1, T, -1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  dorhr_col(0, 0, 1, A, 1, T, 0, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  dorhr_col(4, 3, 2, A, 4, T, 1, D, INFO);
  chkxer('DORHR_COL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
