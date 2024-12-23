// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zpstf2.dart';
import 'package:dart_lapack/src/zpstrf.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void zerrps(final String PATH, final Nout NUNIT) {
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
