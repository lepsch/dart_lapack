// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dpftrf.dart';
import 'package:dart_lapack/src/dpftri.dart';
import 'package:dart_lapack/src/dpftrs.dart';
import 'package:dart_lapack/src/dsfrk.dart';
import 'package:dart_lapack/src/dtfsm.dart';
import 'package:dart_lapack/src/dtftri.dart';
import 'package:dart_lapack/src/dtfttp.dart';
import 'package:dart_lapack/src/dtfttr.dart';
import 'package:dart_lapack/src/dtpttf.dart';
import 'package:dart_lapack/src/dtpttr.dart';
import 'package:dart_lapack/src/dtrttf.dart';
import 'package:dart_lapack/src/dtrttp.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../eig/chkxer.dart';
import 'common.dart';

void derrrfp(final Nout NUNIT) {
  final A = Matrix<double>(1, 1), B = Matrix<double>(1, 1);
  final INFO = Box(0);

  infoc.NOUT = NUNIT;
  infoc.OK.value = true;
  A[1][1] = 1.0;
  B[1][1] = 1.0;
  final ALPHA = 1.0;
  final BETA = 1.0;

  srnamc.SRNAMT = 'DPFTRF';
  infoc.INFOT = 1;
  dpftrf('/', 'U', 0, A.asArray(), INFO);
  chkxer('DPFTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dpftrf('N', '/', 0, A.asArray(), INFO);
  chkxer('DPFTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dpftrf('N', 'U', -1, A.asArray(), INFO);
  chkxer('DPFTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DPFTRS';
  infoc.INFOT = 1;
  dpftrs('/', 'U', 0, 0, A.asArray(), B, 1, INFO);
  chkxer('DPFTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dpftrs('N', '/', 0, 0, A.asArray(), B, 1, INFO);
  chkxer('DPFTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dpftrs('N', 'U', -1, 0, A.asArray(), B, 1, INFO);
  chkxer('DPFTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dpftrs('N', 'U', 0, -1, A.asArray(), B, 1, INFO);
  chkxer('DPFTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dpftrs('N', 'U', 0, 0, A.asArray(), B, 0, INFO);
  chkxer('DPFTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DPFTRI';
  infoc.INFOT = 1;
  dpftri('/', 'U', 0, A.asArray(), INFO);
  chkxer('DPFTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dpftri('N', '/', 0, A.asArray(), INFO);
  chkxer('DPFTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dpftri('N', 'U', -1, A.asArray(), INFO);
  chkxer('DPFTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTFSM';
  infoc.INFOT = 1;
  dtfsm('/', 'L', 'U', 'T', 'U', 0, 0, ALPHA, A.asArray(), B, 1);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtfsm('N', '/', 'U', 'T', 'U', 0, 0, ALPHA, A.asArray(), B, 1);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtfsm('N', 'L', '/', 'T', 'U', 0, 0, ALPHA, A.asArray(), B, 1);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtfsm('N', 'L', 'U', '/', 'U', 0, 0, ALPHA, A.asArray(), B, 1);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dtfsm('N', 'L', 'U', 'T', '/', 0, 0, ALPHA, A.asArray(), B, 1);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dtfsm('N', 'L', 'U', 'T', 'U', -1, 0, ALPHA, A.asArray(), B, 1);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  dtfsm('N', 'L', 'U', 'T', 'U', 0, -1, ALPHA, A.asArray(), B, 1);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  dtfsm('N', 'L', 'U', 'T', 'U', 0, 0, ALPHA, A.asArray(), B, 0);
  chkxer('DTFSM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTFTRI';
  infoc.INFOT = 1;
  dtftri('/', 'L', 'N', 0, A.asArray(), INFO);
  chkxer('DTFTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtftri('N', '/', 'N', 0, A.asArray(), INFO);
  chkxer('DTFTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtftri('N', 'L', '/', 0, A.asArray(), INFO);
  chkxer('DTFTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtftri('N', 'L', 'N', -1, A.asArray(), INFO);
  chkxer('DTFTRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTFTTR';
  infoc.INFOT = 1;
  dtfttr('/', 'U', 0, A.asArray(), B, 1, INFO);
  chkxer('DTFTTR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtfttr('N', '/', 0, A.asArray(), B, 1, INFO);
  chkxer('DTFTTR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtfttr('N', 'U', -1, A.asArray(), B, 1, INFO);
  chkxer('DTFTTR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  dtfttr('N', 'U', 0, A.asArray(), B, 0, INFO);
  chkxer('DTFTTR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTRTTF';
  infoc.INFOT = 1;
  dtrttf('/', 'U', 0, A, 1, B.asArray(), INFO);
  chkxer('DTRTTF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtrttf('N', '/', 0, A, 1, B.asArray(), INFO);
  chkxer('DTRTTF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtrttf('N', 'U', -1, A, 1, B.asArray(), INFO);
  chkxer('DTRTTF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dtrttf('N', 'U', 0, A, 0, B.asArray(), INFO);
  chkxer('DTRTTF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTFTTP';
  infoc.INFOT = 1;
  dtfttp('/', 'U', 0, A.asArray(), B.asArray(), INFO);
  chkxer('DTFTTP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtfttp('N', '/', 0, A.asArray(), B.asArray(), INFO);
  chkxer('DTFTTP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtfttp('N', 'U', -1, A.asArray(), B.asArray(), INFO);
  chkxer('DTFTTP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTPTTF';
  infoc.INFOT = 1;
  dtpttf('/', 'U', 0, A.asArray(), B.asArray(), INFO);
  chkxer('DTPTTF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtpttf('N', '/', 0, A.asArray(), B.asArray(), INFO);
  chkxer('DTPTTF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dtpttf('N', 'U', -1, A.asArray(), B.asArray(), INFO);
  chkxer('DTPTTF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTRTTP';
  infoc.INFOT = 1;
  dtrttp('/', 0, A, 1, B.asArray(), INFO);
  chkxer('DTRTTP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtrttp('U', -1, A, 1, B.asArray(), INFO);
  chkxer('DTRTTP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dtrttp('U', 0, A, 0, B.asArray(), INFO);
  chkxer('DTRTTP', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DTPTTR';
  infoc.INFOT = 1;
  dtpttr('/', 0, A.asArray(), B, 1, INFO);
  chkxer('DTPTTR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dtpttr('U', -1, A.asArray(), B, 1, INFO);
  chkxer('DTPTTR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dtpttr('U', 0, A.asArray(), B, 0, INFO);
  chkxer('DTPTTR', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'DSFRK';
  infoc.INFOT = 1;
  dsfrk('/', 'U', 'N', 0, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('DSFRK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  dsfrk('N', '/', 'N', 0, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('DSFRK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  dsfrk('N', 'U', '/', 0, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('DSFRK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  dsfrk('N', 'U', 'N', -1, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('DSFRK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  dsfrk('N', 'U', 'N', 0, -1, ALPHA, A, 1, BETA, B.asArray());
  chkxer('DSFRK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  dsfrk('N', 'U', 'N', 0, 0, ALPHA, A, 0, BETA, B.asArray());
  chkxer('DSFRK', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT.println(
        ' DOUBLE PRECISION RFP routines passed the tests of the error exits');
  } else {
    infoc.NOUT
        .println(' *** RFP routines failed the tests of the error exits ***');
  }
}
