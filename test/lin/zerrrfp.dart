// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zhfrk.dart';
import 'package:lapack/src/zpftrf.dart';
import 'package:lapack/src/zpftri.dart';
import 'package:lapack/src/zpftrs.dart';
import 'package:lapack/src/ztfsm.dart';
import 'package:lapack/src/ztftri.dart';
import 'package:lapack/src/ztfttp.dart';
import 'package:lapack/src/ztfttr.dart';
import 'package:lapack/src/ztpttf.dart';
import 'package:lapack/src/ztpttr.dart';
import 'package:lapack/src/ztrttf.dart';
import 'package:lapack/src/ztrttp.dart';

import 'chkxer.dart';
import 'common.dart';

void zerrrfp(final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  double ALPHA, BETA;
  Complex CALPHA;
  final A = Matrix<Complex>(1, 1), B = Matrix<Complex>(1, 1);
  final INFO = Box(0);

  final NOUT = NUNIT;
  infoc.OK.value = true;
  A[1][1] = Complex(1.0, 1.0);
  B[1][1] = Complex(1.0, 1.0);
  ALPHA = 1.0;
  CALPHA = Complex(1.0, 1.0);
  BETA = 1.0;

  srnamc.SRNAMT = 'ZPFTRF';
  infoc.INFOT = 1;
  zpftrf('/', 'U', 0, A.asArray(), INFO);
  chkxer('ZPFTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zpftrf('N', '/', 0, A.asArray(), INFO);
  chkxer('ZPFTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zpftrf('N', 'U', -1, A.asArray(), INFO);
  chkxer('ZPFTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZPFTRS';
  infoc.INFOT = 1;
  zpftrs('/', 'U', 0, 0, A.asArray(), B, 1, INFO);
  chkxer('ZPFTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zpftrs('N', '/', 0, 0, A.asArray(), B, 1, INFO);
  chkxer('ZPFTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zpftrs('N', 'U', -1, 0, A.asArray(), B, 1, INFO);
  chkxer('ZPFTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zpftrs('N', 'U', 0, -1, A.asArray(), B, 1, INFO);
  chkxer('ZPFTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  zpftrs('N', 'U', 0, 0, A.asArray(), B, 0, INFO);
  chkxer('ZPFTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZPFTRI';
  infoc.INFOT = 1;
  zpftri('/', 'U', 0, A.asArray(), INFO);
  chkxer('ZPFTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zpftri('N', '/', 0, A.asArray(), INFO);
  chkxer('ZPFTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zpftri('N', 'U', -1, A.asArray(), INFO);
  chkxer('ZPFTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTFSM';
  infoc.INFOT = 1;
  ztfsm('/', 'L', 'U', 'C', 'U', 0, 0, CALPHA, A.asArray(), B, 1);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztfsm('N', '/', 'U', 'C', 'U', 0, 0, CALPHA, A.asArray(), B, 1);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztfsm('N', 'L', '/', 'C', 'U', 0, 0, CALPHA, A.asArray(), B, 1);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztfsm('N', 'L', 'U', '/', 'U', 0, 0, CALPHA, A.asArray(), B, 1);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztfsm('N', 'L', 'U', 'C', '/', 0, 0, CALPHA, A.asArray(), B, 1);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  ztfsm('N', 'L', 'U', 'C', 'U', -1, 0, CALPHA, A.asArray(), B, 1);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 7;
  ztfsm('N', 'L', 'U', 'C', 'U', 0, -1, CALPHA, A.asArray(), B, 1);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 11;
  ztfsm('N', 'L', 'U', 'C', 'U', 0, 0, CALPHA, A.asArray(), B, 0);
  chkxer('ZTFSM', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTFTRI';
  infoc.INFOT = 1;
  ztftri('/', 'L', 'N', 0, A.asArray(), INFO);
  chkxer('ZTFTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztftri('N', '/', 'N', 0, A.asArray(), INFO);
  chkxer('ZTFTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztftri('N', 'L', '/', 0, A.asArray(), INFO);
  chkxer('ZTFTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztftri('N', 'L', 'N', -1, A.asArray(), INFO);
  chkxer('ZTFTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTFTTR';
  infoc.INFOT = 1;
  ztfttr('/', 'U', 0, A.asArray(), B, 1, INFO);
  chkxer('ZTFTTR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztfttr('N', '/', 0, A.asArray(), B, 1, INFO);
  chkxer('ZTFTTR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztfttr('N', 'U', -1, A.asArray(), B, 1, INFO);
  chkxer('ZTFTTR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 6;
  ztfttr('N', 'U', 0, A.asArray(), B, 0, INFO);
  chkxer('ZTFTTR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTRTTF';
  infoc.INFOT = 1;
  ztrttf('/', 'U', 0, A, 1, B.asArray(), INFO);
  chkxer('ZTRTTF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztrttf('N', '/', 0, A, 1, B.asArray(), INFO);
  chkxer('ZTRTTF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztrttf('N', 'U', -1, A, 1, B.asArray(), INFO);
  chkxer('ZTRTTF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztrttf('N', 'U', 0, A, 0, B.asArray(), INFO);
  chkxer('ZTRTTF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTFTTP';
  infoc.INFOT = 1;
  ztfttp('/', 'U', 0, A.asArray(), B.asArray(), INFO);
  chkxer('ZTFTTP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztfttp('N', '/', 0, A.asArray(), B.asArray(), INFO);
  chkxer('ZTFTTP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztfttp('N', 'U', -1, A.asArray(), B.asArray(), INFO);
  chkxer('ZTFTTP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTPTTF';
  infoc.INFOT = 1;
  ztpttf('/', 'U', 0, A.asArray(), B.asArray(), INFO);
  chkxer('ZTPTTF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztpttf('N', '/', 0, A.asArray(), B.asArray(), INFO);
  chkxer('ZTPTTF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  ztpttf('N', 'U', -1, A.asArray(), B.asArray(), INFO);
  chkxer('ZTPTTF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTRTTP';
  infoc.INFOT = 1;
  ztrttp('/', 0, A, 1, B.asArray(), INFO);
  chkxer('ZTRTTP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztrttp('U', -1, A, 1, B.asArray(), INFO);
  chkxer('ZTRTTP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  ztrttp('U', 0, A, 0, B.asArray(), INFO);
  chkxer('ZTRTTP', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZTPTTR';
  infoc.INFOT = 1;
  ztpttr('/', 0, A.asArray(), B, 1, INFO);
  chkxer('ZTPTTR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  ztpttr('U', -1, A.asArray(), B, 1, INFO);
  chkxer('ZTPTTR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  ztpttr('U', 0, A.asArray(), B, 0, INFO);
  chkxer('ZTPTTR', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  srnamc.SRNAMT = 'ZHFRK';
  infoc.INFOT = 1;
  zhfrk('/', 'U', 'N', 0, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('ZHFRK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 2;
  zhfrk('N', '/', 'N', 0, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('ZHFRK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 3;
  zhfrk('N', 'U', '/', 0, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('ZHFRK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 4;
  zhfrk('N', 'U', 'N', -1, 0, ALPHA, A, 1, BETA, B.asArray());
  chkxer('ZHFRK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 5;
  zhfrk('N', 'U', 'N', 0, -1, ALPHA, A, 1, BETA, B.asArray());
  chkxer('ZHFRK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  infoc.INFOT = 8;
  zhfrk('N', 'U', 'N', 0, 0, ALPHA, A, 0, BETA, B.asArray());
  chkxer('ZHFRK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK);

  // Print a summary line.

  if (infoc.OK.value) {
    NOUT.println(' Complex RFP routines passed the tests of the error exits');
  } else {
    NOUT.println(' *** RFP routines failed the tests of the error exits ***');
  }
}
