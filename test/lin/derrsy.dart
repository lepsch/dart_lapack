// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dspcon.dart';
import 'package:dart_lapack/src/dsprfs.dart';
import 'package:dart_lapack/src/dsptrf.dart';
import 'package:dart_lapack/src/dsptri.dart';
import 'package:dart_lapack/src/dsptrs.dart';
import 'package:dart_lapack/src/dsycon.dart';
import 'package:dart_lapack/src/dsycon_3.dart';
import 'package:dart_lapack/src/dsycon_rook.dart';
import 'package:dart_lapack/src/dsyrfs.dart';
import 'package:dart_lapack/src/dsytf2.dart';
import 'package:dart_lapack/src/dsytf2_rk.dart';
import 'package:dart_lapack/src/dsytf2_rook.dart';
import 'package:dart_lapack/src/dsytrf.dart';
import 'package:dart_lapack/src/dsytrf_aa.dart';
import 'package:dart_lapack/src/dsytrf_aa_2stage.dart';
import 'package:dart_lapack/src/dsytrf_rk.dart';
import 'package:dart_lapack/src/dsytrf_rook.dart';
import 'package:dart_lapack/src/dsytri.dart';
import 'package:dart_lapack/src/dsytri2.dart';
import 'package:dart_lapack/src/dsytri2x.dart';
import 'package:dart_lapack/src/dsytri_3.dart';
import 'package:dart_lapack/src/dsytri_3x.dart';
import 'package:dart_lapack/src/dsytri_rook.dart';
import 'package:dart_lapack/src/dsytrs.dart';
import 'package:dart_lapack/src/dsytrs_3.dart';
import 'package:dart_lapack/src/dsytrs_aa.dart';
import 'package:dart_lapack/src/dsytrs_aa_2stage.dart';
import 'package:dart_lapack/src/dsytrs_rook.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';

import '../test_driver.dart';
import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrsy(final String PATH, final Nout NUNIT, final TestDriver test) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4;
  final IP = Array<int>(NMAX), IW = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      E = Array<double>(NMAX),
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
    E[J] = 0.0;
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = 0.0;
    X[J] = 0.0;
    IP[J] = J;
    IW[J] = J;
  }
  ANRM.value = 1.0;
  RCOND.value = 1.0;
  OK.value = true;

  if (lsamen(2, C2, 'SY')) {
    // Test error exits of the routines that use factorization
    // of a symmetric indefinite matrix with partial
    // (Bunch-Kaufman) pivoting.

    test('DSYTRF', () {
      srnamc.SRNAMT = 'DSYTRF';
      infoc.INFOT = 1;
      dsytrf('/', 0, A, 1, IP, W, 1, INFO);
      chkxer('DSYTRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrf('U', -1, A, 1, IP, W, 1, INFO);
      chkxer('DSYTRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrf('U', 2, A, 1, IP, W, 4, INFO);
      chkxer('DSYTRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrf('U', 0, A, 1, IP, W, 0, INFO);
      chkxer('DSYTRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrf('U', 0, A, 1, IP, W, -2, INFO);
      chkxer('DSYTRF', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTF2', () {
      srnamc.SRNAMT = 'DSYTF2';
      infoc.INFOT = 1;
      dsytf2('/', 0, A, 1, IP, INFO);
      chkxer('DSYTF2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytf2('U', -1, A, 1, IP, INFO);
      chkxer('DSYTF2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytf2('U', 2, A, 1, IP, INFO);
      chkxer('DSYTF2', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRI', () {
      srnamc.SRNAMT = 'DSYTRI';
      infoc.INFOT = 1;
      dsytri('/', 0, A, 1, IP, W, INFO);
      chkxer('DSYTRI', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytri('U', -1, A, 1, IP, W, INFO);
      chkxer('DSYTRI', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytri('U', 2, A, 1, IP, W, INFO);
      chkxer('DSYTRI', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRI2', () {
      srnamc.SRNAMT = 'DSYTRI2';
      infoc.INFOT = 1;
      dsytri2('/', 0, A, 1, IP, W, IW[1], INFO);
      chkxer('DSYTRI2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytri2('U', -1, A, 1, IP, W, IW[1], INFO);
      chkxer('DSYTRI2', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytri2('U', 2, A, 1, IP, W, IW[1], INFO);
      chkxer('DSYTRI2', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRI2X', () {
      srnamc.SRNAMT = 'DSYTRI2X';
      infoc.INFOT = 1;
      dsytri2x('/', 0, A, 1, IP, W.asMatrix(), 1, INFO);
      chkxer('DSYTRI2X', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytri2x('U', -1, A, 1, IP, W.asMatrix(), 1, INFO);
      chkxer('DSYTRI2X', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytri2x('U', 2, A, 1, IP, W.asMatrix(), 1, INFO);
      chkxer('DSYTRI2X', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRS', () {
      srnamc.SRNAMT = 'DSYTRS';
      infoc.INFOT = 1;
      dsytrs('/', 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrs('U', -1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrs('U', 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrs('U', 2, 1, A, 1, IP, B.asMatrix(), 2, INFO);
      chkxer('DSYTRS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsytrs('U', 2, 1, A, 2, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYRFS', () {
      srnamc.SRNAMT = 'DSYRFS';
      infoc.INFOT = 1;
      dsyrfs('/', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DSYRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsyrfs('U', -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DSYRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsyrfs('U', 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DSYRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsyrfs('U', 2, 1, A, 1, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DSYRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsyrfs('U', 2, 1, A, 2, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DSYRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsyrfs('U', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2, R1,
          R2, W, IW, INFO);
      chkxer('DSYRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 12;
      dsyrfs('U', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1, R1,
          R2, W, IW, INFO);
      chkxer('DSYRFS', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYCON', () {
      srnamc.SRNAMT = 'DSYCON';
      infoc.INFOT = 1;
      dsycon('/', 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsycon('U', -1, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsycon('U', 2, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsycon('U', 1, A, 1, IP, -1.0, RCOND, W, IW, INFO);
      chkxer('DSYCON', infoc.INFOT, NOUT, LERR, OK);
    });
  } else if (lsamen(2, C2, 'SR')) {
    // Test error exits of the routines that use factorization
    // of a symmetric indefinite matrix with rook
    // (bounded Bunch-Kaufman) pivoting.

    test('DSYTRF_ROOK', () {
      srnamc.SRNAMT = 'DSYTRF_ROOK';
      infoc.INFOT = 1;
      dsytrf_rook('/', 0, A, 1, IP, W, 1, INFO);
      chkxer('DSYTRF_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrf_rook('U', -1, A, 1, IP, W, 1, INFO);
      chkxer('DSYTRF_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrf_rook('U', 2, A, 1, IP, W, 4, INFO);
      chkxer('DSYTRF_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrf_rook('U', 0, A, 1, IP, W, 0, INFO);
      chkxer('DSYTRF_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrf_rook('U', 0, A, 1, IP, W, -2, INFO);
      chkxer('DSYTRF_ROOK', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTF2_ROOK', () {
      srnamc.SRNAMT = 'DSYTF2_ROOK';
      infoc.INFOT = 1;
      dsytf2_rook('/', 0, A, 1, IP, INFO);
      chkxer('DSYTF2_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytf2_rook('U', -1, A, 1, IP, INFO);
      chkxer('DSYTF2_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytf2_rook('U', 2, A, 1, IP, INFO);
      chkxer('DSYTF2_ROOK', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRI_ROOK', () {
      srnamc.SRNAMT = 'DSYTRI_ROOK';
      infoc.INFOT = 1;
      dsytri_rook('/', 0, A, 1, IP, W, INFO);
      chkxer('DSYTRI_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytri_rook('U', -1, A, 1, IP, W, INFO);
      chkxer('DSYTRI_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytri_rook('U', 2, A, 1, IP, W, INFO);
      chkxer('DSYTRI_ROOK', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRS_ROOK', () {
      srnamc.SRNAMT = 'DSYTRS_ROOK';
      infoc.INFOT = 1;
      dsytrs_rook('/', 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrs_rook('U', -1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrs_rook('U', 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrs_rook('U', 2, 1, A, 1, IP, B.asMatrix(), 2, INFO);
      chkxer('DSYTRS_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsytrs_rook('U', 2, 1, A, 2, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_ROOK', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYCON_ROOK', () {
      srnamc.SRNAMT = 'DSYCON_ROOK';
      infoc.INFOT = 1;
      dsycon_rook('/', 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsycon_rook('U', -1, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsycon_rook('U', 2, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON_ROOK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsycon_rook('U', 1, A, 1, IP, -1.0, RCOND, W, IW, INFO);
      chkxer('DSYCON_ROOK', infoc.INFOT, NOUT, LERR, OK);
    });
  } else if (lsamen(2, C2, 'SK')) {
    // Test error exits of the routines that use factorization
    // of a symmetric indefinite matrix with rook
    // (bounded Bunch-Kaufman) pivoting with the new storage
    // format for factors L ( or U) and D.

    // L (or U) is stored in A, diagonal of D is stored on the
    // diagonal of A, subdiagonal of D is stored in a separate array E.

    test('DSYTRF_RK', () {
      srnamc.SRNAMT = 'DSYTRF_RK';
      infoc.INFOT = 1;
      dsytrf_rk('/', 0, A, 1, E, IP, W, 1, INFO);
      chkxer('DSYTRF_RK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrf_rk('U', -1, A, 1, E, IP, W, 1, INFO);
      chkxer('DSYTRF_RK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrf_rk('U', 2, A, 1, E, IP, W, 1, INFO);
      chkxer('DSYTRF_RK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsytrf_rk('U', 0, A, 1, E, IP, W, 0, INFO);
      chkxer('DSYTRF_RK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsytrf_rk('U', 0, A, 1, E, IP, W, -2, INFO);
      chkxer('DSYTRF_RK', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTF2_RK', () {
      srnamc.SRNAMT = 'DSYTF2_RK';
      infoc.INFOT = 1;
      dsytf2_rk('/', 0, A, 1, E, IP, INFO);
      chkxer('DSYTF2_RK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytf2_rk('U', -1, A, 1, E, IP, INFO);
      chkxer('DSYTF2_RK', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytf2_rk('U', 2, A, 1, E, IP, INFO);
      chkxer('DSYTF2_RK', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRI_3', () {
      srnamc.SRNAMT = 'DSYTRI_3';
      infoc.INFOT = 1;
      dsytri_3('/', 0, A, 1, E, IP, W, 1, INFO);
      chkxer('DSYTRI_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytri_3('U', -1, A, 1, E, IP, W, 1, INFO);
      chkxer('DSYTRI_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytri_3('U', 2, A, 1, E, IP, W, 1, INFO);
      chkxer('DSYTRI_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsytri_3('U', 0, A, 1, E, IP, W, 0, INFO);
      chkxer('DSYTRI_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsytri_3('U', 0, A, 1, E, IP, W, -2, INFO);
      chkxer('DSYTRI_3', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRI_3X', () {
      srnamc.SRNAMT = 'DSYTRI_3X';
      infoc.INFOT = 1;
      dsytri_3x('/', 0, A, 1, E, IP, W.asMatrix(), 1, INFO);
      chkxer('DSYTRI_3X', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytri_3x('U', -1, A, 1, E, IP, W.asMatrix(), 1, INFO);
      chkxer('DSYTRI_3X', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytri_3x('U', 2, A, 1, E, IP, W.asMatrix(), 1, INFO);
      chkxer('DSYTRI_3X', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRS_3', () {
      srnamc.SRNAMT = 'DSYTRS_3';
      infoc.INFOT = 1;
      dsytrs_3('/', 0, 0, A, 1, E, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrs_3('U', -1, 0, A, 1, E, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrs_3('U', 0, -1, A, 1, E, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrs_3('U', 2, 1, A, 1, E, IP, B.asMatrix(), 2, INFO);
      chkxer('DSYTRS_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 9;
      dsytrs_3('U', 2, 1, A, 2, E, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_3', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYCON_3', () {
      srnamc.SRNAMT = 'DSYCON_3';
      infoc.INFOT = 1;
      dsycon_3('/', 0, A, 1, E, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsycon_3('U', -1, A, 1, E, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsycon_3('U', 2, A, 1, E, IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSYCON_3', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsycon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, IW, INFO);
      chkxer('DSYCON_3', infoc.INFOT, NOUT, LERR, OK);
    });
  } else if (lsamen(2, C2, 'SA')) {
    // Test error exits of the routines that use factorization
    // of a symmetric indefinite matrix with Aasen's algorithm.

    test('DSYTRF_AA', () {
      srnamc.SRNAMT = 'DSYTRF_AA';
      infoc.INFOT = 1;
      dsytrf_aa('/', 0, A, 1, IP, W, 1, INFO);
      chkxer('DSYTRF_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrf_aa('U', -1, A, 1, IP, W, 1, INFO);
      chkxer('DSYTRF_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrf_aa('U', 2, A, 1, IP, W, 4, INFO);
      chkxer('DSYTRF_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrf_aa('U', 0, A, 1, IP, W, 0, INFO);
      chkxer('DSYTRF_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrf_aa('U', 0, A, 1, IP, W, -2, INFO);
      chkxer('DSYTRF_AA', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRS_AA', () {
      srnamc.SRNAMT = 'DSYTRS_AA';
      infoc.INFOT = 1;
      dsytrs_aa('/', 0, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYTRS_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrs_aa('U', -1, 0, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYTRS_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrs_aa('U', 0, -1, A, 1, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYTRS_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrs_aa('U', 2, 1, A, 1, IP, B.asMatrix(), 2, W, 1, INFO);
      chkxer('DSYTRS_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsytrs_aa('U', 2, 1, A, 2, IP, B.asMatrix(), 1, W, 1, INFO);
      chkxer('DSYTRS_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsytrs_aa('U', 0, 1, A, 2, IP, B.asMatrix(), 1, W, 0, INFO);
      chkxer('DSYTRS_AA', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsytrs_aa('U', 0, 1, A, 2, IP, B.asMatrix(), 1, W, -2, INFO);
      chkxer('DSYTRS_AA', infoc.INFOT, NOUT, LERR, OK);
    });
  } else if (lsamen(2, C2, 'S2')) {
    // Test error exits of the routines that use factorization
    // of a symmetric indefinite matrix with Aasen's algorithm.

    test('DSYTRF_AA_2STAGE', () {
      srnamc.SRNAMT = 'DSYTRF_AA_2STAGE';
      infoc.INFOT = 1;
      dsytrf_aa_2stage('/', 0, A, 1, A.asArray(), 1, IP, IP, W, 1, INFO);
      chkxer('DSYTRF_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrf_aa_2stage('U', -1, A, 1, A.asArray(), 1, IP, IP, W, 1, INFO);
      chkxer('DSYTRF_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 4;
      dsytrf_aa_2stage('U', 2, A, 1, A.asArray(), 2, IP, IP, W, 1, INFO);
      chkxer('DSYTRF_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 6;
      dsytrf_aa_2stage('U', 2, A, 2, A.asArray(), 1, IP, IP, W, 1, INFO);
      chkxer('DSYTRF_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsytrf_aa_2stage('U', 2, A, 2, A.asArray(), 8, IP, IP, W, 0, INFO);
      chkxer('DSYTRF_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSYTRS_AA_2STAGE', () {
      srnamc.SRNAMT = 'DSYTRS_AA_2STAGE';
      infoc.INFOT = 1;
      dsytrs_aa_2stage(
          '/', 0, 0, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsytrs_aa_2stage(
          'U', -1, 0, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsytrs_aa_2stage(
          'U', 0, -1, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dsytrs_aa_2stage(
          'U', 2, 1, A, 1, A.asArray(), 1, IP, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsytrs_aa_2stage(
          'U', 2, 1, A, 2, A.asArray(), 1, IP, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_AA_2STAGE', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 11;
      dsytrs_aa_2stage(
          'U', 2, 1, A, 2, A.asArray(), 8, IP, IP, B.asMatrix(), 1, INFO);
      chkxer('DSYTRS_AA_STAGE', infoc.INFOT, NOUT, LERR, OK);
    });
  } else if (lsamen(2, C2, 'SP')) {
    // Test error exits of the routines that use factorization
    // of a symmetric indefinite packed matrix with partial
    // (Bunch-Kaufman) pivoting.

    test('DSPTRF', () {
      srnamc.SRNAMT = 'DSPTRF';
      infoc.INFOT = 1;
      dsptrf('/', 0, A.asArray(), IP, INFO);
      chkxer('DSPTRF', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsptrf('U', -1, A.asArray(), IP, INFO);
      chkxer('DSPTRF', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSPTRI', () {
      srnamc.SRNAMT = 'DSPTRI';
      infoc.INFOT = 1;
      dsptri('/', 0, A.asArray(), IP, W, INFO);
      chkxer('DSPTRI', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsptri('U', -1, A.asArray(), IP, W, INFO);
      chkxer('DSPTRI', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSPTRS', () {
      srnamc.SRNAMT = 'DSPTRS';
      infoc.INFOT = 1;
      dsptrs('/', 0, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPTRS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsptrs('U', -1, 0, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPTRS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsptrs('U', 0, -1, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPTRS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 7;
      dsptrs('U', 2, 1, A.asArray(), IP, B.asMatrix(), 1, INFO);
      chkxer('DSPTRS', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSPRFS', () {
      srnamc.SRNAMT = 'DSPRFS';
      infoc.INFOT = 1;
      dsprfs('/', 0, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DSPRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dsprfs('U', -1, 0, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DSPRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 3;
      dsprfs('U', 0, -1, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DSPRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 8;
      dsprfs('U', 2, 1, A.asArray(), AF.asArray(), IP, B.asMatrix(), 1,
          X.asMatrix(), 2, R1, R2, W, IW, INFO);
      chkxer('DSPRFS', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 10;
      dsprfs('U', 2, 1, A.asArray(), AF.asArray(), IP, B.asMatrix(), 2,
          X.asMatrix(), 1, R1, R2, W, IW, INFO);
      chkxer('DSPRFS', infoc.INFOT, NOUT, LERR, OK);
    });

    test('DSPCON', () {
      srnamc.SRNAMT = 'DSPCON';
      infoc.INFOT = 1;
      dspcon('/', 0, A.asArray(), IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSPCON', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 2;
      dspcon('U', -1, A.asArray(), IP, ANRM.value, RCOND, W, IW, INFO);
      chkxer('DSPCON', infoc.INFOT, NOUT, LERR, OK);
      infoc.INFOT = 5;
      dspcon('U', 1, A.asArray(), IP, -1.0, RCOND, W, IW, INFO);
      chkxer('DSPCON', infoc.INFOT, NOUT, LERR, OK);
    });
  }

  // Print a summary line.
  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
