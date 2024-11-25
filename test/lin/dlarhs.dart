// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dgbmv.dart';
import 'package:dart_lapack/src/blas/dgemm.dart';
import 'package:dart_lapack/src/blas/dsbmv.dart';
import 'package:dart_lapack/src/blas/dspmv.dart';
import 'package:dart_lapack/src/blas/dsymm.dart';
import 'package:dart_lapack/src/blas/dtbmv.dart';
import 'package:dart_lapack/src/blas/dtpmv.dart';
import 'package:dart_lapack/src/blas/dtrmm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlarhs(
  final String PATH,
  final String XTYPE,
  final String UPLO,
  final String TRANS,
  final int M,
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> B_,
  final int LDB,
  final Array<int> ISEED_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final ISEED = ISEED_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool BAND, GEN, NOTRAN, QRS, SYM, TRAN, TRI;
  String C1, DIAG;
  String C2;
  int J;

  // Test the input parameters.

  INFO.value = 0;
  C1 = PATH[0];
  C2 = PATH.substring(1, 3);
  TRAN = lsame(TRANS, 'T') || lsame(TRANS, 'C');
  NOTRAN = !TRAN;
  GEN = lsame(PATH[1], 'G');
  QRS = lsame(PATH[1], 'Q') || lsame(PATH[2], 'Q');
  SYM = lsame(PATH[1], 'P') || lsame(PATH[1], 'S');
  TRI = lsame(PATH[1], 'T');
  BAND = lsame(PATH[2], 'B');
  if (!lsame(C1, 'Double precision')) {
    INFO.value = -1;
  } else if (!(lsame(XTYPE, 'N') || lsame(XTYPE, 'C'))) {
    INFO.value = -2;
  } else if ((SYM || TRI) && !(lsame(UPLO, 'U') || lsame(UPLO, 'L'))) {
    INFO.value = -3;
  } else if ((GEN || QRS) && !(TRAN || lsame(TRANS, 'N'))) {
    INFO.value = -4;
  } else if (M < 0) {
    INFO.value = -5;
  } else if (N < 0) {
    INFO.value = -6;
  } else if (BAND && KL < 0) {
    INFO.value = -7;
  } else if (BAND && KU < 0) {
    INFO.value = -8;
  } else if (NRHS < 0) {
    INFO.value = -9;
  } else if ((!BAND && LDA < max(1, M)) ||
      (BAND && (SYM || TRI) && LDA < KL + 1) ||
      (BAND && GEN && LDA < KL + KU + 1)) {
    INFO.value = -11;
  } else if ((NOTRAN && LDX < max(1, N)) || (TRAN && LDX < max(1, M))) {
    INFO.value = -13;
  } else if ((NOTRAN && LDB < max(1, M)) || (TRAN && LDB < max(1, N))) {
    INFO.value = -15;
  }
  if (INFO.value != 0) {
    xerbla('DLARHS', -INFO.value);
    return;
  }

  // Initialize X to NRHS random vectors unless XTYPE = 'C'.

  final (NX, MB) = TRAN ? (M, N) : (N, M);
  if (!lsame(XTYPE, 'C')) {
    for (J = 1; J <= NRHS; J++) {
      dlarnv(2, ISEED, N, X(1, J).asArray());
    }
  }

  // Multiply X by op(A) using an appropriate
  // matrix multiply routine.

  if (lsamen(2, C2, 'GE') ||
      lsamen(2, C2, 'QR') ||
      lsamen(2, C2, 'LQ') ||
      lsamen(2, C2, 'QL') ||
      lsamen(2, C2, 'RQ')) {
    // General matrix

    dgemm(TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, ZERO, B, LDB);
  } else if (lsamen(2, C2, 'PO') || lsamen(2, C2, 'SY')) {
    // Symmetric matrix, 2-D storage

    dsymm('Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB);
  } else if (lsamen(2, C2, 'GB')) {
    // General matrix, band storage

    for (J = 1; J <= NRHS; J++) {
      dgbmv(TRANS, MB, NX, KL, KU, ONE, A, LDA, X(1, J).asArray(), 1, ZERO,
          B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'PB')) {
    // Symmetric matrix, band storage

    for (J = 1; J <= NRHS; J++) {
      dsbmv(UPLO, N, KL, ONE, A, LDA, X(1, J).asArray(), 1, ZERO,
          B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'PP') || lsamen(2, C2, 'SP')) {
    // Symmetric matrix, packed storage

    for (J = 1; J <= NRHS; J++) {
      dspmv(UPLO, N, ONE, A.asArray(), X(1, J).asArray(), 1, ZERO,
          B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'TR')) {
    // Triangular matrix.  Note that for triangular matrices,
    //    KU = 1 => non-unit triangular
    //    KU = 2 => unit triangular

    dlacpy('Full', N, NRHS, X, LDX, B, LDB);
    if (KU == 2) {
      DIAG = 'U';
    } else {
      DIAG = 'N';
    }
    dtrmm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB);
  } else if (lsamen(2, C2, 'TP')) {
    // Triangular matrix, packed storage

    dlacpy('Full', N, NRHS, X, LDX, B, LDB);
    if (KU == 2) {
      DIAG = 'U';
    } else {
      DIAG = 'N';
    }
    for (J = 1; J <= NRHS; J++) {
      dtpmv(UPLO, TRANS, DIAG, N, A.asArray(), B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'TB')) {
    // Triangular matrix, banded storage

    dlacpy('Full', N, NRHS, X, LDX, B, LDB);
    if (KU == 2) {
      DIAG = 'U';
    } else {
      DIAG = 'N';
    }
    for (J = 1; J <= NRHS; J++) {
      dtbmv(UPLO, TRANS, DIAG, N, KL, A, LDA, B(1, J).asArray(), 1);
    }
  } else {
    // If PATH is none of the above, return with an error code.

    INFO.value = -1;
    xerbla('DLARHS', -INFO.value);
  }
}
