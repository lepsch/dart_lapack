// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgbmv.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zhbmv.dart';
import 'package:dart_lapack/src/blas/zhemm.dart';
import 'package:dart_lapack/src/blas/zhpmv.dart';
import 'package:dart_lapack/src/blas/zsymm.dart';
import 'package:dart_lapack/src/blas/ztbmv.dart';
import 'package:dart_lapack/src/blas/ztpmv.dart';
import 'package:dart_lapack/src/blas/ztrmm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/lsamen.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlarnv.dart';
import 'package:dart_lapack/src/zspmv.dart';

import 'zsbmv.dart';

void zlarhs(
  final String PATH,
  final String XTYPE,
  final String UPLO,
  final String TRANS,
  final int M,
  final int N,
  final int KL,
  final int KU,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<int> ISEED_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final B = B_.having(ld: LDB);
  final ISEED = ISEED_.having(length: 4);
  bool BAND, GEN, NOTRAN, QRS, SYM, TRAN, TRI;
  String C1, DIAG;
  String C2;
  int J, MB, NX;

  // Test the input parameters.

  INFO.value = 0;
  C1 = PATH.substring(0, 1);
  C2 = PATH.substring(1, 3);
  TRAN = lsame(TRANS, 'T') || lsame(TRANS, 'C');
  NOTRAN = !TRAN;
  GEN = lsame(PATH.substring(1, 2), 'G');
  QRS = lsame(PATH.substring(1, 2), 'Q') || lsame(PATH.substring(2, 3), 'Q');
  SYM = lsame(PATH.substring(1, 2), 'P') ||
      lsame(PATH.substring(1, 2), 'S') ||
      lsame(PATH.substring(1, 2), 'H');
  TRI = lsame(PATH.substring(1, 2), 'T');
  BAND = lsame(PATH.substring(2, 3), 'B');
  if (!lsame(C1, 'Zomplex precision')) {
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
    xerbla('ZLARHS', -INFO.value);
    return;
  }

  // Initialize X to NRHS random vectors unless XTYPE = 'C'.

  if (TRAN) {
    NX = M;
    MB = N;
  } else {
    NX = N;
    MB = M;
  }
  if (!lsame(XTYPE, 'C')) {
    for (J = 1; J <= NRHS; J++) {
      zlarnv(2, ISEED, N, X(1, J).asArray());
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

    zgemm(TRANS, 'N', MB, NRHS, NX, Complex.one, A, LDA, X, LDX, Complex.zero,
        B, LDB);
  } else if (lsamen(2, C2, 'PO') || lsamen(2, C2, 'HE')) {
    // Hermitian matrix, 2-D storage

    zhemm('Left', UPLO, N, NRHS, Complex.one, A, LDA, X, LDX, Complex.zero, B,
        LDB);
  } else if (lsamen(2, C2, 'SY')) {
    // Symmetric matrix, 2-D storage

    zsymm('Left', UPLO, N, NRHS, Complex.one, A, LDA, X, LDX, Complex.zero, B,
        LDB);
  } else if (lsamen(2, C2, 'GB')) {
    // General matrix, band storage

    for (J = 1; J <= NRHS; J++) {
      zgbmv(TRANS, M, N, KL, KU, Complex.one, A, LDA, X(1, J).asArray(), 1,
          Complex.zero, B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'PB') || lsamen(2, C2, 'HB')) {
    // Hermitian matrix, band storage

    for (J = 1; J <= NRHS; J++) {
      zhbmv(UPLO, N, KL, Complex.one, A, LDA, X(1, J).asArray(), 1,
          Complex.zero, B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'SB')) {
    // Symmetric matrix, band storage

    for (J = 1; J <= NRHS; J++) {
      zsbmv(UPLO, N, KL, Complex.one, A, LDA, X(1, J).asArray(), 1,
          Complex.zero, B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'PP') || lsamen(2, C2, 'HP')) {
    // Hermitian matrix, packed storage

    for (J = 1; J <= NRHS; J++) {
      zhpmv(UPLO, N, Complex.one, A.asArray(), X(1, J).asArray(), 1,
          Complex.zero, B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'SP')) {
    // Symmetric matrix, packed storage

    for (J = 1; J <= NRHS; J++) {
      zspmv(UPLO, N, Complex.one, A.asArray(), X(1, J).asArray(), 1,
          Complex.zero, B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'TR')) {
    // Triangular matrix.  Note that for triangular matrices,
    //    KU = 1 => non-unit triangular
    //    KU = 2 => unit triangular

    zlacpy('Full', N, NRHS, X, LDX, B, LDB);
    if (KU == 2) {
      DIAG = 'U';
    } else {
      DIAG = 'N';
    }
    ztrmm('Left', UPLO, TRANS, DIAG, N, NRHS, Complex.one, A, LDA, B, LDB);
  } else if (lsamen(2, C2, 'TP')) {
    // Triangular matrix, packed storage

    zlacpy('Full', N, NRHS, X, LDX, B, LDB);
    if (KU == 2) {
      DIAG = 'U';
    } else {
      DIAG = 'N';
    }
    for (J = 1; J <= NRHS; J++) {
      ztpmv(UPLO, TRANS, DIAG, N, A.asArray(), B(1, J).asArray(), 1);
    }
  } else if (lsamen(2, C2, 'TB')) {
    // Triangular matrix, banded storage

    zlacpy('Full', N, NRHS, X, LDX, B, LDB);
    if (KU == 2) {
      DIAG = 'U';
    } else {
      DIAG = 'N';
    }
    for (J = 1; J <= NRHS; J++) {
      ztbmv(UPLO, TRANS, DIAG, N, KL, A, LDA, B(1, J).asArray(), 1);
    }
  } else {
    // If none of the above, set INFO = -1 and return;

    INFO.value = -1;
    xerbla('ZLARHS', -INFO.value);
  }
}
