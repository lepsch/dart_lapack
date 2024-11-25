// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zggsvp3.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/ztgsja.dart';

void zggsvd3(
  final String JOBU,
  final String JOBV,
  final String JOBQ,
  final int M,
  final int N,
  final int P,
  final Box<int> K,
  final Box<int> L,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Array<double> ALPHA_,
  final Array<double> BETA_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final Q = Q_.having(ld: LDQ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  bool WANTQ, WANTU, WANTV, LQUERY;
  int I, IBND, ISUB, J, LWKOPT;
  double ANORM, BNORM, SMAX, TEMP, TOLA = 0, TOLB = 0, ULP, UNFL;
  final NCYCLE = Box(0);

  // Decode and test the input parameters

  WANTU = lsame(JOBU, 'U');
  WANTV = lsame(JOBV, 'V');
  WANTQ = lsame(JOBQ, 'Q');
  LQUERY = (LWORK == -1);
  LWKOPT = 1;

  // Test the input arguments

  INFO.value = 0;
  if (!(WANTU || lsame(JOBU, 'N'))) {
    INFO.value = -1;
  } else if (!(WANTV || lsame(JOBV, 'N'))) {
    INFO.value = -2;
  } else if (!(WANTQ || lsame(JOBQ, 'N'))) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (P < 0) {
    INFO.value = -6;
  } else if (LDA < max(1, M)) {
    INFO.value = -10;
  } else if (LDB < max(1, P)) {
    INFO.value = -12;
  } else if (LDU < 1 || (WANTU && LDU < M)) {
    INFO.value = -16;
  } else if (LDV < 1 || (WANTV && LDV < P)) {
    INFO.value = -18;
  } else if (LDQ < 1 || (WANTQ && LDQ < N)) {
    INFO.value = -20;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -24;
  }

  // Compute workspace

  if (INFO.value == 0) {
    zggsvp3(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU,
        V, LDV, Q, LDQ, IWORK, RWORK, WORK, WORK, -1, INFO);
    LWKOPT = N + WORK[1].toInt();
    LWKOPT = max(2 * N, LWKOPT);
    LWKOPT = max(1, LWKOPT);
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZGGSVD3', -INFO.value);
    return;
  }
  if (LQUERY) {
    return;
  }

  // Compute the Frobenius norm of matrices A and B

  ANORM = zlange('1', M, N, A, LDA, RWORK);
  BNORM = zlange('1', P, N, B, LDB, RWORK);

  // Get machine precision and set up threshold for determining
  // the effective numerical rank of the matrices A and B.

  ULP = dlamch('Precision');
  UNFL = dlamch('Safe Minimum');
  TOLA = max(M, N) * max(ANORM, UNFL) * ULP;
  TOLB = max(P, N) * max(BNORM, UNFL) * ULP;

  zggsvp3(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU,
      V, LDV, Q, LDQ, IWORK, RWORK, WORK, WORK(N + 1), LWORK - N, INFO);

  // Compute the GSVD of two upper "triangular" matrices

  ztgsja(JOBU, JOBV, JOBQ, M, P, N, K.value, L.value, A, LDA, B, LDB, TOLA,
      TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO);

  // Sort the singular values and store the pivot indices in IWORK
  // Copy ALPHA to RWORK, then sort ALPHA in RWORK

  dcopy(N, ALPHA, 1, RWORK, 1);
  IBND = min(L.value, M - K.value);
  for (I = 1; I <= IBND; I++) {
    // Scan for largest ALPHA(K+I)

    ISUB = I;
    SMAX = RWORK[K.value + I];
    for (J = I + 1; J <= IBND; J++) {
      TEMP = RWORK[K.value + J];
      if (TEMP > SMAX) {
        ISUB = J;
        SMAX = TEMP;
      }
    }
    if (ISUB != I) {
      RWORK[K.value + ISUB] = RWORK[K.value + I];
      RWORK[K.value + I] = SMAX;
      IWORK[K.value + I] = K.value + ISUB;
    } else {
      IWORK[K.value + I] = K.value + I;
    }
  }

  WORK[1] = LWKOPT.toComplex();
}
