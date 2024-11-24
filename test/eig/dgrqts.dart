// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggrqf.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/dorgrq.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dgrqts(
  final int M,
  final int P,
  final int N,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final Matrix<double> Q_,
  final Matrix<double> R_,
  final int LDA,
  final Array<double> TAUA_,
  final Matrix<double> B_,
  final Matrix<double> BF_,
  final Matrix<double> Z_,
  final Matrix<double> T_,
  final Matrix<double> BWK_,
  final int LDB,
  final Array<double> TAUB_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final R = R_.having(ld: LDA);
  final TAUA = TAUA_.having();
  final B = B_.having(ld: LDB);
  final BF = BF_.having(ld: LDB);
  final Z = Z_.having(ld: LDB);
  final T = T_.having(ld: LDB);
  final BWK = BWK_.having(ld: LDB);
  final TAUB = TAUB_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = -1.0e+10;
  final INFO = Box(0);
  double ANORM, BNORM, RESID, ULP, UNFL;

  ULP = dlamch('Precision');
  UNFL = dlamch('Safe minimum');

  // Copy the matrix A to the array AF.

  dlacpy('Full', M, N, A, LDA, AF, LDA);
  dlacpy('Full', P, N, B, LDB, BF, LDB);

  ANORM = max(dlange('1', M, N, A, LDA, RWORK), UNFL);
  BNORM = max(dlange('1', P, N, B, LDB, RWORK), UNFL);

  // Factorize the matrices A and B in the arrays AF and BF.

  dggrqf(M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO);

  // Generate the N-by-N matrix Q

  dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  if (M <= N) {
    if (M > 0 && M < N) dlacpy('Full', M, N - M, AF, LDA, Q(N - M + 1, 1), LDA);
    if (M > 1) {
      dlacpy('Lower', M - 1, M - 1, AF(2, N - M + 1), LDA,
          Q(N - M + 2, N - M + 1), LDA);
    }
  } else {
    if (N > 1) {
      dlacpy('Lower', N - 1, N - 1, AF(M - N + 2, 1), LDA, Q(2, 1), LDA);
    }
  }
  dorgrq(N, N, min(M, N), Q, LDA, TAUA, WORK, LWORK, INFO);

  // Generate the P-by-P matrix Z

  dlaset('Full', P, P, ROGUE, ROGUE, Z, LDB);
  if (P > 1) dlacpy('Lower', P - 1, N, BF(2, 1), LDB, Z(2, 1), LDB);
  dorgqr(P, P, min(P, N), Z, LDB, TAUB, WORK, LWORK, INFO);

  // Copy R

  dlaset('Full', M, N, ZERO, ZERO, R, LDA);
  if (M <= N) {
    dlacpy('Upper', M, M, AF(1, N - M + 1), LDA, R(1, N - M + 1), LDA);
  } else {
    dlacpy('Full', M - N, N, AF, LDA, R, LDA);
    dlacpy('Upper', N, N, AF(M - N + 1, 1), LDA, R(M - N + 1, 1), LDA);
  }

  // Copy T

  dlaset('Full', P, N, ZERO, ZERO, T, LDB);
  dlacpy('Upper', P, N, BF, LDB, T, LDB);

  // Compute R - A*Q'

  dgemm(
      'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA);

  // Compute norm( R - A*Q' ) / ( max(M,N)*norm(A)*ULP ) .

  RESID = dlange('1', M, N, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, max(M, N))) / ANORM) / ULP;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute T*Q - Z'*B

  dgemm('Transpose', 'No transpose', P, N, P, ONE, Z, LDB, B, LDB, ZERO, BWK,
      LDB);
  dgemm('No transpose', 'No transpose', P, N, N, ONE, T, LDB, Q, LDA, -ONE, BWK,
      LDB);

  // Compute norm( T*Q - Z'*B ) / ( max(P,N)*norm(A)*ULP ) .

  RESID = dlange('1', P, N, BWK, LDB, RWORK);
  if (BNORM > ZERO) {
    RESULT[2] = ((RESID / max(1, max(P, M))) / BNORM) / ULP;
  } else {
    RESULT[2] = ZERO;
  }

  // Compute I - Q*Q'

  dlaset('Full', N, N, ZERO, ONE, R, LDA);
  dsyrk('Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q'*Q ) / ( N * ULP ) .

  RESID = dlansy('1', 'Upper', N, R, LDA, RWORK);
  RESULT[3] = (RESID / max(1, N)) / ULP;

  // Compute I - Z'*Z

  dlaset('Full', P, P, ZERO, ONE, T, LDB);
  dsyrk('Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB);

  // Compute norm( I - Z'*Z ) / ( P*ULP ) .

  RESID = dlansy('1', 'Upper', P, T, LDB, RWORK);
  RESULT[4] = (RESID / max(1, P)) / ULP;
}
