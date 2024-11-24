// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zggrqf.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungqr.dart';
import 'package:dart_lapack/src/zungrq.dart';

void zgrqts(
  final int M,
  final int P,
  final int N,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final Matrix<Complex> Q_,
  final Matrix<Complex> R_,
  final int LDA,
  final Array<Complex> TAUA_,
  final Matrix<Complex> B_,
  final Matrix<Complex> BF_,
  final Matrix<Complex> Z_,
  final Matrix<Complex> T_,
  final Matrix<Complex> BWK_,
  final int LDB,
  final Array<Complex> TAUB_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final Q = Q_.having(ld: LDA);
  final R = R_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final BF = BF_.having(ld: LDB);
  final Z = Z_.having(ld: LDB);
  final T = T_.having(ld: LDB);
  final BWK = BWK_.having(ld: LDB);
  final TAUA = TAUA_.having();
  final TAUB = TAUB_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having(length: 4);
  const ZERO = 0.0, ONE = 1.0;
  const CROGUE = Complex(-1.0e+10, 0.0);
  final INFO = Box(0);
  double ANORM, BNORM, RESID, ULP, UNFL;

  ULP = dlamch('Precision');
  UNFL = dlamch('Safe minimum');

  // Copy the matrix A to the array AF.

  zlacpy('Full', M, N, A, LDA, AF, LDA);
  zlacpy('Full', P, N, B, LDB, BF, LDB);

  ANORM = max(zlange('1', M, N, A, LDA, RWORK), UNFL);
  BNORM = max(zlange('1', P, N, B, LDB, RWORK), UNFL);

  // Factorize the matrices A and B in the arrays AF and BF.

  zggrqf(M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO);

  // Generate the N-by-N matrix Q

  zlaset('Full', N, N, CROGUE, CROGUE, Q, LDA);
  if (M <= N) {
    if (M > 0 && M < N) zlacpy('Full', M, N - M, AF, LDA, Q(N - M + 1, 1), LDA);
    if (M > 1) {
      zlacpy('Lower', M - 1, M - 1, AF(2, N - M + 1), LDA,
          Q(N - M + 2, N - M + 1), LDA);
    }
  } else {
    if (N > 1) {
      zlacpy('Lower', N - 1, N - 1, AF(M - N + 2, 1), LDA, Q(2, 1), LDA);
    }
  }
  zungrq(N, N, min(M, N), Q, LDA, TAUA, WORK, LWORK, INFO);

  // Generate the P-by-P matrix Z

  zlaset('Full', P, P, CROGUE, CROGUE, Z, LDB);
  if (P > 1) zlacpy('Lower', P - 1, N, BF(2, 1), LDB, Z(2, 1), LDB);
  zungqr(P, P, min(P, N), Z, LDB, TAUB, WORK, LWORK, INFO);

  // Copy R

  zlaset('Full', M, N, Complex.zero, Complex.zero, R, LDA);
  if (M <= N) {
    zlacpy('Upper', M, M, AF(1, N - M + 1), LDA, R(1, N - M + 1), LDA);
  } else {
    zlacpy('Full', M - N, N, AF, LDA, R, LDA);
    zlacpy('Upper', N, N, AF(M - N + 1, 1), LDA, R(M - N + 1, 1), LDA);
  }

  // Copy T

  zlaset('Full', P, N, Complex.zero, Complex.zero, T, LDB);
  zlacpy('Upper', P, N, BF, LDB, T, LDB);

  // Compute R - A*Q'

  zgemm('No transpose', 'Conjugate transpose', M, N, N, -Complex.one, A, LDA, Q,
      LDA, Complex.one, R, LDA);

  // Compute norm( R - A*Q' ) / ( max(M,N)*norm(A)*ULP ) .

  RESID = zlange('1', M, N, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, max(M, N)) ) / ANORM) / ULP;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute T*Q - Z'*B

  zgemm('Conjugate transpose', 'No transpose', P, N, P, Complex.one, Z, LDB, B,
      LDB, Complex.zero, BWK, LDB);
  zgemm('No transpose', 'No transpose', P, N, N, Complex.one, T, LDB, Q, LDA,
      -Complex.one, BWK, LDB);

  // Compute norm( T*Q - Z'*B ) / ( max(P,N)*norm(A)*ULP ) .

  RESID = zlange('1', P, N, BWK, LDB, RWORK);
  if (BNORM > ZERO) {
    RESULT[2] = ((RESID / max(1, max(P, M)) ) / BNORM) / ULP;
  } else {
    RESULT[2] = ZERO;
  }

  // Compute I - Q*Q'

  zlaset('Full', N, N, Complex.zero, Complex.one, R, LDA);
  zherk('Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q'*Q ) / ( N * ULP ) .

  RESID = zlanhe('1', 'Upper', N, R, LDA, RWORK);
  RESULT[3] = (RESID / max(1, N) ) / ULP;

  // Compute I - Z'*Z

  zlaset('Full', P, P, Complex.zero, Complex.one, T, LDB);
  zherk('Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB);

  // Compute norm( I - Z'*Z ) / ( P*ULP ) .

  RESID = zlanhe('1', 'Upper', P, T, LDB, RWORK);
  RESULT[4] = (RESID / max(1, P) ) / ULP;
}
