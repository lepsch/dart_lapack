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
import 'package:dart_lapack/src/zggqrf.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlange.dart';
import 'package:dart_lapack/src/zlanhe.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zungqr.dart';
import 'package:dart_lapack/src/zungrq.dart';

void zgqrts(
  final int N,
  final int M,
  final int P,
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
  final RESULT = RESULT_.having(length: 4);
  const ZERO = 0.0, ONE = 1.0;
  const CROGUE = Complex(-1.0e+10, 0.0);
  final INFO = Box(0);
  double ANORM, BNORM, RESID, ULP, UNFL;

  ULP = dlamch('Precision');
  UNFL = dlamch('Safe minimum');

  // Copy the matrix A to the array AF.

  zlacpy('Full', N, M, A, LDA, AF, LDA);
  zlacpy('Full', N, P, B, LDB, BF, LDB);

  ANORM = max(zlange('1', N, M, A, LDA, RWORK), UNFL);
  BNORM = max(zlange('1', N, P, B, LDB, RWORK), UNFL);

  // Factorize the matrices A and B in the arrays AF and BF.

  zggqrf(N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO);

  // Generate the N-by-N matrix Q

  zlaset('Full', N, N, CROGUE, CROGUE, Q, LDA);
  zlacpy('Lower', N - 1, M, AF(2, 1), LDA, Q(2, 1), LDA);
  zungqr(N, N, min(N, M), Q, LDA, TAUA, WORK, LWORK, INFO);

  // Generate the P-by-P matrix Z

  zlaset('Full', P, P, CROGUE, CROGUE, Z, LDB);
  if (N <= P) {
    if (N > 0 && N < P) zlacpy('Full', N, P - N, BF, LDB, Z(P - N + 1, 1), LDB);
    if (N > 1) {
      zlacpy('Lower', N - 1, N - 1, BF(2, P - N + 1), LDB,
          Z(P - N + 2, P - N + 1), LDB);
    }
  } else {
    if (P > 1) {
      zlacpy('Lower', P - 1, P - 1, BF(N - P + 2, 1), LDB, Z(2, 1), LDB);
    }
  }
  zungrq(P, P, min(N, P), Z, LDB, TAUB, WORK, LWORK, INFO);

  // Copy R

  zlaset('Full', N, M, Complex.zero, Complex.zero, R, LDA);
  zlacpy('Upper', N, M, AF, LDA, R, LDA);

  // Copy T

  zlaset('Full', N, P, Complex.zero, Complex.zero, T, LDB);
  if (N <= P) {
    zlacpy('Upper', N, N, BF(1, P - N + 1), LDB, T(1, P - N + 1), LDB);
  } else {
    zlacpy('Full', N - P, P, BF, LDB, T, LDB);
    zlacpy('Upper', P, P, BF(N - P + 1, 1), LDB, T(N - P + 1, 1), LDB);
  }

  // Compute R - Q'*A

  zgemm('Conjugate transpose', 'No transpose', N, M, N, -Complex.one, Q, LDA, A,
      LDA, Complex.one, R, LDA);

  // Compute norm( R - Q'*A ) / ( max(M,N)*norm(A)*ULP ) .

  RESID = zlange('1', N, M, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, max(M, N))) / ANORM) / ULP;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute T*Z - Q'*B

  zgemm('No Transpose', 'No transpose', N, P, P, Complex.one, T, LDB, Z, LDB,
      Complex.zero, BWK, LDB);
  zgemm('Conjugate transpose', 'No transpose', N, P, N, -Complex.one, Q, LDA, B,
      LDB, Complex.one, BWK, LDB);

  // Compute norm( T*Z - Q'*B ) / ( max(P,N)*norm(A)*ULP ) .

  RESID = zlange('1', N, P, BWK, LDB, RWORK);
  if (BNORM > ZERO) {
    RESULT[2] = ((RESID / max(1, max(P, N))) / BNORM) / ULP;
  } else {
    RESULT[2] = ZERO;
  }

  // Compute I - Q'*Q

  zlaset('Full', N, N, Complex.zero, Complex.one, R, LDA);
  zherk('Upper', 'Conjugate transpose', N, N, -ONE, Q, LDA, ONE, R, LDA);

  // Compute norm( I - Q'*Q ) / ( N * ULP ) .

  RESID = zlanhe('1', 'Upper', N, R, LDA, RWORK);
  RESULT[3] = (RESID / max(1, N)) / ULP;

  // Compute I - Z'*Z

  zlaset('Full', P, P, Complex.zero, Complex.one, T, LDB);
  zherk('Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB);

  // Compute norm( I - Z'*Z ) / ( P*ULP ) .

  RESID = zlanhe('1', 'Upper', P, T, LDB, RWORK);
  RESULT[4] = (RESID / max(1, P)) / ULP;
}
