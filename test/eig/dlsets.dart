// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dgglse.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/matrix.dart';

import 'dget02.dart';

void dlsets(
  final int M,
  final int P,
  final int N,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> BF_,
  final int LDB,
  final Array<double> C_,
  final Array<double> CF_,
  final Array<double> D_,
  final Array<double> DF_,
  final Array<double> X_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT,
) {
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final BF = BF_.having(ld: LDB);
  final C = C_.having();
  final CF = CF_.having();
  final D = D_.having();
  final DF = DF_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final INFO = Box(0);

  // Copy the matrices A and B to the arrays AF and BF,
  // and the vectors C and D to the arrays CF and DF,

  dlacpy('Full', M, N, A, LDA, AF, LDA);
  dlacpy('Full', P, N, B, LDB, BF, LDB);
  dcopy(M, C, 1, CF, 1);
  dcopy(P, D, 1, DF, 1);

  // Solve LSE problem

  dgglse(M, N, P, AF, LDA, BF, LDB, CF, DF, X, WORK, LWORK, INFO);

  // Test the residual for the solution of LSE

  // Compute RESULT[1] = norm( A*x - c ) / norm(A)*norm(X)*EPS

  dcopy(M, C, 1, CF, 1);
  dcopy(P, D, 1, DF, 1);
  dget02('No transpose', M, N, 1, A, LDA, X.asMatrix(N), N, CF.asMatrix(M), M,
      RWORK, RESULT.box(1));

  // Compute result[2] = norm( B*x - d ) / norm(B)*norm(X)*EPS

  dget02('No transpose', P, N, 1, B, LDB, X.asMatrix(N), N, DF.asMatrix(P), P,
      RWORK, RESULT.box(2));
}
