import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggqrf.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/dorgrq.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dgqrts(
  final int N,
  final int M,
  final int P,
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

  dlacpy('Full', N, M, A, LDA, AF, LDA);
  dlacpy('Full', N, P, B, LDB, BF, LDB);

  ANORM = max(dlange('1', N, M, A, LDA, RWORK), UNFL);
  BNORM = max(dlange('1', N, P, B, LDB, RWORK), UNFL);

  // Factorize the matrices A and B in the arrays AF and BF.

  dggqrf(N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO);

  // Generate the N-by-N matrix Q

  dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA);
  dlacpy('Lower', N - 1, M, AF(2, 1), LDA, Q(2, 1), LDA);
  dorgqr(N, N, min(N, M), Q, LDA, TAUA, WORK, LWORK, INFO);

  // Generate the P-by-P matrix Z

  dlaset('Full', P, P, ROGUE, ROGUE, Z, LDB);
  if (N <= P) {
    if (N > 0 && N < P) dlacpy('Full', N, P - N, BF, LDB, Z(P - N + 1, 1), LDB);
    if (N > 1) {
      dlacpy('Lower', N - 1, N - 1, BF(2, P - N + 1), LDB,
          Z(P - N + 2, P - N + 1), LDB);
    }
  } else {
    if (P > 1) {
      dlacpy('Lower', P - 1, P - 1, BF(N - P + 2, 1), LDB, Z(2, 1), LDB);
    }
  }
  dorgrq(P, P, min(N, P), Z, LDB, TAUB, WORK, LWORK, INFO);

  // Copy R

  dlaset('Full', N, M, ZERO, ZERO, R, LDA);
  dlacpy('Upper', N, M, AF, LDA, R, LDA);

  // Copy T

  dlaset('Full', N, P, ZERO, ZERO, T, LDB);
  if (N <= P) {
    dlacpy('Upper', N, N, BF(1, P - N + 1), LDB, T(1, P - N + 1), LDB);
  } else {
    dlacpy('Full', N - P, P, BF, LDB, T, LDB);
    dlacpy('Upper', P, P, BF(N - P + 1, 1), LDB, T(N - P + 1, 1), LDB);
  }

  // Compute R - Q'*A

  dgemm(
      'Transpose', 'No transpose', N, M, N, -ONE, Q, LDA, A, LDA, ONE, R, LDA);

  // Compute norm( R - Q'*A ) / ( max(M,N)*norm(A)*ULP ) .

  RESID = dlange('1', N, M, R, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / max(1, max(M, N))) / ANORM) / ULP;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute T*Z - Q'*B

  dgemm('No Transpose', 'No transpose', N, P, P, ONE, T, LDB, Z, LDB, ZERO, BWK,
      LDB);
  dgemm('Transpose', 'No transpose', N, P, N, -ONE, Q, LDA, B, LDB, ONE, BWK,
      LDB);

  // Compute norm( T*Z - Q'*B ) / ( max(P,N)*norm(A)*ULP ) .

  RESID = dlange('1', N, P, BWK, LDB, RWORK);
  if (BNORM > ZERO) {
    RESULT[2] = ((RESID / max(1, max(P, N))) / BNORM) / ULP;
  } else {
    RESULT[2] = ZERO;
  }

  // Compute I - Q'*Q

  dlaset('Full', N, N, ZERO, ONE, R, LDA);
  dsyrk('Upper', 'Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA);

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
