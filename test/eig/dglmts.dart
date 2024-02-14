import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggglm.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dglmts(
  final int N,
  final int M,
  final int P,
  final Matrix<double> A_,
  final Matrix<double> AF_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> BF_,
  final int LDB,
  final Array<double> D_,
  final Array<double> DF_,
  final Array<double> X_,
  final Array<double> U_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDA);
  final B = B_.dim(LDB);
  final BF = BF_.dim(LDB);
  final D = D_.dim();
  final DF = DF_.dim();
  final X = X_.dim();
  final U = U_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  final INFO = Box(0);
  double ANORM, BNORM, DNORM, EPS, UNFL, XNORM, YNORM;

  EPS = dlamch('Epsilon');
  UNFL = dlamch('Safe minimum');
  ANORM = max(dlange('1', N, M, A, LDA, RWORK), UNFL);
  BNORM = max(dlange('1', N, P, B, LDB, RWORK), UNFL);

  // Copy the matrices A and B to the arrays AF and BF,
  // and the vector D the array DF.

  dlacpy('Full', N, M, A, LDA, AF, LDA);
  dlacpy('Full', N, P, B, LDB, BF, LDB);
  dcopy(N, D, 1, DF, 1);

  // Solve GLM problem

  dggglm(N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK, INFO);

  // Test the residual for the solution of LSE
  //
  //                   norm( d - A*x - B*u )
  //   RESULT = -----------------------------------------
  //            (norm(A)+norm(B))*(norm(x)+norm(u))*EPS

  dcopy(N, D, 1, DF, 1);
  dgemv('No transpose', N, M, -ONE, A, LDA, X, 1, ONE, DF, 1);

  dgemv('No transpose', N, P, -ONE, B, LDB, U, 1, ONE, DF, 1);

  DNORM = dasum(N, DF, 1);
  XNORM = dasum(M, X, 1) + dasum(P, U, 1);
  YNORM = ANORM + BNORM;

  if (XNORM <= ZERO) {
    RESULT.value = ZERO;
  } else {
    RESULT.value = ((DNORM / YNORM) / XNORM) / EPS;
  }
}
