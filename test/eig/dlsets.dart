import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgglse.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/matrix.dart';

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
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDA);
  final B = B_.dim(LDB);
  final BF = BF_.dim(LDB);
  final C = C_.dim();
  final CF = CF_.dim();
  final D = D_.dim();
  final DF = DF_.dim();
  final X = X_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
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
