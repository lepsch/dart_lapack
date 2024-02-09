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
  final Matrix<double> A,
  final Matrix<double> AF,
  final int LDA,
  final Matrix<double> B,
  final Matrix<double> BF,
  final int LDB,
  final Array<double> C,
  final Array<double> CF,
  final Array<double> D,
  final Array<double> DF,
  final Array<double> X,
  final Array<double> WORK,
  final int LWORK,
  final Array<double> RWORK,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
  dget02('No transpose', M, N, 1, A, LDA, X, N, CF, M, RWORK, RESULT[1]);

  // Compute result[2] = norm( B*x - d ) / norm(B)*norm(X)*EPS

  dget02('No transpose', P, N, 1, B, LDB, X, N, DF, P, RWORK, RESULT[2]);
}
