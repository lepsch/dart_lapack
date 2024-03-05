import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgglse.dart';
import 'package:lapack/src/zlacpy.dart';

import 'zget02.dart';

void zlsets(
  final int M,
  final int P,
  final int N,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> BF_,
  final int LDB,
  final Array<Complex> C_,
  final Array<Complex> CF_,
  final Array<Complex> D_,
  final Array<Complex> DF_,
  final Array<Complex> X_,
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
  final B = B_.having(ld: LDB);
  final BF = BF_.having(ld: LDB);
  final C = C_.having();
  final CF = CF_.having();
  final D = D_.having();
  final DF = DF_.having();
  final X = X_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having(length: 2);

  final INFO = Box(0);
  // ..
  // .. External Subroutines ..
  // EXTERNAL ZCOPY, ZGET02, ZGGLSE, ZLACPY

  // Copy the matrices A and B to the arrays AF and BF,
  // and the vectors C and D to the arrays CF and DF,

  zlacpy('Full', M, N, A, LDA, AF, LDA);
  zlacpy('Full', P, N, B, LDB, BF, LDB);
  zcopy(M, C, 1, CF, 1);
  zcopy(P, D, 1, DF, 1);

  // Solve LSE problem

  zgglse(M, N, P, AF, LDA, BF, LDB, CF, DF, X, WORK, LWORK, INFO);

  // Test the residual for the solution of LSE

  // Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS

  zcopy(M, C, 1, CF, 1);
  zcopy(P, D, 1, DF, 1);
  zget02('No transpose', M, N, 1, A, LDA, X.asMatrix(), N, CF.asMatrix(), M,
      RWORK, RESULT(1));

  // Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS

  zget02('No transpose', P, N, 1, B, LDB, X.asMatrix(), N, DF.asMatrix(), P,
      RWORK, RESULT(2));
}
