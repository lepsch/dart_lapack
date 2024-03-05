import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dsymm.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dsgt01(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int M,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> D_,
  final Array<double> WORK_,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Z = Z_.having(ld: LDZ);
  final D = D_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I;
  double ANORM, ULP;

  RESULT[1] = ZERO;
  if (N <= 0) return;

  ULP = dlamch('Epsilon');

  // Compute product of 1-norms of A and Z.

  ANORM = dlansy('1', UPLO, N, A, LDA, WORK) * dlange('1', N, M, Z, LDZ, WORK);
  if (ANORM == ZERO) ANORM = ONE;

  if (ITYPE == 1) {
    // Norm of AZ - BZD

    dsymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK.asMatrix(N), N);
    for (I = 1; I <= M; I++) {
      dscal(N, D[I], Z(1, I).asArray(), 1);
    }
    dsymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, -ONE, WORK.asMatrix(N), N);

    RESULT[1] =
        (dlange('1', N, M, WORK.asMatrix(N), N, WORK) / ANORM) / (N * ULP);
  } else if (ITYPE == 2) {
    // Norm of ABZ - ZD

    dsymm('Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO, WORK.asMatrix(N), N);
    for (I = 1; I <= M; I++) {
      dscal(N, D[I], Z(1, I).asArray(), 1);
    }
    dsymm('Left', UPLO, N, M, ONE, A, LDA, WORK.asMatrix(N), N, -ONE, Z, LDZ);

    RESULT[1] = (dlange('1', N, M, Z, LDZ, WORK) / ANORM) / (N * ULP);
  } else if (ITYPE == 3) {
    // Norm of BAZ - ZD

    dsymm('Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO, WORK.asMatrix(N), N);
    for (I = 1; I <= M; I++) {
      dscal(N, D[I], Z(1, I).asArray(), 1);
    }
    dsymm('Left', UPLO, N, M, ONE, B, LDB, WORK.asMatrix(N), N, -ONE, Z, LDZ);

    RESULT[1] = (dlange('1', N, M, Z, LDZ, WORK) / ANORM) / (N * ULP);
  }
}
