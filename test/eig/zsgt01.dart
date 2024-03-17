import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zhemm.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhe.dart';

void zsgt01(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int M,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<double> D_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final D = D_.having();
  final RESULT = RESULT_.having();
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int I;
  double ANORM, ULP;

  RESULT[1] = ZERO;
  if (N <= 0) return;

  ULP = dlamch('Epsilon');

  // Compute product of 1-norms of A and Z.

  ANORM =
      zlanhe('1', UPLO, N, A, LDA, RWORK) * zlange('1', N, M, Z, LDZ, RWORK);
  if (ANORM == ZERO) ANORM = ONE;

  if (ITYPE == 1) {
    // Norm of AZ - BZD

    zhemm('Left', UPLO, N, M, Complex.one, A, LDA, Z, LDZ, Complex.zero,
        WORK.asMatrix(), N);
    for (I = 1; I <= M; I++) {
      zdscal(N, D[I], Z(1, I).asArray(), 1);
    }
    zhemm('Left', UPLO, N, M, Complex.one, B, LDB, Z, LDZ, -Complex.one,
        WORK.asMatrix(), N);

    RESULT[1] =
        (zlange('1', N, M, WORK.asMatrix(), N, RWORK) / ANORM) / (N * ULP);
  } else if (ITYPE == 2) {
    // Norm of ABZ - ZD

    zhemm('Left', UPLO, N, M, Complex.one, B, LDB, Z, LDZ, Complex.zero,
        WORK.asMatrix(), N);
    for (I = 1; I <= M; I++) {
      zdscal(N, D[I], Z(1, I).asArray(), 1);
    }
    zhemm('Left', UPLO, N, M, Complex.one, A, LDA, WORK.asMatrix(), N,
        -Complex.one, Z, LDZ);

    RESULT[1] = (zlange('1', N, M, Z, LDZ, RWORK) / ANORM) / (N * ULP);
  } else if (ITYPE == 3) {
    // Norm of BAZ - ZD

    zhemm('Left', UPLO, N, M, Complex.one, A, LDA, Z, LDZ, Complex.zero,
        WORK.asMatrix(), N);
    for (I = 1; I <= M; I++) {
      zdscal(N, D[I], Z(1, I).asArray(), 1);
    }
    zhemm('Left', UPLO, N, M, Complex.one, B, LDB, WORK.asMatrix(), N,
        -Complex.one, Z, LDZ);

    RESULT[1] = (zlange('1', N, M, Z, LDZ, RWORK) / ANORM) / (N * ULP);
  }
}
