import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dsptrs.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dspcon(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<int> IPIV_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  int I, IP;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (ANORM < ZERO) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DSPCON', -INFO.value);
    return;
  }

  // Quick return if possible

  RCOND.value = ZERO;
  if (N == 0) {
    RCOND.value = ONE;
    return;
  } else if (ANORM <= ZERO) {
    return;
  }

  // Check that the diagonal matrix D is nonsingular.

  if (UPPER) {
    // Upper triangular storage: examine D from bottom to top

    IP = N * (N + 1) ~/ 2;
    for (I = N; I >= 1; I--) {
      if (IPIV[I] > 0 && AP[IP] == ZERO) return;
      IP -= I;
    }
  } else {
    // Lower triangular storage: examine D from top to bottom.

    IP = 1;
    for (I = 1; I <= N; I++) {
      if (IPIV[I] > 0 && AP[IP] == ZERO) return;
      IP += N - I + 1;
    }
  }

  // Estimate the 1-norm of the inverse.

  KASE.value = 0;
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    // Multiply by inv(L*D*L**T) or inv(U*D*U**T).

    dsptrs(UPLO, N, 1, AP, IPIV, WORK.asMatrix(N), N, INFO);
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
