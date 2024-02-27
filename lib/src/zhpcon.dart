import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhptrs.dart';
import 'package:lapack/src/zlacn2.dart';

void zhpcon(
  final String UPLO,
  final int N,
  final Array<Complex> AP,
  final Array<int> IPIV_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final IPIV = IPIV_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
    xerbla('ZHPCON', -INFO.value);
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
      // 10
      if (IPIV[I] > 0 && AP[IP] == Complex.zero) return;
      IP = IP - I;
    } // 10
  } else {
    // Lower triangular storage: examine D from top to bottom.

    IP = 1;
    for (I = 1; I <= N; I++) {
      // 20
      if (IPIV[I] > 0 && AP[IP] == Complex.zero) return;
      IP = IP + N - I + 1;
    } // 20
  }

  // Estimate the 1-norm of the inverse.

  KASE.value = 0;
  while (true) {
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    // Multiply by inv(L*D*L**H) or inv(U*D*U**H).

    zhptrs(UPLO, N, 1, AP, IPIV, WORK.asMatrix(), N, INFO);
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
