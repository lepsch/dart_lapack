import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgttrs.dart';
import 'package:lapack/src/zlacn2.dart';

void zgtcon(
  final String NORM,
  final int N,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Array<Complex> DU2_,
  final Array<int> IPIV_,
  final double ANORM,
  final Box<double> RCOND,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final IPIV = IPIV_.having();
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  final WORK = WORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  bool ONENRM;
  int I, KASE1;
  final KASE = Box(0);
  final AINVNM = Box(0.0);
  final ISAVE = Array<int>(3);

  // Test the input arguments.

  INFO.value = 0;
  ONENRM = NORM == '1' || lsame(NORM, 'O');
  if (!ONENRM && !lsame(NORM, 'I')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (ANORM < ZERO) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZGTCON', -INFO.value);
    return;
  }

  // Quick return if possible

  RCOND.value = ZERO;
  if (N == 0) {
    RCOND.value = ONE;
    return;
  } else if (ANORM == ZERO) {
    return;
  }

  // Check that D(1:N) is non-zero.

  for (I = 1; I <= N; I++) {
    // 10
    if (D[I] == Complex.zero) return;
  } // 10

  AINVNM.value = ZERO;
  if (ONENRM) {
    KASE1 = 1;
  } else {
    KASE1 = 2;
  }
  KASE.value = 0;
  while (true) {
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;
    if (KASE.value == KASE1) {
      // Multiply by inv(U)*inv(L).

      zgttrs('No transpose', N, 1, DL, D, DU, DU2, IPIV, WORK.asMatrix(N), N,
          INFO);
    } else {
      // Multiply by inv(L**H)*inv(U**H).

      zgttrs('Conjugate transpose', N, 1, DL, D, DU, DU2, IPIV,
          WORK.asMatrix(N), N, INFO);
    }
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
