import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dlatbs.dart';
import 'package:lapack/src/drscl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpbcon(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<double> AB_,
  final int LDAB,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool UPPER;
  String NORMIN;
  int IX;
  double SCALE, SMLNUM;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0), SCALEL = Box(0.0), SCALEU = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KD < 0) {
    INFO.value = -3;
  } else if (LDAB < KD + 1) {
    INFO.value = -5;
  } else if (ANORM < ZERO) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DPBCON', -INFO.value);
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

  SMLNUM = dlamch('Safe minimum');

  // Estimate the 1-norm of the inverse.

  KASE.value = 0;
  NORMIN = 'N';
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    if (UPPER) {
      // Multiply by inv(U**T).

      dlatbs('Upper', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK,
          SCALEL, WORK(2 * N + 1), INFO);
      NORMIN = 'Y';

      // Multiply by inv(U).

      dlatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK,
          SCALEU, WORK(2 * N + 1), INFO);
    } else {
      // Multiply by inv(L).

      dlatbs('Lower', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK,
          SCALEL, WORK(2 * N + 1), INFO);
      NORMIN = 'Y';

      // Multiply by inv(L**T).

      dlatbs('Lower', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK,
          SCALEU, WORK(2 * N + 1), INFO);
    }

    // Multiply by 1/SCALE if doing so will not cause overflow.

    SCALE = SCALEL.value * SCALEU.value;
    if (SCALE != ONE) {
      IX = idamax(N, WORK, 1);
      if (SCALE < WORK[IX].abs() * SMLNUM || SCALE == ZERO) return;
      drscl(N, SCALE, WORK, 1);
    }
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
