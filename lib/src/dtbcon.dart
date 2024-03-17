import 'dart:math';

import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dlantb.dart';
import 'package:lapack/src/dlatbs.dart';
import 'package:lapack/src/drscl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtbcon(
  final String NORM,
  final String UPLO,
  final String DIAG,
  final int N,
  final int KD,
  final Matrix<double> AB_,
  final int LDAB,
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
  bool NOUNIT, ONENRM, UPPER;
  String NORMIN = '';
  int IX, KASE1 = 0;
  double ANORM, SMLNUM, XNORM;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0), SCALE = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  ONENRM = NORM == '1' || lsame(NORM, 'O');
  NOUNIT = lsame(DIAG, 'N');

  if (!ONENRM && !lsame(NORM, 'I')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (KD < 0) {
    INFO.value = -5;
  } else if (LDAB < KD + 1) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DTBCON', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) {
    RCOND.value = ONE;
    return;
  }

  RCOND.value = ZERO;
  SMLNUM = dlamch('Safe minimum') * (max(1, N)).toDouble();

  // Compute the norm of the triangular matrix A.

  ANORM = dlantb(NORM, UPLO, DIAG, N, KD, AB, LDAB, WORK);

  // Continue only if ANORM > 0.

  if (ANORM > ZERO) {
    // Estimate the norm of the inverse of A.

    AINVNM.value = ZERO;
    NORMIN = 'N';
    if (ONENRM) {
      KASE1 = 1;
    } else {
      KASE1 = 2;
    }
    KASE.value = 0;
    while (true) {
      dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
      if (KASE.value == 0) break;
      if (KASE.value == KASE1) {
        // Multiply by inv(A).

        dlatbs(UPLO, 'No transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE,
            WORK(2 * N + 1), INFO);
      } else {
        // Multiply by inv(A**T).

        dlatbs(UPLO, 'Transpose', DIAG, NORMIN, N, KD, AB, LDAB, WORK, SCALE,
            WORK(2 * N + 1), INFO);
      }
      NORMIN = 'Y';

      // Multiply by 1/SCALE if doing so will not cause overflow.

      if (SCALE.value != ONE) {
        IX = idamax(N, WORK, 1);
        XNORM = WORK[IX].abs();
        if (SCALE.value < XNORM * SMLNUM || SCALE.value == ZERO) return;
        drscl(N, SCALE.value, WORK, 1);
      }
    }

    // Compute the estimate of the reciprocal condition number.

    if (AINVNM.value != ZERO) RCOND.value = (ONE / ANORM) / AINVNM.value;
  }
}
