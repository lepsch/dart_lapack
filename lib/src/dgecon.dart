import 'dart:math';

import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dlatrs.dart';
import 'package:lapack/src/drscl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgecon(
  final String NORM,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final double ANORM,
  final Box<double> RCOND,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool ONENRM;
  String NORMIN;
  int IX, KASE1;
  double SCALE, SMLNUM, HUGEVAL;
  final ISAVE = Array<int>(3);
  final KASE = Box(0);
  final AINVNM = Box(0.0), SL = Box(0.0), SU = Box(0.0);

  HUGEVAL = dlamch('Overflow');

  // Test the input parameters.

  INFO.value = 0;
  ONENRM = NORM == '1' || lsame(NORM, 'O');
  if (!ONENRM && !lsame(NORM, 'I')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (ANORM < ZERO) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DGECON', -INFO.value);
    return;
  }

  // Quick return if possible

  RCOND.value = ZERO;
  if (N == 0) {
    RCOND.value = ONE;
    return;
  } else if (ANORM == ZERO) {
    return;
  } else if (disnan(ANORM)) {
    RCOND.value = ANORM;
    INFO.value = -5;
    return;
  } else if (ANORM > HUGEVAL) {
    INFO.value = -5;
    return;
  }

  SMLNUM = dlamch('Safe minimum');

  // Estimate the norm of inv(A).

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
      // Multiply by inv(L).

      dlatrs('Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL,
          WORK(2 * N + 1), INFO);

      // Multiply by inv(U).

      dlatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU,
          WORK(3 * N + 1), INFO);
    } else {
      // Multiply by inv(U**T).

      dlatrs('Upper', 'Transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU,
          WORK(3 * N + 1), INFO);

      // Multiply by inv(L**T).

      dlatrs('Lower', 'Transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL,
          WORK(2 * N + 1), INFO);
    }

    // Divide X by 1/(SL*SU) if doing so will not cause overflow.

    SCALE = SL.value * SU.value;
    NORMIN = 'Y';
    if (SCALE != ONE) {
      IX = idamax(N, WORK, 1);
      if (SCALE < WORK[IX].abs() * SMLNUM || SCALE == ZERO) return;
      drscl(N, SCALE, WORK, 1);
    }
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) {
    RCOND.value = (ONE / AINVNM.value) / ANORM;
  } else {
    INFO.value = 1;
    return;
  }

  // Check for NaNs and Infs

  if (disnan(RCOND.value) || RCOND.value > HUGEVAL) INFO.value = 1;
}
