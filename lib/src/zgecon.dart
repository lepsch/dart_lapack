import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zdrscl.dart';
import 'package:lapack/src/zlacn2.dart';
import 'package:lapack/src/zlatrs.dart';

void zgecon(
  final String NORM,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final double ANORM,
  final Box<double> RCOND,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  const ONE = 1.0, ZERO = 0.0;
  bool ONENRM;
  String NORMIN;
  int IX, KASE1;
  double SCALE, SMLNUM, HUGEVAL;
  final ISAVE = Array<int>(3);
  final AINVNM = Box(0.0), SL = Box(0.0), SU = Box(0.0);
  final KASE = Box(0);

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

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
    xerbla('ZGECON', -INFO.value);
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
    zlacn2(N, WORK(N + 1), WORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;
    if (KASE.value == KASE1) {
      // Multiply by inv(L).

      zlatrs('Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL,
          RWORK, INFO);

      // Multiply by inv(U).

      zlatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU,
          RWORK(N + 1), INFO);
    } else {
      // Multiply by inv(U**H).

      zlatrs('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, A, LDA,
          WORK, SU, RWORK(N + 1), INFO);

      // Multiply by inv(L**H).

      zlatrs('Lower', 'Conjugate transpose', 'Unit', NORMIN, N, A, LDA, WORK,
          SL, RWORK, INFO);
    }

    // Divide X by 1/(SL*SU) if doing so will not cause overflow.

    SCALE = SL.value * SU.value;
    NORMIN = 'Y';
    if (SCALE != ONE) {
      IX = izamax(N, WORK, 1);
      if (SCALE < CABS1(WORK[IX]) * SMLNUM || SCALE == ZERO) return;
      zdrscl(N, SCALE, WORK, 1);
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
