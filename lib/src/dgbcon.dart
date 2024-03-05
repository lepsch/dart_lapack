import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dlatbs.dart';
import 'package:lapack/src/drscl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgbcon(
  final String NORM,
  final int N,
  final int KL,
  final int KU,
  final Matrix<double> AB_,
  final int LDAB,
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
  final AB = AB_.having(ld: LDAB);
  final IPIV = IPIV_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool LNOTI, ONENRM;
  String NORMIN;
  int IX, J, JP, KASE1, KD, LM;
  double SMLNUM, T;
  final ISAVE = Array<int>(3);
  final SCALE = Box(0.0), AINVNM = Box(0.0);
  final KASE = Box(0);

  // Test the input parameters.

  INFO.value = 0;
  ONENRM = NORM == '1' || lsame(NORM, 'O');
  if (!ONENRM && !lsame(NORM, 'I')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0) {
    INFO.value = -3;
  } else if (KU < 0) {
    INFO.value = -4;
  } else if (LDAB < 2 * KL + KU + 1) {
    INFO.value = -6;
  } else if (ANORM < ZERO) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('DGBCON', -INFO.value);
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

  // Estimate the norm of inv(A).

  AINVNM.value = ZERO;
  NORMIN = 'N';
  if (ONENRM) {
    KASE1 = 1;
  } else {
    KASE1 = 2;
  }
  KD = KL + KU + 1;
  LNOTI = KL > 0;
  KASE.value = 0;
  while (true) {
    dlacn2(N, WORK(N + 1), WORK, IWORK, AINVNM, KASE, ISAVE);
    if (KASE.value == 0) break;

    if (KASE.value == KASE1) {
      // Multiply by inv(L).

      if (LNOTI) {
        for (J = 1; J <= N - 1; J++) {
          LM = min(KL, N - J);
          JP = IPIV[J];
          T = WORK[JP];
          if (JP != J) {
            WORK[JP] = WORK[J];
            WORK[J] = T;
          }
          daxpy(LM, -T, AB(KD + 1, J).asArray(), 1, WORK(J + 1), 1);
        }
      }

      // Multiply by inv(U).

      dlatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KL + KU, AB, LDAB,
          WORK, SCALE, WORK(2 * N + 1), INFO);
    } else {
      // Multiply by inv(U**T).

      dlatbs('Upper', 'Transpose', 'Non-unit', NORMIN, N, KL + KU, AB, LDAB,
          WORK, SCALE, WORK(2 * N + 1), INFO);

      // Multiply by inv(L**T).

      if (LNOTI) {
        for (J = N - 1; J >= 1; J--) {
          LM = min(KL, N - J);
          WORK[J] =
              WORK[J] - ddot(LM, AB(KD + 1, J).asArray(), 1, WORK(J + 1), 1);
          JP = IPIV[J];
          if (JP != J) {
            T = WORK[JP];
            WORK[JP] = WORK[J];
            WORK[J] = T;
          }
        }
      }
    }

    // Divide X by 1/SCALE if doing so will not cause overflow.

    NORMIN = 'Y';
    if (SCALE.value != ONE) {
      IX = idamax(N, WORK, 1);
      if (SCALE.value < WORK[IX].abs() * SMLNUM || SCALE.value == ZERO) return;
      drscl(N, SCALE.value, WORK, 1);
    }
  }

  // Compute the estimate of the reciprocal condition number.

  if (AINVNM.value != ZERO) RCOND.value = (ONE / AINVNM.value) / ANORM;
}
