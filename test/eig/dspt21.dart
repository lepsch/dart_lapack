import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dspmv.dart';
import 'package:lapack/src/blas/dspr.dart';
import 'package:lapack/src/blas/dspr2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dopmtr.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dspt21(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int KBAND,
  final Array<double> AP_,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> U_,
  final int LDU,
  final Array<double> VP_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final D = D_.having();
  final E = E_.having();
  final U = U_.having(ld: LDU);
  final VP = VP_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  const HALF = 1.0 / 2.0;
  bool LOWER;
  String CUPLO;
  int J, JP, JP1, JR, LAP;
  double ANORM, TEMP, ULP, UNFL, VSAVE, WNORM = 0;
  final IINFO = Box(0);

  // 1)      Constants

  RESULT[1] = ZERO;
  if (ITYPE == 1) RESULT[2] = ZERO;
  if (N <= 0) return;

  LAP = (N * (N + 1)) ~/ 2;

  if (lsame(UPLO, 'U')) {
    LOWER = false;
    CUPLO = 'U';
  } else {
    LOWER = true;
    CUPLO = 'L';
  }

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');

  // Some Error Checks

  if (ITYPE < 1 || ITYPE > 3) {
    RESULT[1] = TEN / ULP;
    return;
  }

  // Do Test 1

  // Norm of A:

  if (ITYPE == 3) {
    ANORM = ONE;
  } else {
    ANORM = max(dlansp('1', CUPLO, N, AP, WORK), UNFL);
  }

  // Compute error matrix:

  if (ITYPE == 1) {
    // ITYPE=1: error = A - U S U**T

    dlaset('Full', N, N, ZERO, ZERO, WORK.asMatrix(N), N);
    dcopy(LAP, AP, 1, WORK, 1);

    for (J = 1; J <= N; J++) {
      dspr(CUPLO, N, -D[J], U(1, J).asArray(), 1, WORK);
    }

    if (N > 1 && KBAND == 1) {
      for (J = 1; J <= N - 1; J++) {
        dspr2(CUPLO, N, -E[J], U(1, J).asArray(), 1, U(1, J + 1).asArray(), 1,
            WORK);
      }
    }
    WNORM = dlansp('1', CUPLO, N, WORK, WORK(pow(N, 2).toInt() + 1));
  } else if (ITYPE == 2) {
    // ITYPE=2: error = V S V**T - A

    dlaset('Full', N, N, ZERO, ZERO, WORK.asMatrix(N), N);

    if (LOWER) {
      WORK[LAP] = D[N];
      for (J = N - 1; J >= 1; J--) {
        JP = ((2 * N - J) * (J - 1)) ~/ 2;
        JP1 = JP + N - J;
        if (KBAND == 1) {
          WORK[JP + J + 1] = (ONE - TAU[J]) * E[J];
          for (JR = J + 2; JR <= N; JR++) {
            WORK[JP + JR] = -TAU[J] * E[J] * VP[JP + JR];
          }
        }

        if (TAU[J] != ZERO) {
          VSAVE = VP[JP + J + 1];
          VP[JP + J + 1] = ONE;
          dspmv('L', N - J, ONE, WORK(JP1 + J + 1), VP(JP + J + 1), 1, ZERO,
              WORK(LAP + 1), 1);
          TEMP =
              -HALF * TAU[J] * ddot(N - J, WORK(LAP + 1), 1, VP(JP + J + 1), 1);
          daxpy(N - J, TEMP, VP(JP + J + 1), 1, WORK(LAP + 1), 1);
          dspr2('L', N - J, -TAU[J], VP(JP + J + 1), 1, WORK(LAP + 1), 1,
              WORK(JP1 + J + 1));
          VP[JP + J + 1] = VSAVE;
        }
        WORK[JP + J] = D[J];
      }
    } else {
      WORK[1] = D[1];
      for (J = 1; J <= N - 1; J++) {
        JP = (J * (J - 1)) ~/ 2;
        JP1 = JP + J;
        if (KBAND == 1) {
          WORK[JP1 + J] = (ONE - TAU[J]) * E[J];
          for (JR = 1; JR <= J - 1; JR++) {
            WORK[JP1 + JR] = -TAU[J] * E[J] * VP[JP1 + JR];
          }
        }

        if (TAU[J] != ZERO) {
          VSAVE = VP[JP1 + J];
          VP[JP1 + J] = ONE;
          dspmv('U', J, ONE, WORK, VP(JP1 + 1), 1, ZERO, WORK(LAP + 1), 1);
          TEMP = -HALF * TAU[J] * ddot(J, WORK(LAP + 1), 1, VP(JP1 + 1), 1);
          daxpy(J, TEMP, VP(JP1 + 1), 1, WORK(LAP + 1), 1);
          dspr2('U', J, -TAU[J], VP(JP1 + 1), 1, WORK(LAP + 1), 1, WORK);
          VP[JP1 + J] = VSAVE;
        }
        WORK[JP1 + J + 1] = D[J + 1];
      }
    }

    for (J = 1; J <= LAP; J++) {
      WORK[J] -= AP[J];
    }
    WNORM = dlansp('1', CUPLO, N, WORK, WORK(LAP + 1));
  } else if (ITYPE == 3) {
    // ITYPE=3: error = U V**T - I

    if (N < 2) return;
    dlacpy(' ', N, N, U, LDU, WORK.asMatrix(N), N);
    dopmtr('R', CUPLO, 'T', N, N, VP, TAU, WORK.asMatrix(N), N,
        WORK(pow(N, 2).toInt() + 1), IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = TEN / ULP;
      return;
    }

    for (J = 1; J <= N; J++) {
      WORK[(N + 1) * (J - 1) + 1] -= ONE;
    }

    WNORM = dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1));
  }

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (N * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, N) / (N * ULP);
    }
  }

  // Do Test 2

  // Compute  U U**T - I

  if (ITYPE == 1) {
    dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK.asMatrix(N), N);

    for (J = 1; J <= N; J++) {
      WORK[(N + 1) * (J - 1) + 1] -= ONE;
    }

    RESULT[2] = min(
            dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1)),
            N) /
        (N * ULP);
  }
}
