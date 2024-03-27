import 'dart:math';

import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zhpmv.dart';
import 'package:lapack/src/blas/zhpr.dart';
import 'package:lapack/src/blas/zhpr2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhp.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zupmtr.dart';

void zhpt21(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int KBAND,
  final Array<Complex> AP_,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> U_,
  final int LDU,
  final Array<Complex> VP_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final U = U_.having(ld: LDU);
  final AP = AP_.having();
  final VP = VP_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final D = D_.having();
  final E = E_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  const HALF = 1.0 / 2.0;
  bool LOWER;
  String CUPLO;
  int J, JP, JP1, JR, LAP;
  double ANORM, ULP, UNFL, WNORM = 0;
  Complex TEMP, VSAVE;
  final IINFO = Box(0);

  // Constants

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
    ANORM = max(zlanhp('1', CUPLO, N, AP, RWORK), UNFL);
  }

  // Compute error matrix:

  if (ITYPE == 1) {
    // ITYPE=1: error = A - U S U**H

    zlaset('Full', N, N, Complex.zero, Complex.zero, WORK.asMatrix(), N);
    zcopy(LAP, AP, 1, WORK, 1);

    for (J = 1; J <= N; J++) {
      zhpr(CUPLO, N, -D[J], U(1, J).asArray(), 1, WORK);
    }

    if (N > 1 && KBAND == 1) {
      for (J = 2; J <= N - 1; J++) {
        zhpr2(CUPLO, N, -E[J].toComplex(), U(1, J).asArray(), 1,
            U(1, J - 1).asArray(), 1, WORK);
      }
    }
    WNORM = zlanhp('1', CUPLO, N, WORK, RWORK);
  } else if (ITYPE == 2) {
    // ITYPE=2: error = V S V**H - A

    zlaset('Full', N, N, Complex.zero, Complex.zero, WORK.asMatrix(), N);

    if (LOWER) {
      WORK[LAP] = D[N].toComplex();
      for (J = N - 1; J >= 1; J--) {
        JP = ((2 * N - J) * (J - 1)) ~/ 2;
        JP1 = JP + N - J;
        if (KBAND == 1) {
          WORK[JP + J + 1] = (Complex.one - TAU[J]) * E[J].toComplex();
          for (JR = J + 2; JR <= N; JR++) {
            WORK[JP + JR] = -TAU[J] * E[J].toComplex() * VP[JP + JR];
          }
        }

        if (TAU[J] != Complex.zero) {
          VSAVE = VP[JP + J + 1];
          VP[JP + J + 1] = Complex.one;
          zhpmv('L', N - J, Complex.one, WORK(JP1 + J + 1), VP(JP + J + 1), 1,
              Complex.zero, WORK(LAP + 1), 1);
          TEMP = -HALF.toComplex() *
              TAU[J] *
              zdotc(N - J, WORK(LAP + 1), 1, VP(JP + J + 1), 1);
          zaxpy(N - J, TEMP, VP(JP + J + 1), 1, WORK(LAP + 1), 1);
          zhpr2('L', N - J, -TAU[J], VP(JP + J + 1), 1, WORK(LAP + 1), 1,
              WORK(JP1 + J + 1));

          VP[JP + J + 1] = VSAVE;
        }
        WORK[JP + J] = D[J].toComplex();
      }
    } else {
      WORK[1] = D[1].toComplex();
      for (J = 1; J <= N - 1; J++) {
        JP = (J * (J - 1)) ~/ 2;
        JP1 = JP + J;
        if (KBAND == 1) {
          WORK[JP1 + J] = (Complex.one - TAU[J]) * E[J].toComplex();
          for (JR = 1; JR <= J - 1; JR++) {
            WORK[JP1 + JR] = -TAU[J] * E[J].toComplex() * VP[JP1 + JR];
          }
        }

        if (TAU[J] != Complex.zero) {
          VSAVE = VP[JP1 + J];
          VP[JP1 + J] = Complex.one;
          zhpmv('U', J, Complex.one, WORK, VP(JP1 + 1), 1, Complex.zero,
              WORK(LAP + 1), 1);
          TEMP = -HALF.toComplex() *
              TAU[J] *
              zdotc(J, WORK(LAP + 1), 1, VP(JP1 + 1), 1);
          zaxpy(J, TEMP, VP(JP1 + 1), 1, WORK(LAP + 1), 1);
          zhpr2('U', J, -TAU[J], VP(JP1 + 1), 1, WORK(LAP + 1), 1, WORK);
          VP[JP1 + J] = VSAVE;
        }
        WORK[JP1 + J + 1] = D[J + 1].toComplex();
      }
    }

    for (J = 1; J <= LAP; J++) {
      WORK[J] -= AP[J];
    }
    WNORM = zlanhp('1', CUPLO, N, WORK, RWORK);
  } else if (ITYPE == 3) {
    // ITYPE=3: error = U V**H - I

    if (N < 2) return;
    zlacpy(' ', N, N, U, LDU, WORK.asMatrix(), N);
    zupmtr('R', CUPLO, 'C', N, N, VP, TAU, WORK.asMatrix(), N,
        WORK(pow(N, 2).toInt() + 1), IINFO);
    if (IINFO.value != 0) {
      RESULT[1] = TEN / ULP;
      return;
    }

    for (J = 1; J <= N; J++) {
      WORK[(N + 1) * (J - 1) + 1] -= Complex.one;
    }

    WNORM = zlange('1', N, N, WORK.asMatrix(), N, RWORK);
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

  // Compute  U U**H - I

  if (ITYPE == 1) {
    zgemm('N', 'C', N, N, N, Complex.one, U, LDU, U, LDU, Complex.zero,
        WORK.asMatrix(), N);

    for (J = 1; J <= N; J++) {
      WORK[(N + 1) * (J - 1) + 1] -= Complex.one;
    }

    RESULT[2] =
        min(zlange('1', N, N, WORK.asMatrix(), N, RWORK), N) / (N * ULP);
  }
}
