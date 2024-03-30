import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zher.dart';
import 'package:lapack/src/blas/zher2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlarfy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zunm2l.dart';
import 'package:lapack/src/zunm2r.dart';

void zhet21(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int KBAND,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final D = D_.having();
  final E = E_.having();
  final RESULT = RESULT_.having(length: 2);
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  bool LOWER;
  String CUPLO;
  int J, JCOL, JR, JROW;
  double ANORM, ULP, UNFL, WNORM = 0;
  Complex VSAVE;
  final IINFO = Box(0);

  RESULT[1] = ZERO;
  if (ITYPE == 1) RESULT[2] = ZERO;
  if (N <= 0) return;

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
    ANORM = max(zlanhe('1', CUPLO, N, A, LDA, RWORK), UNFL);
  }

  // Compute error matrix:

  if (ITYPE == 1) {
    // ITYPE=1: error = A - U S U**H

    zlaset('Full', N, N, Complex.zero, Complex.zero, WORK.asMatrix(), N);
    zlacpy(CUPLO, N, N, A, LDA, WORK.asMatrix(), N);

    for (J = 1; J <= N; J++) {
      zher(CUPLO, N, -D[J], U(1, J).asArray(), 1, WORK.asMatrix(), N);
    }

    if (N > 1 && KBAND == 1) {
      for (J = 1; J <= N - 1; J++) {
        zher2(CUPLO, N, -E[J].toComplex(), U(1, J).asArray(), 1,
            U(1, J + 1).asArray(), 1, WORK.asMatrix(), N);
      }
    }
    WNORM = zlanhe('1', CUPLO, N, WORK.asMatrix(), N, RWORK);
  } else if (ITYPE == 2) {
    // ITYPE=2: error = V S V**H - A

    zlaset('Full', N, N, Complex.zero, Complex.zero, WORK.asMatrix(), N);

    if (LOWER) {
      WORK[pow(N, 2).toInt()] = D[N].toComplex();
      for (J = N - 1; J >= 1; J--) {
        if (KBAND == 1) {
          WORK[(N + 1) * (J - 1) + 2] =
              (Complex.one - TAU[J]) * E[J].toComplex();
          for (JR = J + 2; JR <= N; JR++) {
            WORK[(J - 1) * N + JR] = -TAU[J] * E[J].toComplex() * V[JR][J];
          }
        }

        VSAVE = V[J + 1][J];
        V[J + 1][J] = Complex.one;
        zlarfy('L', N - J, V(J + 1, J).asArray(), 1, TAU[J],
            WORK((N + 1) * J + 1).asMatrix(), N, WORK(pow(N, 2).toInt() + 1));
        V[J + 1][J] = VSAVE;
        WORK[(N + 1) * (J - 1) + 1] = D[J].toComplex();
      }
    } else {
      WORK[1] = D[1].toComplex();
      for (J = 1; J <= N - 1; J++) {
        if (KBAND == 1) {
          WORK[(N + 1) * J] = (Complex.one - TAU[J]) * E[J].toComplex();
          for (JR = 1; JR <= J - 1; JR++) {
            WORK[J * N + JR] = -TAU[J] * E[J].toComplex() * V[JR][J + 1];
          }
        }

        VSAVE = V[J][J + 1];
        V[J][J + 1] = Complex.one;
        zlarfy('U', J, V(1, J + 1).asArray(), 1, TAU[J], WORK.asMatrix(), N,
            WORK(pow(N, 2).toInt() + 1));
        V[J][J + 1] = VSAVE;
        WORK[(N + 1) * J + 1] = D[J + 1].toComplex();
      }
    }

    for (JCOL = 1; JCOL <= N; JCOL++) {
      if (LOWER) {
        for (JROW = JCOL; JROW <= N; JROW++) {
          WORK[JROW + N * (JCOL - 1)] =
              WORK[JROW + N * (JCOL - 1)] - A[JROW][JCOL];
        }
      } else {
        for (JROW = 1; JROW <= JCOL; JROW++) {
          WORK[JROW + N * (JCOL - 1)] =
              WORK[JROW + N * (JCOL - 1)] - A[JROW][JCOL];
        }
      }
    }
    WNORM = zlanhe('1', CUPLO, N, WORK.asMatrix(), N, RWORK);
  } else if (ITYPE == 3) {
    // ITYPE=3: error = U V**H - I

    if (N < 2) return;
    zlacpy(' ', N, N, U, LDU, WORK.asMatrix(), N);
    if (LOWER) {
      zunm2r('R', 'C', N, N - 1, N - 1, V(2, 1), LDV, TAU,
          WORK(N + 1).asMatrix(), N, WORK(pow(N, 2).toInt() + 1), IINFO);
    } else {
      zunm2l('R', 'C', N, N - 1, N - 1, V(1, 2), LDV, TAU, WORK.asMatrix(), N,
          WORK(pow(N, 2).toInt() + 1), IINFO);
    }
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
