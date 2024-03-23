import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyr.dart';
import 'package:lapack/src/blas/dsyr2.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dstt21(
  final int N,
  final int KBAND,
  final Array<double> AD_,
  final Array<double> AE_,
  final Array<double> SD_,
  final Array<double> SE_,
  final Matrix<double> U_,
  final int LDU,
  final Array<double> WORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AD = AD_.having();
  final AE = AE_.having();
  final SD = SD_.having();
  final SE = SE_.having();
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  int J;
  double ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM;

  // 1)      Constants

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');

  // Do Test 1

  // Copy A & Compute its 1-Norm:

  dlaset('Full', N, N, ZERO, ZERO, WORK.asMatrix(N), N);

  ANORM = ZERO;
  TEMP1 = ZERO;

  for (J = 1; J <= N - 1; J++) {
    WORK[(N + 1) * (J - 1) + 1] = AD[J];
    WORK[(N + 1) * (J - 1) + 2] = AE[J];
    TEMP2 = AE[J].abs();
    ANORM = max(ANORM, AD[J].abs() + TEMP1 + TEMP2);
    TEMP1 = TEMP2;
  }

  WORK[pow(N, 2).toInt()] = AD[N];
  ANORM = max(ANORM, max(AD[N].abs() + TEMP1, UNFL));

  // Norm of A - USU'

  for (J = 1; J <= N; J++) {
    dsyr('L', N, -SD[J], U(1, J).asArray(), 1, WORK.asMatrix(N), N);
  }

  if (N > 1 && KBAND == 1) {
    for (J = 1; J <= N - 1; J++) {
      dsyr2('L', N, -SE[J], U(1, J).asArray(), 1, U(1, J + 1).asArray(), 1,
          WORK.asMatrix(N), N);
    }
  }

  WNORM = dlansy('1', 'L', N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1));

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (N * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, N * ANORM) / ANORM) / (N * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, N.toDouble()) / (N * ULP);
    }
  }

  // Do Test 2

  // Compute  UU' - I

  dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK.asMatrix(N), N);

  for (J = 1; J <= N; J++) {
    WORK[(N + 1) * (J - 1) + 1] -= ONE;
  }

  RESULT[2] = min(N.toDouble(),
          dlange('1', N, N, WORK.asMatrix(N), N, WORK(pow(N, 2).toInt() + 1))) /
      (N * ULP);
}
