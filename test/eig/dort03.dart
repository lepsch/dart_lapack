import 'dart:math';

import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

import 'dort01.dart';

void dort03(
  final String RC,
  final int MU,
  final int MV,
  final int N,
  final int K,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> V_,
  final int LDV,
  final Array<double> WORK_,
  final int LWORK,
  final Box<double> RESULT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();
  // List<String>       RC;
  // int                INFO.value, K, LDU, LDV, LWORK, MU, MV, N;
  // double             RESULT.value;
  // double             U( LDU, * ), V( LDV, * ), WORK( * );
  // // ..

  // double             ZERO, ONE;
  const ZERO = 0.0, ONE = 1.0;
  int I, IRC, J, LMX;
  double RES1, S, ULP;
  final RES2 = Box(0.0);

  // Check inputs

  INFO.value = 0;
  if (lsame(RC, 'R')) {
    IRC = 0;
  } else if (lsame(RC, 'C')) {
    IRC = 1;
  } else {
    IRC = -1;
  }
  if (IRC == -1) {
    INFO.value = -1;
  } else if (MU < 0) {
    INFO.value = -2;
  } else if (MV < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > max(MU, MV)) {
    INFO.value = -5;
  } else if ((IRC == 0 && LDU < max(1, MU)) || (IRC == 1 && LDU < max(1, N))) {
    INFO.value = -7;
  } else if ((IRC == 0 && LDV < max(1, MV)) || (IRC == 1 && LDV < max(1, N))) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DORT03', -INFO.value);
    return;
  }

  // Initialize result

  RESULT.value = ZERO;
  if (MU == 0 || MV == 0 || N == 0) return;

  // Machine constants

  ULP = dlamch('Precision');

  if (IRC == 0) {
    // Compare rows

    RES1 = ZERO;
    for (I = 1; I <= K; I++) {
      LMX = idamax(N, U(I, 1).asArray(), LDU);
      S = sign(ONE, U[I][LMX]) * sign(ONE, V[I][LMX]);
      for (J = 1; J <= N; J++) {
        RES1 = max(RES1, (U[I][J] - S * V[I][J]).abs());
      }
    }
    RES1 /= N * ULP;

    // Compute orthogonality of rows of V.

    dort01('Rows', MV, N, V, LDV, WORK, LWORK, RES2);
  } else {
    // Compare columns

    RES1 = ZERO;
    for (I = 1; I <= K; I++) {
      LMX = idamax(N, U(1, I).asArray(), 1);
      S = sign(ONE, U[LMX][I]) * sign(ONE, V[LMX][I]);
      for (J = 1; J <= N; J++) {
        RES1 = max(RES1, (U[J][I] - S * V[J][I]).abs());
      }
    }
    RES1 /= N * ULP;

    // Compute orthogonality of columns of V.

    dort01('Columns', N, MV, V, LDV, WORK, LWORK, RES2);
  }

  RESULT.value = min(max(RES1, RES2.value), ONE / ULP);
}
