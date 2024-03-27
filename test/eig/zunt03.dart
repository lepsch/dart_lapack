import 'dart:math';

import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

import 'zunt01.dart';

void zunt03(
  final String RC,
  final int MU,
  final int MV,
  final int N,
  final int K,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Box<double> RESULT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDV);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, IRC, J, LMX;
  double RES1, ULP;
  Complex S, SU, SV;
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
    xerbla('ZUNT03', -INFO.value);
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
      LMX = izamax(N, U(I, 1).asArray(), LDU);
      if (V[I][LMX] == Complex.zero) {
        SV = Complex.one;
      } else {
        SV = V[I][LMX].abs().toComplex() / V[I][LMX];
      }
      if (U[I][LMX] == Complex.zero) {
        SU = Complex.one;
      } else {
        SU = U[I][LMX].abs().toComplex() / U[I][LMX];
      }
      S = SV / SU;
      for (J = 1; J <= N; J++) {
        RES1 = max(RES1, (U[I][J] - S * V[I][J]).abs());
      }
    }
    RES1 /= N * ULP;

    // Compute orthogonality of rows of V.

    zunt01('Rows', MV, N, V, LDV, WORK, LWORK, RWORK, RES2);
  } else {
    // Compare columns

    RES1 = ZERO;
    for (I = 1; I <= K; I++) {
      LMX = izamax(N, U(1, I).asArray(), 1);
      if (V[LMX][I] == Complex.zero) {
        SV = Complex.one;
      } else {
        SV = V[LMX][I].abs().toComplex() / V[LMX][I];
      }
      if (U[LMX][I] == Complex.zero) {
        SU = Complex.one;
      } else {
        SU = U[LMX][I].abs().toComplex() / U[LMX][I];
      }
      S = SV / SU;
      for (J = 1; J <= N; J++) {
        RES1 = max(RES1, (U[J][I] - S * V[J][I]).abs());
      }
    }
    RES1 /= N * ULP;

    // Compute orthogonality of columns of V.

    zunt01('Columns', N, MV, V, LDV, WORK, LWORK, RWORK, RES2);
  }

  RESULT.value = min(max(RES1, RES2.value), ONE / ULP);
}
