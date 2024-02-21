import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgtsv(
  final int N,
  final int NRHS,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.dim();
  final D = D_.dim();
  final DU = DU_.dim();
  final B = B_.dim(LDB);
  const ZERO = 0.0;
  int I = 0, J = 0;
  double FACT = 0, TEMP = 0;

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('DGTSV ', -INFO.value);
    return;
  }

  if (N == 0) return;

  if (NRHS == 1) {
    for (I = 1; I <= N - 2; I++) {
      if ((D[I]).abs() >= (DL[I]).abs()) {
        // No row interchange required

        if (D[I] != ZERO) {
          FACT = DL[I] / D[I];
          D[I + 1] = D[I + 1] - FACT * DU[I];
          B[I + 1][1] = B[I + 1][1] - FACT * B[I][1];
        } else {
          INFO.value = I;
          return;
        }
        DL[I] = ZERO;
      } else {
        // Interchange rows I and I+1

        FACT = D[I] / DL[I];
        D[I] = DL[I];
        TEMP = D[I + 1];
        D[I + 1] = DU[I] - FACT * TEMP;
        DL[I] = DU[I + 1];
        DU[I + 1] = -FACT * DL[I];
        DU[I] = TEMP;
        TEMP = B[I][1];
        B[I][1] = B[I + 1][1];
        B[I + 1][1] = TEMP - FACT * B[I + 1][1];
      }
    }
    if (N > 1) {
      I = N - 1;
      if ((D[I]).abs() >= (DL[I]).abs()) {
        if (D[I] != ZERO) {
          FACT = DL[I] / D[I];
          D[I + 1] = D[I + 1] - FACT * DU[I];
          B[I + 1][1] = B[I + 1][1] - FACT * B[I][1];
        } else {
          INFO.value = I;
          return;
        }
      } else {
        FACT = D[I] / DL[I];
        D[I] = DL[I];
        TEMP = D[I + 1];
        D[I + 1] = DU[I] - FACT * TEMP;
        DU[I] = TEMP;
        TEMP = B[I][1];
        B[I][1] = B[I + 1][1];
        B[I + 1][1] = TEMP - FACT * B[I + 1][1];
      }
    }
    if (D[N] == ZERO) {
      INFO.value = N;
      return;
    }
  } else {
    for (I = 1; I <= N - 2; I++) {
      if ((D[I]).abs() >= (DL[I]).abs()) {
        // No row interchange required

        if (D[I] != ZERO) {
          FACT = DL[I] / D[I];
          D[I + 1] = D[I + 1] - FACT * DU[I];
          for (J = 1; J <= NRHS; J++) {
            B[I + 1][J] = B[I + 1][J] - FACT * B[I][J];
          }
        } else {
          INFO.value = I;
          return;
        }
        DL[I] = ZERO;
      } else {
        // Interchange rows I and I+1

        FACT = D[I] / DL[I];
        D[I] = DL[I];
        TEMP = D[I + 1];
        D[I + 1] = DU[I] - FACT * TEMP;
        DL[I] = DU[I + 1];
        DU[I + 1] = -FACT * DL[I];
        DU[I] = TEMP;
        for (J = 1; J <= NRHS; J++) {
          TEMP = B[I][J];
          B[I][J] = B[I + 1][J];
          B[I + 1][J] = TEMP - FACT * B[I + 1][J];
        }
      }
    }
    if (N > 1) {
      I = N - 1;
      if ((D[I]).abs() >= (DL[I]).abs()) {
        if (D[I] != ZERO) {
          FACT = DL[I] / D[I];
          D[I + 1] = D[I + 1] - FACT * DU[I];
          for (J = 1; J <= NRHS; J++) {
            B[I + 1][J] = B[I + 1][J] - FACT * B[I][J];
          }
        } else {
          INFO.value = I;
          return;
        }
      } else {
        FACT = D[I] / DL[I];
        D[I] = DL[I];
        TEMP = D[I + 1];
        D[I + 1] = DU[I] - FACT * TEMP;
        DU[I] = TEMP;
        for (J = 1; J <= NRHS; J++) {
          TEMP = B[I][J];
          B[I][J] = B[I + 1][J];
          B[I + 1][J] = TEMP - FACT * B[I + 1][J];
        }
      }
    }
    if (D[N] == ZERO) {
      INFO.value = N;
      return;
    }
  }

  // Back solve with the matrix U from the factorization.

  if (NRHS <= 2) {
    J = 1;
    while (true) {
      B[N][J] = B[N][J] / D[N];
      if (N > 1) B[N - 1][J] = (B[N - 1][J] - DU[N - 1] * B[N][J]) / D[N - 1];
      for (I = N - 2; I >= 1; I--) {
        B[I][J] = (B[I][J] - DU[I] * B[I + 1][J] - DL[I] * B[I + 2][J]) / D[I];
      }
      if (J < NRHS) {
        J = J + 1;
        continue;
      }
      break;
    }
  } else {
    for (J = 1; J <= NRHS; J++) {
      B[N][J] = B[N][J] / D[N];
      if (N > 1) B[N - 1][J] = (B[N - 1][J] - DU[N - 1] * B[N][J]) / D[N - 1];
      for (I = N - 2; I >= 1; I--) {
        B[I][J] = (B[I][J] - DU[I] * B[I + 1][J] - DL[I] * B[I + 2][J]) / D[I];
      }
    }
  }
}
