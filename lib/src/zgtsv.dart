import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgtsv(
  final int N,
  final int NRHS,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final B = B_.having(ld: LDB);
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  int J, K;
  Complex MULT, TEMP;

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZGTSV', -INFO.value);
    return;
  }

  if (N == 0) return;

  for (K = 1; K <= N - 1; K++) {
    if (DL[K] == Complex.zero) {
      // Subdiagonal is zero, no elimination is required.

      if (D[K] == Complex.zero) {
        // Diagonal is zero: set INFO = K and return; a unique
        // solution can not be found.

        INFO.value = K;
        return;
      }
    } else if (D[K].cabs1() >= DL[K].cabs1()) {
      // No row interchange required

      MULT = DL[K] / D[K];
      D[K + 1] -= MULT * DU[K];
      for (J = 1; J <= NRHS; J++) {
        B[K + 1][J] -= MULT * B[K][J];
      }
      if (K < (N - 1)) DL[K] = Complex.zero;
    } else {
      // Interchange rows K and K+1

      MULT = D[K] / DL[K];
      D[K] = DL[K];
      TEMP = D[K + 1];
      D[K + 1] = DU[K] - MULT * TEMP;
      if (K < (N - 1)) {
        DL[K] = DU[K + 1];
        DU[K + 1] = -MULT * DL[K];
      }
      DU[K] = TEMP;
      for (J = 1; J <= NRHS; J++) {
        TEMP = B[K][J];
        B[K][J] = B[K + 1][J];
        B[K + 1][J] = TEMP - MULT * B[K + 1][J];
      }
    }
  }
  if (D[N] == Complex.zero) {
    INFO.value = N;
    return;
  }

  // Back solve with the matrix U from the factorization.

  for (J = 1; J <= NRHS; J++) {
    B[N][J] /= D[N];
    if (N > 1) B[N - 1][J] = (B[N - 1][J] - DU[N - 1] * B[N][J]) / D[N - 1];
    for (K = N - 2; K >= 1; K--) {
      B[K][J] = (B[K][J] - DU[K] * B[K + 1][J] - DL[K] * B[K + 2][J]) / D[K];
    }
  }
}
