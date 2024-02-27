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
  final B = B_.dim(LDB);
  final DL = DL_.dim();
  final D = D_.dim();
  final DU = DU_.dim();
  int J, K;
  Complex MULT, TEMP;

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (NRHS < 0) {
    INFO.value = -2;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  }
  if (INFO.value != 0) {
    xerbla('ZGTSV ', -INFO.value);
    return;
  }

  if (N == 0) return;

  for (K = 1; K <= N - 1; K++) {
    // 30
    if (DL[K] == Complex.zero) {
      // Subdiagonal is zero, no elimination is required.

      if (D[K] == Complex.zero) {
        // Diagonal is zero: set INFO.value = K and return; a unique
        // solution can not be found.

        INFO.value = K;
        return;
      }
    } else if (CABS1(D[K]) >= CABS1(DL[K])) {
      // No row interchange required

      MULT = DL[K] / D[K];
      D[K + 1] = D[K + 1] - MULT * DU[K];
      for (J = 1; J <= NRHS; J++) {
        // 10
        B[K + 1][J] = B[K + 1][J] - MULT * B[K][J];
      } // 10
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
        // 20
        TEMP = B[K][J];
        B[K][J] = B[K + 1][J];
        B[K + 1][J] = TEMP - MULT * B[K + 1][J];
      } // 20
    }
  } // 30
  if (D[N] == Complex.zero) {
    INFO.value = N;
    return;
  }

  // Back solve with the matrix U from the factorization.

  for (J = 1; J <= NRHS; J++) {
    // 50
    B[N][J] = B[N][J] / D[N];
    if (N > 1) B[N - 1][J] = (B[N - 1][J] - DU[N - 1] * B[N][J]) / D[N - 1];
    for (K = N - 2; K >= 1; K--) {
      // 40
      B[K][J] = (B[K][J] - DU[K] * B[K + 1][J] - DL[K] * B[K + 2][J]) / D[K];
    } // 40
  } // 50
}
