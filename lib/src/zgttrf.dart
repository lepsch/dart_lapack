import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgttrf(
  final int N,
  final Array<Complex> DL_,
  final Array<Complex> D_,
  final Array<Complex> DU_,
  final Array<Complex> DU2_,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  final IPIV = IPIV_.having();
  const ZERO = 0.0;
  int I;
  Complex FACT, TEMP;

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
    xerbla('ZGTTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Initialize IPIV(i) = i and DU2(i) = 0

  for (I = 1; I <= N; I++) {
    // 10
    IPIV[I] = I;
  } // 10
  for (I = 1; I <= N - 2; I++) {
    // 20
    DU2[I] = Complex.zero;
  } // 20

  for (I = 1; I <= N - 2; I++) {
    // 30
    if (CABS1(D[I]) >= CABS1(DL[I])) {
      // No row interchange required, eliminate DL(I)

      if (CABS1(D[I]) != ZERO) {
        FACT = DL[I] / D[I];
        DL[I] = FACT;
        D[I + 1] = D[I + 1] - FACT * DU[I];
      }
    } else {
      // Interchange rows I and I+1, eliminate DL(I)

      FACT = D[I] / DL[I];
      D[I] = DL[I];
      DL[I] = FACT;
      TEMP = DU[I];
      DU[I] = D[I + 1];
      D[I + 1] = TEMP - FACT * D[I + 1];
      DU2[I] = DU[I + 1];
      DU[I + 1] = -FACT * DU[I + 1];
      IPIV[I] = I + 1;
    }
  } // 30
  if (N > 1) {
    I = N - 1;
    if (CABS1(D[I]) >= CABS1(DL[I])) {
      if (CABS1(D[I]) != ZERO) {
        FACT = DL[I] / D[I];
        DL[I] = FACT;
        D[I + 1] = D[I + 1] - FACT * DU[I];
      }
    } else {
      FACT = D[I] / DL[I];
      D[I] = DL[I];
      DL[I] = FACT;
      TEMP = DU[I];
      DU[I] = D[I + 1];
      D[I + 1] = TEMP - FACT * D[I + 1];
      IPIV[I] = I + 1;
    }
  }

  // Check for a zero on the diagonal of U.

  for (I = 1; I <= N; I++) {
    if (CABS1(D[I]) == ZERO) {
      INFO.value = I;
      return;
    }
  }
}
