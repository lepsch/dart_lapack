import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

double zla_gbrpvgrw(
  final int N,
  final int KL,
  final int KU,
  final int NCOLS,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> AFB_,
  final int LDAFB,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final AFB = AFB_.having(ld: LDAB);
  int I, J, KD;
  double AMAX, UMAX, RPVGRW;

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  RPVGRW = 1.0;

  KD = KU + 1;
  for (J = 1; J <= NCOLS; J++) {
    AMAX = 0.0;
    UMAX = 0.0;
    for (I = max(J - KU, 1); I <= min(J + KL, N); I++) {
      AMAX = max(CABS1(AB[KD + I - J][J]), AMAX);
    }
    for (I = max(J - KU, 1); I <= J; I++) {
      UMAX = max(CABS1(AFB[KD + I - J][J]), UMAX);
    }
    if (UMAX != 0.0) {
      RPVGRW = min(AMAX / UMAX, RPVGRW);
    }
  }
  return RPVGRW;
}
