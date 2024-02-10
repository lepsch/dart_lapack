import 'package:lapack/src/intrinsics/sign.dart';

bool dlctes(final double ZR, final double ZI, final double D) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;

  if (D == ZERO) {
    return (ZR < ZERO);
  } else {
    return (sign(ONE, ZR) != sign(ONE, D));
  }
}
