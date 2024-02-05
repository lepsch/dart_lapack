import 'package:lapack/src/dlapy2.dart';

import 'common.dart';

bool dslect(final double ZR, final double ZI) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I;
  double RMIN, X;
  const ZERO = 0.0;

  if (sslct.SELOPT == 0) {
    return (ZR < ZERO);
  }
  RMIN = dlapy2(ZR - sslct.SELWR[1], ZI - sslct.SELWI[1]);
  var result = sslct.SELVAL[1];
  for (I = 2; I <= sslct.SELDIM; I++) {
    X = dlapy2(ZR - sslct.SELWR[I], ZI - sslct.SELWI[I]);
    if (X <= RMIN) {
      RMIN = X;
      result = sslct.SELVAL[I];
    }
  }

  return result;
}
