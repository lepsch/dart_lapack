import 'package:lapack/src/complex.dart';

import 'common.dart';

bool zlctsx(final Complex ALPHA, final Complex BETA) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  bool result;
  if (mn.FS) {
    mn.I = mn.I + 1;
    if (mn.I <= mn.M) {
      result = false;
    } else {
      result = true;
    }
    if (mn.I == mn.MPLUSN) {
      mn.FS = false;
      mn.I = 0;
    }
  } else {
    mn.I = mn.I + 1;
    if (mn.I <= mn.N) {
      result = true;
    } else {
      result = false;
    }
    if (mn.I == mn.MPLUSN) {
      mn.FS = true;
      mn.I = 0;
    }
  }

  // IF( BETA == CZERO ) THEN
  // ZLCTSX = ( ALPHA.toDouble() > ZERO )
  // ELSE
  // ZLCTSX = ( (ALPHA/BETA).toDouble() > ZERO )
  // END IF

  return result;
}
