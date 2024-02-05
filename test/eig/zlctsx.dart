import 'package:lapack/src/complex.dart';

import 'common.dart';

bool zlctsx(ALPHA, BETA) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  Complex ALPHA, BETA;
  // ..

// =====================================================================

  // .. Parameters ..
  // double                         ZERO;
  // PARAMETER          ( ZERO = 0.0 )
  // Complex            CZERO
  // PARAMETER          ( CZERO = ( 0.0, 0.0 ) )
  // ..
  // .. Scalars in Common ..
  // bool mn.FS;
  // int mn.I, mn.M, mn.MPLUSN, mn.N;
  // ..
  // .. Common blocks ..
  // COMMON / MN / mn.M, mn.N, mn.MPLUSN, mn.I, mn.FS
  // ..
  // .. Save statement ..
  // SAVE;
  // ..
  // .. Executable Statements ..

  if (mn.FS) {
    mn.I = mn.I + 1;
    if (mn.I <= mn.M) {
      ZLCTSX = false;
    } else {
      ZLCTSX = true;
    }
    if (mn.I == mn.MPLUSN) {
      mn.FS = false;
      mn.I = 0;
    }
  } else {
    mn.I = mn.I + 1;
    if (mn.I <= mn.N) {
      ZLCTSX = true;
    } else {
      ZLCTSX = false;
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

  return;
}
