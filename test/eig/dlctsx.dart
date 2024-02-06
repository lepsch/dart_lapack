import 'common.dart';

bool dlctsx(AR, AI, BETA) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double AI, AR, BETA;
  // ..

// =====================================================================

  // .. Scalars in Common ..
  // bool               mn.FS;
  // int                mn.I, mn.M, mn.MPLUSN, mn.N;
  // ..
  // .. Common blocks ..
  // COMMON / mn / mn.M, mn.N, mn.MPLUSN, mn.I, mn.FS
  // ..
  // .. Save statement ..
  SAVE;

  if (mn.FS) {
    mn.I = mn.I + 1;
    if (mn.I <= mn.M) {
      DLCTSX = false;
    } else {
      DLCTSX = true;
    }
    if (mn.I == mn.MPLUSN) {
      mn.FS = false;
      mn.I = 0;
    }
  } else {
    mn.I = mn.I + 1;
    if (mn.I <= mn.N) {
      DLCTSX = true;
    } else {
      DLCTSX = false;
    }
    if (mn.I == mn.MPLUSN) {
      mn.FS = true;
      mn.I = 0;
    }
  }

  // IF( AR/BETA > 0.0 )THEN
  // DLCTSX = true;
  // ELSE
  // DLCTSX = false;
  // END IF
}
