import 'common.dart';

bool dlctsx(
  double AR,
  double AI,
  double BETA,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final bool result;

  if (mn.FS) {
    mn.I++;
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
    mn.I++;
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

  // IF( AR/BETA > 0.0 )THEN
  // DLCTSX = true;
  // ELSE
  // DLCTSX = false;
  // END IF

  return result;
}
