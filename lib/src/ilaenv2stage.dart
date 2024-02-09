import 'package:lapack/src/iparam2stage.dart';

int ilaenv2stage(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int N1,
  final int N2,
  final int N3,
  final int N4,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int IISPEC;
  switch (ISPEC) {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:

      // 2stage eigenvalues and SVD or related subroutines.

      IISPEC = 16 + ISPEC;
      return iparam2stage(IISPEC, NAME, OPTS, N1, N2, N3, N4);
    default:
      // Invalid value for ISPEC

      return -1;
  }
}
