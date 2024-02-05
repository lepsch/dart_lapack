import 'common.dart';

void xlaenv(final int ISPEC, final int NVALUE) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  if (ISPEC >= 1 && ISPEC <= 16) {
    claenv.IPARMS[ISPEC] = NVALUE;
  }
}
