import 'package:lapack/src/format_extensions.dart';

import 'common.dart';

void xerbla(final String SRNAME, final int INFO) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  infoc.LERR.value = true;
  if (INFO != infoc.INFOT) {
    if (infoc.INFOT != 0) {
      infoc.NOUT.println(
          ' *** XERBLA was called from ${srnamc.SRNAMT} with INFO = ${INFO.i6} instead of ${infoc.INFOT.i2} ***');
    } else {
      infoc.NOUT.println(
          ' *** On entry to $SRNAME parameter number ${INFO.i6} had an illegal value ***');
    }
    infoc.OK.value = false;
  }
  if (SRNAME != srnamc.SRNAMT) {
    infoc.NOUT.println(
        ' *** XERBLA was called with SRNAME = $SRNAME instead of ${srnamc.SRNAMT.a6} ***');
    infoc.OK.value = false;
  }
}
