import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/nio.dart';

void chkxer(
  final String SRNAMT,
  final int INFOT,
  final Nout NOUT,
  final Box<bool> LERR,
  final Box<bool> OK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  if (!LERR.value) {
    NOUT.println(
        ' *** Illegal value of parameter number ${INFOT.i2} not detected by ${SRNAMT.trim().a6} ***');
    OK.value = false;
  }
  LERR.value = false;
}
