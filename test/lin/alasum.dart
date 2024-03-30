import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/nio.dart';

void alasum(
  final String TYPE,
  final Nout NOUT,
  final int NFAIL,
  final int NRUN,
  final int NERRS,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  if (NFAIL > 0) {
    NOUT.println(
        ' ${TYPE.a3}: ${NFAIL.i6} out of ${NRUN.i6} tests failed to pass the threshold');
  } else {
    NOUT.println(
        '\n All tests for ${TYPE.a3} routines passed the threshold ( ${NRUN.i6} tests run)');
  }
  if (NERRS > 0) {
    NOUT.println('      ${NERRS.i6} error messages recorded');
  }
}
