import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dlaord(
  final String JOB,
  final int N,
  final Array<double> X_,
  final int INCX,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();

  final INC = INCX.abs();
  if (lsame(JOB, 'I')) {
    // Sort in increasing order

    sortLoop:
    for (var I = 2; I <= N; I++) {
      var IX = 1 + (I - 1) * INC;
      while (true) {
        if (IX == 1) continue sortLoop;
        final IXNEXT = IX - INC;
        if (X[IX] > X[IXNEXT]) {
          continue sortLoop;
        } else {
          final TEMP = X[IX];
          X[IX] = X[IXNEXT];
          X[IXNEXT] = TEMP;
        }
        IX = IXNEXT;
      }
    }
  } else if (lsame(JOB, 'D')) {
    // Sort in decreasing order

    sortLoop:
    for (var I = 2; I <= N; I++) {
      var IX = 1 + (I - 1) * INC;
      while (true) {
        if (IX == 1) continue sortLoop;
        final IXNEXT = IX - INC;
        if (X[IX] < X[IXNEXT]) {
          continue sortLoop;
        } else {
          final TEMP = X[IX];
          X[IX] = X[IXNEXT];
          X[IXNEXT] = TEMP;
        }
        IX = IXNEXT;
      }
    }
  }
}
