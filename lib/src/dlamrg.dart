import 'package:lapack/src/matrix.dart';

void dlamrg(
  final int N1,
  final int N2,
  final Array<double> A_,
  final int DTRD1,
  final int DTRD2,
  final Array<int> INDEX,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  int I, IND1, IND2, N1SV, N2SV;

  N1SV = N1;
  N2SV = N2;
  if (DTRD1 > 0) {
    IND1 = 1;
  } else {
    IND1 = N1;
  }
  if (DTRD2 > 0) {
    IND2 = 1 + N1;
  } else {
    IND2 = N1 + N2;
  }
  I = 1;
  while (N1SV > 0 && N2SV > 0) {
    if (A[IND1] <= A[IND2]) {
      INDEX[I] = IND1;
      I = I + 1;
      IND1 = IND1 + DTRD1;
      N1SV = N1SV - 1;
    } else {
      INDEX[I] = IND2;
      I = I + 1;
      IND2 = IND2 + DTRD2;
      N2SV = N2SV - 1;
    }
  }
  // end while
  if (N1SV == 0) {
    for (N1SV = 1; N1SV <= N2SV; N1SV++) {
      INDEX[I] = IND2;
      I = I + 1;
      IND2 = IND2 + DTRD2;
    }
  } else {
    // N2SV == 0
    for (N2SV = 1; N2SV <= N1SV; N2SV++) {
      INDEX[I] = IND1;
      I = I + 1;
      IND1 = IND1 + DTRD1;
    }
  }
}
