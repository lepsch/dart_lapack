import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlarot(
  final bool LROWS,
  final bool LLEFT,
  final bool LRIGHT,
  final int NL,
  final double C,
  final double S,
  final Array<double> A_,
  final int LDA,
  final Box<double> XLEFT,
  final Box<double> XRIGHT,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  int IINC, INEXT, IX, IY, IYT = 0, NT;
  final XT = Array<double>(2), YT = Array<double>(2);

  // Set up indices, arrays for ends

  if (LROWS) {
    IINC = LDA;
    INEXT = 1;
  } else {
    IINC = 1;
    INEXT = LDA;
  }

  if (LLEFT) {
    NT = 1;
    IX = 1 + IINC;
    IY = 2 + LDA;
    XT[1] = A[1];
    YT[1] = XLEFT.value;
  } else {
    NT = 0;
    IX = 1;
    IY = 1 + INEXT;
  }

  if (LRIGHT) {
    IYT = 1 + INEXT + (NL - 1) * IINC;
    NT++;
    XT[NT] = XRIGHT.value;
    YT[NT] = A[IYT];
  }

  // Check for errors

  if (NL < NT) {
    xerbla('DLAROT', 4);
    return;
  }
  if (LDA <= 0 || (!LROWS && LDA < NL - NT)) {
    xerbla('DLAROT', 8);
    return;
  }

  // Rotate

  drot(NL - NT, A(IX), IINC, A(IY), IINC, C, S);
  drot(NT, XT, 1, YT, 1, C, S);

  // Stuff values back into XLEFT.value, XRIGHT.value, etc.

  if (LLEFT) {
    A[1] = XT[1];
    XLEFT.value = YT[1];
  }

  if (LRIGHT) {
    XRIGHT.value = XT[NT];
    A[IYT] = YT[NT];
  }
}
