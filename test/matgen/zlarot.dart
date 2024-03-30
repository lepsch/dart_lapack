import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zlarot(
  final bool LROWS,
  final bool LLEFT,
  final bool LRIGHT,
  final int NL,
  final Complex C,
  final Complex S,
  final Array<Complex> A_,
  final int LDA,
  final Box<Complex> XLEFT,
  final Box<Complex> XRIGHT,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  int IINC, INEXT, IX, IY, IYT = 0, J, NT;
  Complex TEMPX;
  final XT = Array<Complex>(2), YT = Array<Complex>(2);

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
    xerbla('ZLAROT', 4);
    return;
  }
  if (LDA <= 0 || (!LROWS && LDA < NL - NT)) {
    xerbla('ZLAROT', 8);
    return;
  }

  // Rotate

  // ZROT( NL-NT, A(IX),IINC, A(IY),IINC, C, S ) with complex C, S

  for (J = 0; J <= NL - NT - 1; J++) {
    TEMPX = C * A[IX + J * IINC] + S * A[IY + J * IINC];
    A[IY + J * IINC] =
        -S.conjugate() * A[IX + J * IINC] + C.conjugate() * A[IY + J * IINC];
    A[IX + J * IINC] = TEMPX;
  }

  // ZROT( NT, XT,1, YT,1, C, S ) with complex C, S

  for (J = 1; J <= NT; J++) {
    TEMPX = C * XT[J] + S * YT[J];
    YT[J] = -S.conjugate() * XT[J] + C.conjugate() * YT[J];
    XT[J] = TEMPX;
  }

  // Stuff values back into XLEFT, XRIGHT, etc.

  if (LLEFT) {
    A[1] = XT[1];
    XLEFT.value = YT[1];
  }

  if (LRIGHT) {
    XRIGHT.value = XT[NT];
    A[IYT] = YT[NT];
  }
}
