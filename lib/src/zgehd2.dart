import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfg.dart';

void zgehd2(
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  int I;
  final ALPHA = Box(Complex.zero);

  // Test the input parameters

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -2;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZGEHD2', -INFO.value);
    return;
  }

  for (I = ILO; I <= IHI - 1; I++) {
    // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)

    ALPHA.value = A[I + 1][I];
    zlarfg(IHI - I, ALPHA, A(min(I + 2, N), I).asArray(), 1, TAU(I));
    A[I + 1][I] = Complex.one;

    // Apply H(i) to A(1:ihi,i+1:ihi) from the right

    zlarf('Right', IHI, IHI - I, A(I + 1, I).asArray(), 1, TAU[I], A(1, I + 1),
        LDA, WORK);

    // Apply H(i)**H to A(i+1:ihi,i+1:n) from the left

    zlarf('Left', IHI - I, N - I, A(I + 1, I).asArray(), 1, TAU[I].conjugate(),
        A(I + 1, I + 1), LDA, WORK);

    A[I + 1][I] = ALPHA.value;
  }
}
