import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarfg.dart';
import 'package:lapack/src/zlarz.dart';

void zlatrz(
  final int M,
  final int N,
  final int L,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  int I;
  final ALPHA = Box(Complex.zero);

  // Quick return if possible

  if (M == 0) {
    return;
  } else if (M == N) {
    for (I = 1; I <= N; I++) {
      TAU[I] = Complex.zero;
    }
    return;
  }

  for (I = M; I >= 1; I--) {
    // Generate elementary reflector H(i) to annihilate
    // [ A(i,i) A(i,n-l+1:n) ]

    zlacgv(L, A(I, N - L + 1).asArray(), LDA);
    ALPHA.value = A[I][I].conjugate();
    zlarfg(L + 1, ALPHA, A(I, N - L + 1).asArray(), LDA, TAU(I));
    TAU[I] = TAU[I].conjugate();

    // Apply H(i) to A(1:i-1,i:n) from the right

    zlarz('Right', I - 1, N - I + 1, L, A(I, N - L + 1).asArray(), LDA,
        TAU[I].conjugate(), A(1, I), LDA, WORK);
    A[I][I] = ALPHA.value.conjugate();
  }
}
