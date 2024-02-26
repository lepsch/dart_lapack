import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';

void zlarzt(
  final String DIRECT,
  final String STOREV,
  final int N,
  final int K,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> TAU_,
  final Matrix<Complex> T_,
  final int LDT,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.dim(LDV);
  final T = T_.dim(LDT);
  final TAU = TAU_.dim();
  int I, J;
  final INFO = Box(0);

  // Check for currently supported options

  INFO.value = 0;
  if (!lsame(DIRECT, 'B')) {
    INFO.value = -1;
  } else if (!lsame(STOREV, 'R')) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('ZLARZT', -INFO.value);
    return;
  }

  for (I = K; I >= 1; I--) {
    // 20
    if (TAU[I] == Complex.zero) {
      // H(i)  =  I

      for (J = I; J <= K; J++) {
        // 10
        T[J][I] = Complex.zero;
      } // 10
    } else {
      // general case

      if (I < K) {
        // T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**H

        zlacgv(N, V(I, 1).asArray(), LDV);
        zgemv('No transpose', K - I, N, -TAU[I], V(I + 1, 1), LDV,
            V(I, 1).asArray(), LDV, Complex.zero, T(I + 1, I).asArray(), 1);
        zlacgv(N, V(I, 1).asArray(), LDV);

        // T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)

        ztrmv('Lower', 'No transpose', 'Non-unit', K - I, T(I + 1, I + 1), LDT,
            T(I + 1, I).asArray(), 1);
      }
      T[I][I] = TAU[I];
    }
  } // 20
}
