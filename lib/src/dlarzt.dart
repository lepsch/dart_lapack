import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlarzt(
  final String DIRECT,
  final String STOREV,
  final int N,
  final int K,
  final Matrix<double> V_,
  final int LDV,
  final Array<double> TAU_,
  final Matrix<double> T_,
  final int LDT,
) {
  final V = V_.dim(LDV);
  final TAU = TAU_.dim();
  final T = T_.dim(LDT);

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;

  // Check for currently supported options

  var INFO = 0;
  if (!lsame(DIRECT, 'B')) {
    INFO = -1;
  } else if (!lsame(STOREV, 'R')) {
    INFO = -2;
  }
  if (INFO != 0) {
    xerbla('DLARZT', -INFO);
    return;
  }

  for (var I = K; I >= 1; I--) {
    if (TAU[I] == ZERO) {
      // H(i)  =  I

      for (var J = I; J <= K; J++) {
        T[J][I] = ZERO;
      }
    } else {
      // general case

      if (I < K) {
        // T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**T

        dgemv('No transpose', K - I, N, -TAU[I], V(I + 1, 1), LDV,
            V(I, 1).asArray(), LDV, ZERO, T(I + 1, I).asArray(), 1);

        // T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)

        dtrmv('Lower', 'No transpose', 'Non-unit', K - I, T(I + 1, I + 1), LDT,
            T(I + 1, I).asArray(), 1);
      }
      T[I][I] = TAU[I];
    }
  }
}
