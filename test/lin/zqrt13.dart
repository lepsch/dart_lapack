import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlarnv.dart';
import 'package:lapack/src/zlascl.dart';

void zqrt13(
  final int SCALE,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<double> NORMA,
  final Array<int> ISEED_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final ISEED = ISEED_.having(length: 4);

  const ONE = 1.0;
  final DUMMY = Array<double>(1);
  final INFO = Box(0);

  if (M <= 0 || N <= 0) return;

  // benign matrix

  for (var J = 1; J <= N; J++) {
    // 10
    zlarnv(2, ISEED, M, A(1, J).asArray());
    if (J <= M) {
      A[J][J] += Complex(
          sign(dzasum(M, A(1, J).asArray(), 1), A[J][J].real).toDouble());
    }
  } // 10

  // scaled versions

  if (SCALE != 1) {
    NORMA.value = zlange('Max', M, N, A, LDA, DUMMY);
    var SMLNUM = dlamch('Safe minimum');
    var BIGNUM = ONE / SMLNUM;
    SMLNUM = SMLNUM / dlamch('Epsilon');
    BIGNUM = ONE / SMLNUM;

    if (SCALE == 2) {
      // matrix scaled up

      zlascl('General', 0, 0, NORMA.value, BIGNUM, M, N, A, LDA, INFO);
    } else if (SCALE == 3) {
      // matrix scaled down

      zlascl('General', 0, 0, NORMA.value, SMLNUM, M, N, A, LDA, INFO);
    }
  }

  NORMA.value = zlange('One-norm', M, N, A, LDA, DUMMY);
}
