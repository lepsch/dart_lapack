import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void ztpttr(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final A = A_.having(ld: LDA);
  bool LOWER;
  int I, J, K;

  // Test the input parameters.

  INFO.value = 0;
  LOWER = lsame(UPLO, 'L');
  if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZTPTTR', -INFO.value);
    return;
  }

  if (LOWER) {
    K = 0;
    for (J = 1; J <= N; J++) {
      for (I = J; I <= N; I++) {
        K++;
        A[I][J] = AP[K];
      }
    }
  } else {
    K = 0;
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J; I++) {
        K++;
        A[I][J] = AP[K];
      }
    }
  }
}
