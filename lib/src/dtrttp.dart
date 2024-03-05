import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtrttp(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> AP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AP = AP_.having();
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
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('DTRTTP', -INFO.value);
    return;
  }

  if (LOWER) {
    K = 0;
    for (J = 1; J <= N; J++) {
      for (I = J; I <= N; I++) {
        K = K + 1;
        AP[K] = A[I][J];
      }
    }
  } else {
    K = 0;
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= J; I++) {
        K = K + 1;
        AP[K] = A[I][J];
      }
    }
  }
}
