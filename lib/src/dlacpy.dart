import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dlacpy(
  final String UPLO,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= min(J, M); I++) {
        B[I][J] = A[I][J];
      }
    }
  } else if (lsame(UPLO, 'L')) {
    for (var J = 1; J <= N; J++) {
      for (var I = J; I <= M; I++) {
        B[I][J] = A[I][J];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= M; I++) {
        B[I][J] = A[I][J];
      }
    }
  }
}
