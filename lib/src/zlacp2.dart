import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zlacp2(
  final String UPLO,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  int I, J;

  if (lsame(UPLO, 'U')) {
    for (J = 1; J <= N; J++) {
      // 20
      for (I = 1; I <= min(J, M); I++) {
        // 10
        B[I][J] = A[I][J].toComplex();
      } // 10
    } // 20
  } else if (lsame(UPLO, 'L')) {
    for (J = 1; J <= N; J++) {
      // 40
      for (I = J; I <= M; I++) {
        // 30
        B[I][J] = A[I][J].toComplex();
      } // 30
    } // 40
  } else {
    for (J = 1; J <= N; J++) {
      // 60
      for (I = 1; I <= M; I++) {
        // 50
        B[I][J] = A[I][J].toComplex();
      } // 50
    } // 60
  }
}
