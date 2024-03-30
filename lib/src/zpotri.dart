import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlauum.dart';
import 'package:lapack/src/ztrtri.dart';

void zpotri(
  final String UPLO,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);

  // Test the input parameters.

  INFO.value = 0;
  if (!lsame(UPLO, 'U') && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZPOTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Invert the triangular Cholesky factor U or L.

  ztrtri(UPLO, 'Non-unit', N, A, LDA, INFO);
  if (INFO.value > 0) return;

  // Form inv(U) * inv(U)**H or inv(L)**H * inv(L).

  zlauum(UPLO, N, A, LDA, INFO);
}
