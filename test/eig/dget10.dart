import 'dart:math';

import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget10(
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> WORK_,
  final Box<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
  const ONE = 1.0, ZERO = 0.0;
  int J;
  double ANORM, EPS, UNFL, WNORM;

  // Quick return if possible

  if (M <= 0 || N <= 0) {
    RESULT.value = ZERO;
    return;
  }

  UNFL = dlamch('Safe minimum');
  EPS = dlamch('Precision');

  WNORM = ZERO;
  for (J = 1; J <= N; J++) {
    dcopy(M, A(1, J).asArray(), 1, WORK, 1);
    daxpy(M, -ONE, B(1, J).asArray(), 1, WORK, 1);
    WNORM = max(WNORM, dasum(N, WORK, 1));
  }

  ANORM = max(dlange('1', M, N, A, LDA, WORK), UNFL);

  if (ANORM > WNORM) {
    RESULT.value = (WNORM / ANORM) / (M * EPS);
  } else {
    if (ANORM < ONE) {
      RESULT.value = (min(WNORM, M * ANORM) / ANORM) / (M * EPS);
    } else {
      RESULT.value = min(WNORM / ANORM, M.toDouble()) / (M * EPS);
    }
  }
}
