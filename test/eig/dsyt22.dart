import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsymm.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'dort01.dart';

void dsyt22(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int M,
  final int KBAND,
  final Matrix<double> A,
  final int LDA,
  final Array<double> D,
  final Array<double> E,
  final Matrix<double> U,
  final int LDU,
  final Matrix<double> V,
  final int LDV,
  final Array<double> TAU,
  final Array<double> WORK,
  final Array<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  int J, JJ, JJ1, JJ2, NN, NNP1;
  double ANORM, ULP, UNFL, WNORM;

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0 || M <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');

  // Do Test 1

  // Norm of A:

  ANORM = max(dlansy('1', UPLO, N, A, LDA, WORK), UNFL);

  // Compute error matrix:

  // ITYPE=1: error = U**T A U - S

  dsymm('L', UPLO, N, M, ONE, A, LDA, U, LDU, ZERO, WORK.asMatrix(N), N);
  NN = N * N;
  NNP1 = NN + 1;
  dgemm('T', 'N', M, M, N, ONE, U, LDU, WORK.asMatrix(N), N, ZERO,
      WORK(NNP1).asMatrix(N), N);
  for (J = 1; J <= M; J++) {
    JJ = NN + (J - 1) * N + J;
    WORK[JJ] = WORK[JJ] - D[J];
  }
  if (KBAND == 1 && N > 1) {
    for (J = 2; J <= M; J++) {
      JJ1 = NN + (J - 1) * N + J - 1;
      JJ2 = NN + (J - 2) * N + J;
      WORK[JJ1] = WORK[JJ1] - E[J - 1];
      WORK[JJ2] = WORK[JJ2] - E[J - 1];
    }
  }
  WNORM = dlansy('1', UPLO, M, WORK(NNP1).asMatrix(N), N, WORK(1));

  if (ANORM > WNORM) {
    RESULT[1] = (WNORM / ANORM) / (M * ULP);
  } else {
    if (ANORM < ONE) {
      RESULT[1] = (min(WNORM, M * ANORM) / ANORM) / (M * ULP);
    } else {
      RESULT[1] = min(WNORM / ANORM, M.toDouble()) / (M * ULP);
    }
  }

  // Do Test 2

  // Compute  U**T U - I

  if (ITYPE == 1) dort01('Columns', N, M, U, LDU, WORK, 2 * N * N, RESULT.box(2));
}
