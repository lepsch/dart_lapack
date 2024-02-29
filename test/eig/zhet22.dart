import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zhemm.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlanhe.dart';

import 'zunt01.dart';

void zhet22(
  final int ITYPE,
  final String UPLO,
  final int N,
  final int M,
  final int KBAND,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
  final A = A_.dim(LDA);
  final U = U_.dim(LDU);
  // final V = V_.dim(LDV);
  // final TAU = TAU_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final D = D_.dim();
  final E = E_.dim();
  final RESULT = RESULT_.dim();
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

  ANORM = max(zlanhe('1', UPLO, N, A, LDA, RWORK), UNFL);

  // Compute error matrix:

  // ITYPE=1: error = U**H A U - S

  zhemm('L', UPLO, N, M, Complex.one, A, LDA, U, LDU, Complex.zero,
      WORK.asMatrix(), N);
  NN = N * N;
  NNP1 = NN + 1;
  zgemm('C', 'N', M, M, N, Complex.one, U, LDU, WORK.asMatrix(), N,
      Complex.zero, WORK(NNP1).asMatrix(), N);
  for (J = 1; J <= M; J++) {
    // 10
    JJ = NN + (J - 1) * N + J;
    WORK[JJ] = WORK[JJ] - D[J].toComplex();
  } // 10
  if (KBAND == 1 && N > 1) {
    for (J = 2; J <= M; J++) {
      // 20
      JJ1 = NN + (J - 1) * N + J - 1;
      JJ2 = NN + (J - 2) * N + J;
      WORK[JJ1] = WORK[JJ1] - E[J - 1].toComplex();
      WORK[JJ2] = WORK[JJ2] - E[J - 1].toComplex();
    } // 20
  }
  WNORM = zlanhe('1', UPLO, M, WORK(NNP1).asMatrix(), N, RWORK);

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

  // Compute  U**H U - I

  if (ITYPE == 1)
    zunt01('Columns', N, M, U, LDU, WORK, 2 * N * N, RWORK, RESULT(2));
}
