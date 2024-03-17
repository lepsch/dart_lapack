import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlansy.dart';

void zstt22(
  final int N,
  final int M,
  final int KBAND,
  final Array<double> AD_,
  final Array<double> AE_,
  final Array<double> SD_,
  final Array<double> SE_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> WORK_,
  final int LDWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having(ld: LDWORK);
  final RWORK = RWORK_.having();
  final AD = AD_.having();
  final AE = AE_.having();
  final SD = SD_.having();
  final SE = SE_.having();
  final RESULT = RESULT_.having(length: 2);
  const ZERO = 0.0, ONE = 1.0;
  int I, J, K;
  double ANORM, ULP, UNFL, WNORM;
  Complex AUKJ;

  RESULT[1] = ZERO;
  RESULT[2] = ZERO;
  if (N <= 0 || M <= 0) return;

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon');

  // Do Test 1

  // Compute the 1-norm of A.

  if (N > 1) {
    ANORM = AD[1].abs() + AE[1].abs();
    for (J = 2; J <= N - 1; J++) {
      // 10
      ANORM = max(ANORM, AD[J].abs() + AE[J].abs() + AE[J - 1].abs());
    } // 10
    ANORM = max(ANORM, AD[N].abs() + AE[N - 1].abs());
  } else {
    ANORM = (AD[1]).abs();
  }
  ANORM = max(ANORM, UNFL);

  // Norm of U*AU - S

  for (I = 1; I <= M; I++) {
    // 40
    for (J = 1; J <= M; J++) {
      // 30
      WORK[I][J] = Complex.zero;
      for (K = 1; K <= N; K++) {
        // 20
        AUKJ = AD[K].toComplex() * U[K][J];
        if (K != N) AUKJ = AUKJ + AE[K].toComplex() * U[K + 1][J];
        if (K != 1) AUKJ = AUKJ + AE[K - 1].toComplex() * U[K - 1][J];
        WORK[I][J] += U[K][I] * AUKJ;
      } // 20
    } // 30
    WORK[I][I] -= SD[I].toComplex();
    if (KBAND == 1) {
      if (I != 1) WORK[I][I - 1] -= SE[I - 1].toComplex();
      if (I != N) WORK[I][I + 1] -= SE[I].toComplex();
    }
  } // 40

  WNORM = zlansy('1', 'L', M, WORK, M, RWORK);

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

  // Compute  U*U - I

  zgemm('T', 'N', M, M, N, Complex.one, U, LDU, U, LDU, Complex.zero, WORK, M);

  for (J = 1; J <= M; J++) {
    // 50
    WORK[J][J] -= Complex.one;
  } // 50

  RESULT[2] = min(M.toDouble(), zlange('1', M, M, WORK, M, RWORK)) / (M * ULP);
}
