import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zggsvd3.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlanhe.dart';
import 'package:lapack/src/zlaset.dart';

void zgsvts3(
  final int M,
  final int P,
  final int N,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> BF_,
  final int LDB,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<double> ALPHA_,
  final Array<double> BETA_,
  final Matrix<Complex> R_,
  final int LDR,
  final Array<int> IWORK_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final AF = AF_.dim(LDA);
  final B = B_.dim(LDB);
  final BF = BF_.dim(LDB);
  final U = U_.dim(LDU);
  final V = V_.dim(LDV);
  final Q = Q_.dim(LDQ);
  final R = R_.dim(LDR);
  final IWORK = IWORK_.dim();
  final WORK = WORK_.dim(LWORK);
  final RWORK = RWORK_.dim();
  final ALPHA = ALPHA_.dim();
  final BETA = BETA_.dim();
  final RESULT = RESULT_.dim();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL;
  final INFO = Box(0), K = Box(0), L = Box(0);

  ULP = dlamch('Precision');
  ULPINV = ONE / ULP;
  UNFL = dlamch('Safe minimum');

  // Copy the matrix A to the array AF.

  zlacpy('Full', M, N, A, LDA, AF, LDA);
  zlacpy('Full', P, N, B, LDB, BF, LDB);

  ANORM = max(zlange('1', M, N, A, LDA, RWORK), UNFL);
  BNORM = max(zlange('1', P, N, B, LDB, RWORK), UNFL);

  // Factorize the matrices A and B in the arrays AF and BF.

  zggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU,
      V, LDV, Q, LDQ, WORK, LWORK, RWORK, IWORK, INFO);

  // Copy R

  for (I = 1; I <= min(K.value + L.value, M); I++) {
    // 20
    for (J = I; J <= K.value + L.value; J++) {
      // 10
      R[I][J] = AF[I][N - K.value - L.value + J];
    } // 10
  } // 20

  if (M - K.value - L.value < 0) {
    for (I = M + 1; I <= K.value + L.value; I++) {
      // 40
      for (J = I; J <= K.value + L.value; J++) {
        // 30
        R[I][J] = BF[I - K.value][N - K.value - L.value + J];
      } // 30
    } // 40
  }

  // Compute A:= U'*A*Q - D1*R

  zgemm('No transpose', 'No transpose', M, N, N, Complex.one, A, LDA, Q, LDQ,
      Complex.zero, WORK.asMatrix(), LDA);

  zgemm('Conjugate transpose', 'No transpose', M, N, M, Complex.one, U, LDU,
      WORK.asMatrix(), LDA, Complex.zero, A, LDA);

  for (I = 1; I <= K.value; I++) {
    // 60
    for (J = I; J <= K.value + L.value; J++) {
      // 50
      A[I][N - K.value - L.value + J] =
          A[I][N - K.value - L.value + J] - R[I][J];
    } // 50
  } // 60

  for (I = K.value + 1; I <= min(K.value + L.value, M); I++) {
    // 80
    for (J = I; J <= K.value + L.value; J++) {
      // 70
      A[I][N - K.value - L.value + J] =
          A[I][N - K.value - L.value + J] - ALPHA[I].toComplex() * R[I][J];
    } // 70
  } // 80

  // Compute norm( U'*A*Q - D1*R ) / ( max(1,M,N)*norm(A)*ULP ) .

  RESID = zlange('1', M, N, A, LDA, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = ((RESID / (max(1, max(M, N))).toDouble()) / ANORM) / ULP;
  } else {
    RESULT[1] = ZERO;
  }

  // Compute B := V'*B*Q - D2*R

  zgemm('No transpose', 'No transpose', P, N, N, Complex.one, B, LDB, Q, LDQ,
      Complex.zero, WORK.asMatrix(), LDB);

  zgemm('Conjugate transpose', 'No transpose', P, N, P, Complex.one, V, LDV,
      WORK.asMatrix(), LDB, Complex.zero, B, LDB);

  for (I = 1; I <= L.value; I++) {
    // 100
    for (J = I; J <= L.value; J++) {
      // 90
      B[I][N - L.value + J] = B[I][N - L.value + J] -
          BETA[K.value + I].toComplex() * R[K.value + I][K.value + J];
    } // 90
  } // 100

  // Compute norm( V'*B*Q - D2*R ) / ( max(P,N)*norm(B)*ULP ) .

  RESID = zlange('1', P, N, B, LDB, RWORK);
  if (BNORM > ZERO) {
    RESULT[2] = ((RESID / (max(1, max(P, N))).toDouble()) / BNORM) / ULP;
  } else {
    RESULT[2] = ZERO;
  }

  // Compute I - U'*U

  zlaset('Full', M, M, Complex.zero, Complex.one, WORK.asMatrix(), LDQ);
  zherk('Upper', 'Conjugate transpose', M, M, -ONE, U, LDU, ONE,
      WORK.asMatrix(), LDU);

  // Compute norm( I - U'*U ) / ( M * ULP ) .

  RESID = zlanhe('1', 'Upper', M, WORK.asMatrix(), LDU, RWORK);
  RESULT[3] = (RESID / (max(1, M)).toDouble()) / ULP;

  // Compute I - V'*V

  zlaset('Full', P, P, Complex.zero, Complex.one, WORK.asMatrix(), LDV);
  zherk('Upper', 'Conjugate transpose', P, P, -ONE, V, LDV, ONE,
      WORK.asMatrix(), LDV);

  // Compute norm( I - V'*V ) / ( P * ULP ) .

  RESID = zlanhe('1', 'Upper', P, WORK.asMatrix(), LDV, RWORK);
  RESULT[4] = (RESID / (max(1, P)).toDouble()) / ULP;

  // Compute I - Q'*Q

  zlaset('Full', N, N, Complex.zero, Complex.one, WORK.asMatrix(), LDQ);
  zherk('Upper', 'Conjugate transpose', N, N, -ONE, Q, LDQ, ONE,
      WORK.asMatrix(), LDQ);

  // Compute norm( I - Q'*Q ) / ( N * ULP ) .

  RESID = zlanhe('1', 'Upper', N, WORK.asMatrix(), LDQ, RWORK);
  RESULT[5] = (RESID / (max(1, N)).toDouble()) / ULP;

  // Check sorting

  dcopy(N, ALPHA, 1, RWORK, 1);
  for (I = K.value + 1; I <= min(K.value + L.value, M); I++) {
    // 110
    J = IWORK[I];
    if (I != J) {
      TEMP = RWORK[I];
      RWORK[I] = RWORK[J];
      RWORK[J] = TEMP;
    }
  } // 110

  RESULT[6] = ZERO;
  for (I = K.value + 1; I <= min(K.value + L.value, M) - 1; I++) {
    // 120
    if (RWORK[I] < RWORK[I + 1]) RESULT[6] = ULPINV;
  } // 120
}
