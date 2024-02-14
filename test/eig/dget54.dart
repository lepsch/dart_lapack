import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget54(
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> S_,
  final int LDS,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> V_,
  final int LDV,
  final Array<double> WORK_,
  final Box<double> RESULT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final S = S_.dim(LDS);
  final T = T_.dim(LDT);
  final U = U_.dim(LDU);
  final V = V_.dim(LDV);
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  double ABNORM, ULP, UNFL, WNORM;
  final DUM = Array<double>(1);

  RESULT.value = ZERO;
  if (N <= 0) return;

  // Constants

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');

  // compute the norm of (A,B)

  dlacpy('Full', N, N, A, LDA, WORK.asMatrix(N), N);
  dlacpy('Full', N, N, B, LDB, WORK(N * N + 1).asMatrix(N), N);
  ABNORM = max(dlange('1', N, 2 * N, WORK.asMatrix(N), N, DUM), UNFL);

  // Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)

  dlacpy(' ', N, N, A, LDA, WORK.asMatrix(N), N);
  dgemm('N', 'N', N, N, N, ONE, U, LDU, S, LDS, ZERO,
      WORK(N * N + 1).asMatrix(N), N);

  dgemm('N', 'C', N, N, N, -ONE, WORK(N * N + 1).asMatrix(N), N, V, LDV, ONE,
      WORK.asMatrix(N), N);

  // Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)

  dlacpy(' ', N, N, B, LDB, WORK(N * N + 1).asMatrix(N), N);
  dgemm('N', 'N', N, N, N, ONE, U, LDU, T, LDT, ZERO,
      WORK(2 * N * N + 1).asMatrix(N), N);

  dgemm('N', 'C', N, N, N, -ONE, WORK(2 * N * N + 1).asMatrix(N), N, V, LDV,
      ONE, WORK(N * N + 1).asMatrix(N), N);

  // Compute norm(W)/ ( ulp*norm((A,B)) )

  WNORM = dlange('1', N, 2 * N, WORK.asMatrix(N), N, DUM);

  if (ABNORM > WNORM) {
    RESULT.value = (WNORM / ABNORM) / (2 * N * ULP);
  } else {
    if (ABNORM < ONE) {
      RESULT.value = (min(WNORM, 2 * N * ABNORM) / ABNORM) / (2 * N * ULP);
    } else {
      RESULT.value = min(WNORM / ABNORM, (2 * N).toDouble()) / (2 * N * ULP);
    }
  }
}
