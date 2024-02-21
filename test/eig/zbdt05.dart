import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';

void zbdt05(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> S_,
  final int NS,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> VT_,
  final int LDVT,
  final Array<Complex> WORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final S = S_.dim();
  final U = U_.dim(LDU);
  final VT = VT_.dim(LDVT);
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  int I, J;
  double ANORM, EPS;
  final DUM = Array<double>(1);

  // Quick return if possible.

  RESID.value = ZERO;
  if (min(M, N) <= 0 || NS <= 0) return;

  EPS = dlamch('Precision');
  ANORM = zlange('M', M, N, A, LDA, DUM);

  // Compute U' * A * V.

  zgemm('N', 'C', M, NS, N, Complex.one, A, LDA, VT, LDVT, Complex.zero,
      WORK(1 + NS * NS).asMatrix(M), M);
  zgemm('C', 'N', NS, NS, M, -Complex.one, U, LDU,
      WORK(1 + NS * NS).asMatrix(M), M, Complex.zero, WORK.asMatrix(NS), NS);

  // norm(S - U' * B * V)

  J = 0;
  for (I = 1; I <= NS; I++) {
    // 10
    WORK[J + I] = WORK[J + I] + Complex(S[I], ZERO);
    RESID.value = max(RESID.value, dzasum(NS, WORK(J + 1), 1));
    J = J + NS;
  } // 10

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    if (ANORM >= RESID.value) {
      RESID.value = (RESID.value / ANORM) / (N.toDouble() * EPS);
    } else {
      if (ANORM < ONE) {
        RESID.value = (min(RESID.value, (N).toDouble() * ANORM) / ANORM) /
            (N.toDouble() * EPS);
      } else {
        RESID.value =
            min(RESID.value / ANORM, (N).toDouble()) / (N.toDouble() * EPS);
      }
    }
  }
}
