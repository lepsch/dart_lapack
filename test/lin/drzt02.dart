import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dormrz.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

double drzt02(
  final int M,
  final int N,
  final Matrix<double> AF_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AF = AF_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);
  final INFO = Box(0);

  if (LWORK < N * N + N) {
    xerbla('DRZT02', 7);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0 || N <= 0) return ZERO;

  // Q := I

  dlaset('Full', N, N, ZERO, ONE, WORK.asMatrix(), N);

  // Q := P(1) * ... * P(m) * Q

  dormrz('Left', 'No transpose', N, N, M, N - M, AF, LDA, TAU, WORK.asMatrix(),
      N, WORK(N * N + 1), LWORK - N * N, INFO);

  // Q := P(m) * ... * P(1) * Q

  dormrz('Left', 'Transpose', N, N, M, N - M, AF, LDA, TAU, WORK.asMatrix(), N,
      WORK(N * N + 1), LWORK - N * N, INFO);

  // Q := Q - I

  for (var I = 1; I <= N; I++) {
    WORK[(I - 1) * N + I] -= ONE;
  }

  return dlange('One-norm', N, N, WORK.asMatrix(), N, RWORK) /
      (dlamch('Epsilon') * max(M, N));
}
