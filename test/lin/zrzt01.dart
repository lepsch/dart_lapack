import 'dart:math';

import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zunmrz.dart';

import 'xerbla.dart';

double zrzt01(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final Matrix<Complex> AF_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AF = AF_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having(length: LWORK);
  const ZERO = 0.0, ONE = 1.0;
  final RWORK = Array<double>(1);
  final INFO = Box(0);

  if (LWORK < M * N + M) {
    xerbla('ZRZT01', 8);
    return ZERO;
  }

  // Quick return if possible

  if (M <= 0 || N <= 0) return ZERO;

  final NORMA = zlange('One-norm', M, N, A, LDA, RWORK);

  // Copy upper triangle R

  zlaset('Full', M, N, Complex.zero, Complex.zero, WORK.asMatrix(), M);
  for (var J = 1; J <= M; J++) {
    for (var I = 1; I <= J; I++) {
      WORK[(J - 1) * M + I] = AF[I][J];
    }
  }

  // R = R * P(1) * ... *P(m)

  zunmrz('Right', 'No transpose', M, N, M, N - M, AF, LDA, TAU, WORK.asMatrix(),
      M, WORK(M * N + 1), LWORK - M * N, INFO);

  // R = R - A

  for (var I = 1; I <= N; I++) {
    zaxpy(M, Complex(-ONE), A(1, I).asArray(), 1, WORK((I - 1) * M + 1), 1);
  }

  var result = zlange('One-norm', M, N, WORK.asMatrix(), M, RWORK);

  result = result / (dlamch('Epsilon') * max(M, N));
  return NORMA != ZERO ? result / NORMA : result;
}
