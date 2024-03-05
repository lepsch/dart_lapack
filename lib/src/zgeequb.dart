import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgeequb(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> R_,
  final Array<double> C_,
  final Box<double> ROWCND,
  final Box<double> COLCND,
  final Box<double> AMAX,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final R = R_.having();
  final C = C_.having();
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  double BIGNUM, RCMAX, RCMIN, SMLNUM, RADIX, LOGRDX;

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  // Test the input parameters.

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  }
  if (INFO.value != 0) {
    xerbla('ZGEEQUB', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (M == 0 || N == 0) {
    ROWCND.value = ONE;
    COLCND.value = ONE;
    AMAX.value = ZERO;
    return;
  }

  // Get machine constants.  Assume SMLNUM is a power of the radix.

  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;
  RADIX = dlamch('B');
  LOGRDX = log(RADIX);

  // Compute row scale factors.

  for (I = 1; I <= M; I++) {
    // 10
    R[I] = ZERO;
  } // 10

  // Find the maximum element in each row.

  for (J = 1; J <= N; J++) {
    // 30
    for (I = 1; I <= M; I++) {
      // 20
      R[I] = max(R[I], CABS1(A[I][J]));
    } // 20
  } // 30
  for (I = 1; I <= M; I++) {
    if (R[I] > ZERO) {
      R[I] = pow(RADIX, log(R[I]) ~/ LOGRDX).toDouble();
    }
  }

  // Find the maximum and minimum scale factors.

  RCMIN = BIGNUM;
  RCMAX = ZERO;
  for (I = 1; I <= M; I++) {
    // 40
    RCMAX = max(RCMAX, R[I]);
    RCMIN = min(RCMIN, R[I]);
  } // 40
  AMAX.value = RCMAX;

  if (RCMIN == ZERO) {
    // Find the first zero scale factor and return an error code.

    for (I = 1; I <= M; I++) {
      // 50
      if (R[I] == ZERO) {
        INFO.value = I;
        return;
      }
    } // 50
  } else {
    // Invert the scale factors.

    for (I = 1; I <= M; I++) {
      // 60
      R[I] = ONE / min(max(R[I], SMLNUM), BIGNUM);
    } // 60

    // Compute ROWCND.value = min(R(I)) / max(R(I)).

    ROWCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
  }

  // Compute column scale factors.

  for (J = 1; J <= N; J++) {
    // 70
    C[J] = ZERO;
  } // 70

  // Find the maximum element in each column,
  // assuming the row scaling computed above.

  for (J = 1; J <= N; J++) {
    // 90
    for (I = 1; I <= M; I++) {
      // 80
      C[J] = max(C[J], CABS1(A[I][J]) * R[I]);
    } // 80
    if (C[J] > ZERO) {
      C[J] = pow(RADIX, log(C[J]) ~/ LOGRDX).toDouble();
    }
  } // 90

  // Find the maximum and minimum scale factors.

  RCMIN = BIGNUM;
  RCMAX = ZERO;
  for (J = 1; J <= N; J++) {
    // 100
    RCMIN = min(RCMIN, C[J]);
    RCMAX = max(RCMAX, C[J]);
  } // 100

  if (RCMIN == ZERO) {
    // Find the first zero scale factor and return an error code.

    for (J = 1; J <= N; J++) {
      // 110
      if (C[J] == ZERO) {
        INFO.value = M + J;
        return;
      }
    } // 110
  } else {
    // Invert the scale factors.

    for (J = 1; J <= N; J++) {
      // 120
      C[J] = ONE / min(max(C[J], SMLNUM), BIGNUM);
    } // 120

    // Compute COLCND.value = min(C(J)) / max(C(J)).

    COLCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
  }
}
