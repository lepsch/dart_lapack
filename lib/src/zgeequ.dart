import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgeequ(
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
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final R = R_.having();
  final C = C_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, J;
  double BIGNUM, RCMAX, RCMIN, SMLNUM;

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
    xerbla('ZGEEQU', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    ROWCND.value = ONE;
    COLCND.value = ONE;
    AMAX.value = ZERO;
    return;
  }

  // Get machine constants.

  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;

  // Compute row scale factors.

  for (I = 1; I <= M; I++) {
    R[I] = ZERO;
  }

  // Find the maximum element in each row.

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      R[I] = max(R[I], A[I][J].cabs1());
    }
  }

  // Find the maximum and minimum scale factors.

  RCMIN = BIGNUM;
  RCMAX = ZERO;
  for (I = 1; I <= M; I++) {
    RCMAX = max(RCMAX, R[I]);
    RCMIN = min(RCMIN, R[I]);
  }
  AMAX.value = RCMAX;

  if (RCMIN == ZERO) {
    // Find the first zero scale factor and return an error code.

    for (I = 1; I <= M; I++) {
      if (R[I] == ZERO) {
        INFO.value = I;
        return;
      }
    }
  } else {
    // Invert the scale factors.

    for (I = 1; I <= M; I++) {
      R[I] = ONE / min(max(R[I], SMLNUM), BIGNUM);
    }

    // Compute ROWCND = min(R(I)) / max(R(I))

    ROWCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
  }

  // Compute column scale factors

  for (J = 1; J <= N; J++) {
    C[J] = ZERO;
  }

  // Find the maximum element in each column,
  // assuming the row scaling computed above.

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      C[J] = max(C[J], A[I][J].cabs1() * R[I]);
    }
  }

  // Find the maximum and minimum scale factors.

  RCMIN = BIGNUM;
  RCMAX = ZERO;
  for (J = 1; J <= N; J++) {
    RCMIN = min(RCMIN, C[J]);
    RCMAX = max(RCMAX, C[J]);
  }

  if (RCMIN == ZERO) {
    // Find the first zero scale factor and return an error code.

    for (J = 1; J <= N; J++) {
      if (C[J] == ZERO) {
        INFO.value = M + J;
        return;
      }
    }
  } else {
    // Invert the scale factors.

    for (J = 1; J <= N; J++) {
      C[J] = ONE / min(max(C[J], SMLNUM), BIGNUM);
    }

    // Compute COLCND = min(C(J)) / max(C(J))

    COLCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
  }
}
