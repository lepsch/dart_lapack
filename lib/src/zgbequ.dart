import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zgbequ(
  final int M,
  final int N,
  final int KL,
  final int KU,
  final Matrix<Complex> AB_,
  final int LDAB,
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
  final AB = AB_.dim(LDAB);
  final C = C_.dim();
  final R = R_.dim();
  const ONE = 1.0, ZERO = 0.0;
  int I, J, KD;
  double BIGNUM, RCMAX, RCMIN, SMLNUM;

  double CABS1(Complex ZDUM) => ZDUM.toDouble().abs() + ZDUM.imaginary.abs();

  // Test the input parameters

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (KL < 0) {
    INFO.value = -3;
  } else if (KU < 0) {
    INFO.value = -4;
  } else if (LDAB < KL + KU + 1) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZGBEQU', -INFO.value);
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
    // 10
    R[I] = ZERO;
  } // 10

  // Find the maximum element in each row.

  KD = KU + 1;
  for (J = 1; J <= N; J++) {
    // 30
    for (I = max(J - KU, 1); I <= min(J + KL, M); I++) {
      // 20
      R[I] = max(R[I], CABS1(AB[KD + I - J][J]));
    } // 20
  } // 30

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

    // Compute ROWCND.value = min(R(I)) / max(R(I))

    ROWCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
  }

  // Compute column scale factors

  for (J = 1; J <= N; J++) {
    // 70
    C[J] = ZERO;
  } // 70

  // Find the maximum element in each column,
  // assuming the row scaling computed above.

  KD = KU + 1;
  for (J = 1; J <= N; J++) {
    // 90
    for (I = max(J - KU, 1); I <= min(J + KL, M); I++) {
      // 80
      C[J] = max(C[J], CABS1(AB[KD + I - J][J]) * R[I]);
    } // 80
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

    // Compute COLCND.value = min(C(J)) / max(C(J))

    COLCND.value = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
  }
}
