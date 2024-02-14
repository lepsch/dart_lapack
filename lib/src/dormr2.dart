import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dormr2(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final TAU = TAU_.dim();
  final C = C_.dim(LDC);
  final WORK = WORK_.dim();
  const ONE = 1.0;
  bool LEFT, NOTRAN;
  int I, I1, I2, I3, MI = 0, NI = 0, NQ;
  double AII;

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');

  // NQ is the order of Q

  if (LEFT) {
    NQ = M;
  } else {
    NQ = N;
  }
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > NQ) {
    INFO.value = -5;
  } else if (LDA < max(1, K)) {
    INFO.value = -7;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DORMR2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0 || K == 0) return;

  if ((LEFT && !NOTRAN) || (!LEFT && NOTRAN)) {
    I1 = 1;
    I2 = K;
    I3 = 1;
  } else {
    I1 = K;
    I2 = 1;
    I3 = -1;
  }

  if (LEFT) {
    NI = N;
  } else {
    MI = M;
  }

  for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
    if (LEFT) {
      // H(i) is applied to C(1:m-k+i,1:n)

      MI = M - K + I;
    } else {
      // H(i) is applied to C(1:m,1:n-k+i)

      NI = N - K + I;
    }

    // Apply H(i)

    AII = A[I][NQ - K + I];
    A[I][NQ - K + I] = ONE;
    dlarf(SIDE, MI, NI, A(I, 1).asArray(), LDA, TAU[I], C, LDC, WORK);
    A[I][NQ - K + I] = AII;
  }
}
