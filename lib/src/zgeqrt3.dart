import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarfg.dart';

void zgeqrt3(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> T_,
  final int LDT,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final T = T_.dim(LDT);

  int I, I1, J, J1, N1, N2;
  final IINFO = Box(0);

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -2;
  } else if (M < N) {
    INFO.value = -1;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZGEQRT3', -INFO.value);
    return;
  }

  if (N == 1) {
    // Compute Householder transform when N=1

    zlarfg(M, A(1, 1), A(min(2, M), 1).asArray(), 1, T(1, 1));
  } else {
    // Otherwise, split A into blocks...

    N1 = N ~/ 2;
    N2 = N - N1;
    J1 = min(N1 + 1, N);
    I1 = min(N + 1, M);

    // Compute A(1:M,1:N1) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H

    zgeqrt3(M, N1, A, LDA, T, LDT, IINFO);

    // Compute A(1:M,J1:N) = Q1^H A(1:M,J1:N) [workspace: T(1:N1,J1:N)]

    for (J = 1; J <= N2; J++) {
      for (I = 1; I <= N1; I++) {
        T[I][J + N1] = A[I][J + N1];
      }
    }
    ztrmm('L', 'L', 'C', 'U', N1, N2, Complex.one, A, LDA, T(1, J1), LDT);

    zgemm('C', 'N', N1, N2, M - N1, Complex.one, A(J1, 1), LDA, A(J1, J1), LDA,
        Complex.one, T(1, J1), LDT);

    ztrmm('L', 'U', 'C', 'N', N1, N2, Complex.one, T, LDT, T(1, J1), LDT);

    zgemm('N', 'N', M - N1, N2, N1, -Complex.one, A(J1, 1), LDA, T(1, J1), LDT,
        Complex.one, A(J1, J1), LDA);

    ztrmm('L', 'L', 'N', 'U', N1, N2, Complex.one, A, LDA, T(1, J1), LDT);

    for (J = 1; J <= N2; J++) {
      for (I = 1; I <= N1; I++) {
        A[I][J + N1] = A[I][J + N1] - T[I][J + N1];
      }
    }

    // Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H

    zgeqrt3(M - N1, N2, A(J1, J1), LDA, T(J1, J1), LDT, IINFO);

    // Compute T3 = T(1:N1,J1:N) = -T1 Y1^H Y2 T2

    for (I = 1; I <= N1; I++) {
      for (J = 1; J <= N2; J++) {
        T[I][J + N1] = A[J + N1][I].conjugate();
      }
    }

    ztrmm(
        'R', 'L', 'N', 'U', N1, N2, Complex.one, A(J1, J1), LDA, T(1, J1), LDT);

    zgemm('C', 'N', N1, N2, M - N, Complex.one, A(I1, 1), LDA, A(I1, J1), LDA,
        Complex.one, T(1, J1), LDT);

    ztrmm('L', 'U', 'N', 'N', N1, N2, -Complex.one, T, LDT, T(1, J1), LDT);

    ztrmm(
        'R', 'U', 'N', 'N', N1, N2, Complex.one, T(J1, J1), LDT, T(1, J1), LDT);

    // Y = (Y1,Y2); R = [ R1  A(1:N1,J1:N) ];  T = [T1 T3]
    //                  [  0        R2     ]       [ 0 T2]
  }
}
