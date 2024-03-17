import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtpqrt2(
  final int M,
  final int N,
  final int L,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> T_,
  final int LDT,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final T = T_.having(ld: LDT);
  const ONE = 1.0, ZERO = 0.0;
  int I, J, P, MP, NP;
  double ALPHA;

  // Test the input arguments

  INFO.value = 0;
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (L < 0 || L > min(M, N)) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, M)) {
    INFO.value = -7;
  } else if (LDT < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DTPQRT2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || M == 0) return;

  for (I = 1; I <= N; I++) {
    // Generate elementary reflector H(I) to annihilate B(:,I)

    P = M - L + min(L, I);
    dlarfg(P + 1, A.box(I, I), B(1, I).asArray(), 1, T.box(I, 1));
    if (I < N) {
      // W(1:N-I) := C(I:M,I+1:N)^H * C(I:M,I) [use W = T(:,N)]

      for (J = 1; J <= N - I; J++) {
        T[J][N] = (A[I][I + J]);
      }
      dgemv('T', P, N - I, ONE, B(1, I + 1), LDB, B(1, I).asArray(), 1, ONE,
          T(1, N).asArray(), 1);

      // C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)^H

      ALPHA = -(T[I][1]);
      for (J = 1; J <= N - I; J++) {
        A[I][I + J] += ALPHA * (T[J][N]);
      }
      dger(P, N - I, ALPHA, B(1, I).asArray(), 1, T(1, N).asArray(), 1,
          B(1, I + 1), LDB);
    }
  }

  for (I = 2; I <= N; I++) {
    // T(1:I-1,I) := C(I:M,1:I-1)^H * (alpha * C(I:M,I))

    ALPHA = -T[I][1];

    for (J = 1; J <= I - 1; J++) {
      T[J][I] = ZERO;
    }
    P = min(I - 1, L);
    MP = min(M - L + 1, M);
    NP = min(P + 1, N);

    // Triangular part of B2

    for (J = 1; J <= P; J++) {
      T[J][I] = ALPHA * B[M - L + J][I];
    }
    dtrmv('U', 'T', 'N', P, B(MP, 1), LDB, T(1, I).asArray(), 1);

    // Rectangular part of B2

    dgemv('T', L, I - 1 - P, ALPHA, B(MP, NP), LDB, B(MP, I).asArray(), 1, ZERO,
        T(NP, I).asArray(), 1);

    // B1

    dgemv('T', M - L, I - 1, ALPHA, B, LDB, B(1, I).asArray(), 1, ONE,
        T(1, I).asArray(), 1);

    // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

    dtrmv('U', 'N', 'N', I - 1, T, LDT, T(1, I).asArray(), 1);

    // T(I,I) = tau(I)

    T[I][I] = T[I][1];
    T[I][1] = ZERO;
  }
}
