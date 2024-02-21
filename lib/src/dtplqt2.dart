import 'dart:math';

import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtplqt2(
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
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final T = T_.dim(LDT);
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
  } else if (LDA < max(1, M)) {
    INFO.value = -5;
  } else if (LDB < max(1, M)) {
    INFO.value = -7;
  } else if (LDT < max(1, M)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DTPLQT2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || M == 0) return;

  for (I = 1; I <= M; I++) {
    // Generate elementary reflector H(I) to annihilate B(I,:)

    P = N - L + min(L, I);
    dlarfg(P + 1, A.box(I, I), B(I, 1).asArray(), LDB, T.box(1, I));
    if (I < M) {
      // W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]

      for (J = 1; J <= M - I; J++) {
        T[M][J] = A[I + J][I];
      }
      dgemv('N', M - I, P, ONE, B(I + 1, 1), LDB, B(I, 1).asArray(), LDB, ONE,
          T(M, 1).asArray(), LDT);

      // C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H

      ALPHA = -(T[1][I]);
      for (J = 1; J <= M - I; J++) {
        A[I + J][I] = A[I + J][I] + ALPHA * (T[M][J]);
      }
      dger(M - I, P, ALPHA, T(M, 1).asArray(), LDT, B(I, 1).asArray(), LDB,
          B(I + 1, 1), LDB);
    }
  }

  for (I = 2; I <= M; I++) {
    // T(I,1:I-1) := C(I:I-1,1:N) * (alpha * C(I,I:N)^H)

    ALPHA = -T[1][I];

    for (J = 1; J <= I - 1; J++) {
      T[I][J] = ZERO;
    }
    P = min(I - 1, L);
    NP = min(N - L + 1, N);
    MP = min(P + 1, M);

    // Triangular part of B2

    for (J = 1; J <= P; J++) {
      T[I][J] = ALPHA * B[I][N - L + J];
    }
    dtrmv('L', 'N', 'N', P, B(1, NP), LDB, T(I, 1).asArray(), LDT);

    // Rectangular part of B2

    dgemv('N', I - 1 - P, L, ALPHA, B(MP, NP), LDB, B(I, NP).asArray(), LDB,
        ZERO, T(I, MP).asArray(), LDT);

    // B1

    dgemv('N', I - 1, N - L, ALPHA, B, LDB, B(I, 1).asArray(), LDB, ONE,
        T(I, 1).asArray(), LDT);

    // T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)

    dtrmv('L', 'T', 'N', I - 1, T, LDT, T(I, 1).asArray(), LDT);

    // T(I,I) = tau(I)

    T[I][I] = T[1][I];
    T[1][I] = ZERO;
  }
  for (I = 1; I <= M; I++) {
    for (J = I + 1; J <= M; J++) {
      T[I][J] = T[J][I];
      T[J][I] = ZERO;
    }
  }
}
