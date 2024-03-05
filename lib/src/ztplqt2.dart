import 'dart:math';

import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgerc.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlarfg.dart';

void ztplqt2(
  final int M,
  final int N,
  final int L,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> T_,
  final int LDT,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final T = T_.having(ld: LDT);

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int I, J, P, MP, NP;
  Complex ALPHA;

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
    xerbla('ZTPLQT2', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || M == 0) return;

  for (I = 1; I <= M; I++) {
    // Generate elementary reflector H(I) to annihilate B(I,:)

    P = N - L + min(L, I);
    zlarfg(P + 1, A(I, I), B(I, 1).asArray(), LDB, T(1, I));
    T[1][I] = T[1][I].conjugate();
    if (I < M) {
      for (J = 1; J <= P; J++) {
        B[I][J] = B[I][J].conjugate();
      }

      // W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]

      for (J = 1; J <= M - I; J++) {
        T[M][J] = A[I + J][I];
      }
      zgemv('N', M - I, P, Complex.one, B(I + 1, 1), LDB, B(I, 1).asArray(),
          LDB, Complex.one, T(M, 1).asArray(), LDT);

      // C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H

      ALPHA = -T[1][I];
      for (J = 1; J <= M - I; J++) {
        A[I + J][I] = A[I + J][I] + ALPHA * T[M][J];
      }
      zgerc(M - I, P, (ALPHA), T(M, 1).asArray(), LDT, B(I, 1).asArray(), LDB,
          B(I + 1, 1), LDB);
      for (J = 1; J <= P; J++) {
        B[I][J] = B[I][J].conjugate();
      }
    }
  }

  for (I = 2; I <= M; I++) {
    // T(I,1:I-1) := C(I:I-1,1:N)**H * (alpha * C(I,I:N))

    ALPHA = -(T[1][I]);
    for (J = 1; J <= I - 1; J++) {
      T[I][J] = Complex.zero;
    }
    P = min(I - 1, L);
    NP = min(N - L + 1, N);
    MP = min(P + 1, M);
    for (J = 1; J <= N - L + P; J++) {
      B[I][J] = B[I][J].conjugate();
    }

    // Triangular part of B2

    for (J = 1; J <= P; J++) {
      T[I][J] = (ALPHA * B[I][N - L + J]);
    }
    ztrmv('L', 'N', 'N', P, B(1, NP), LDB, T(I, 1).asArray(), LDT);

    // Rectangular part of B2

    zgemv('N', I - 1 - P, L, ALPHA, B(MP, NP), LDB, B(I, NP).asArray(), LDB,
        Complex.zero, T(I, MP).asArray(), LDT);

    // B1

    zgemv('N', I - 1, N - L, ALPHA, B, LDB, B(I, 1).asArray(), LDB, Complex.one,
        T(I, 1).asArray(), LDT);

    // T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)

    for (J = 1; J <= I - 1; J++) {
      T[I][J] = T[I][J].conjugate();
    }
    ztrmv('L', 'C', 'N', I - 1, T, LDT, T(I, 1).asArray(), LDT);
    for (J = 1; J <= I - 1; J++) {
      T[I][J] = T[I][J].conjugate();
    }
    for (J = 1; J <= N - L + P; J++) {
      B[I][J] = B[I][J].conjugate();
    }

    // T(I,I) = tau(I)

    T[I][I] = T[1][I];
    T[1][I] = Complex.zero;
  }
  for (I = 1; I <= M; I++) {
    for (J = I + 1; J <= M; J++) {
      T[I][J] = T[J][I];
      T[J][I] = Complex.zero;
    }
  }
}
