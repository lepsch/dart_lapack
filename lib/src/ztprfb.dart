// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/ztrmm.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void ztprfb(
  final String SIDE,
  final String TRANS,
  final String DIRECT,
  final String STOREV,
  final int M,
  final int N,
  final int K,
  final int L,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> WORK_,
  final int LDWORK,
) {
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final WORK = WORK_.having(ld: LDWORK);
  int I, J, MP, NP, KP;
  bool LEFT, FORWARD, COLUMN, RIGHT, BACKWARD, ROW;

  // Quick return if possible
  if (M <= 0 || N <= 0 || K <= 0 || L < 0) return;

  if (lsame(STOREV, 'C')) {
    COLUMN = true;
    ROW = false;
  } else if (lsame(STOREV, 'R')) {
    COLUMN = false;
    ROW = true;
  } else {
    COLUMN = false;
    ROW = false;
  }

  if (lsame(SIDE, 'L')) {
    LEFT = true;
    RIGHT = false;
  } else if (lsame(SIDE, 'R')) {
    LEFT = false;
    RIGHT = true;
  } else {
    LEFT = false;
    RIGHT = false;
  }

  if (lsame(DIRECT, 'F')) {
    FORWARD = true;
    BACKWARD = false;
  } else if (lsame(DIRECT, 'B')) {
    FORWARD = false;
    BACKWARD = true;
  } else {
    FORWARD = false;
    BACKWARD = false;
  }

  if (COLUMN && FORWARD && LEFT) {
    // Let  W =  [ I ]    (K-by-K)
    //           [ V ]    (M-by-K)
    //
    // Form  H C  or  H**H C  where  C = [ A ]  (K-by-N)
    //                                   [ B ]  (M-by-N)
    //
    // H = I - W T W**H          or  H**H = I - W T**H W**H
    //
    // A -=   T (A + V**H B)  or  A -=   T**H (A + V**H B)
    // B -= V T (A + V**H B)  or  B -= V T**H (A + V**H B)

    MP = min(M - L + 1, M);
    KP = min(L + 1, K);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        WORK[I][J] = B[M - L + I][J];
      }
    }
    ztrmm('L', 'U', 'C', 'N', L, N, Complex.one, V(MP, 1), LDV, WORK, LDWORK);
    zgemm('C', 'N', L, N, M - L, Complex.one, V, LDV, B, LDB, Complex.one, WORK,
        LDWORK);
    zgemm('C', 'N', K - L, N, M, Complex.one, V(1, KP), LDV, B, LDB,
        Complex.zero, WORK(KP, 1), LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('L', 'U', TRANS, 'N', K, N, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('N', 'N', M - L, N, K, -Complex.one, V, LDV, WORK, LDWORK,
        Complex.one, B, LDB);
    zgemm('N', 'N', L, N, K - L, -Complex.one, V(MP, KP), LDV, WORK(KP, 1),
        LDWORK, Complex.one, B(MP, 1), LDB);
    ztrmm('L', 'U', 'N', 'N', L, N, Complex.one, V(MP, 1), LDV, WORK, LDWORK);
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        B[M - L + I][J] -= WORK[I][J];
      }
    }
  } else if (COLUMN && FORWARD && RIGHT) {
    // Let  W =  [ I ]    (K-by-K)
    //           [ V ]    (N-by-K)
    //
    // Form  C H or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N)
    //
    // H = I - W T W**H          or  H**H = I - W T**H W**H
    //
    // A -= (A + B V) T      or  A -= (A + B V) T**H
    // B -= (A + B V) T V**H  or  B -= (A + B V) T**H V**H

    NP = min(N - L + 1, N);
    KP = min(L + 1, K);

    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][J] = B[I][N - L + J];
      }
    }
    ztrmm('R', 'U', 'N', 'N', M, L, Complex.one, V(NP, 1), LDV, WORK, LDWORK);
    zgemm('N', 'N', M, L, N - L, Complex.one, B, LDB, V, LDV, Complex.one, WORK,
        LDWORK);
    zgemm('N', 'N', M, K - L, N, Complex.one, B, LDB, V(1, KP), LDV,
        Complex.zero, WORK(1, KP), LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('R', 'U', TRANS, 'N', M, K, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('N', 'C', M, N - L, K, -Complex.one, WORK, LDWORK, V, LDV,
        Complex.one, B, LDB);
    zgemm('N', 'C', M, L, K - L, -Complex.one, WORK(1, KP), LDWORK, V(NP, KP),
        LDV, Complex.one, B(1, NP), LDB);
    ztrmm('R', 'U', 'C', 'N', M, L, Complex.one, V(NP, 1), LDV, WORK, LDWORK);
    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        B[I][N - L + J] -= WORK[I][J];
      }
    }
  } else if (COLUMN && BACKWARD && LEFT) {
    // Let  W =  [ V ]    (M-by-K)
    //           [ I ]    (K-by-K)
    //
    // Form  H C  or  H**H C  where  C = [ B ]  (M-by-N)
    //                                   [ A ]  (K-by-N)
    //
    // H = I - W T W**H          or  H**H = I - W T**H W**H
    //
    // A -=   T (A + V**H B)  or  A -=   T**H (A + V**H B)
    // B -= V T (A + V**H B)  or  B -= V T**H (A + V**H B)

    MP = min(L + 1, M);
    KP = min(K - L + 1, K);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        WORK[K - L + I][J] = B[I][J];
      }
    }

    ztrmm('L', 'L', 'C', 'N', L, N, Complex.one, V(1, KP), LDV, WORK(KP, 1),
        LDWORK);
    zgemm('C', 'N', L, N, M - L, Complex.one, V(MP, KP), LDV, B(MP, 1), LDB,
        Complex.one, WORK(KP, 1), LDWORK);
    zgemm('C', 'N', K - L, N, M, Complex.one, V, LDV, B, LDB, Complex.zero,
        WORK, LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('L', 'L', TRANS, 'N', K, N, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('N', 'N', M - L, N, K, -Complex.one, V(MP, 1), LDV, WORK, LDWORK,
        Complex.one, B(MP, 1), LDB);
    zgemm('N', 'N', L, N, K - L, -Complex.one, V, LDV, WORK, LDWORK,
        Complex.one, B, LDB);
    ztrmm('L', 'L', 'N', 'N', L, N, Complex.one, V(1, KP), LDV, WORK(KP, 1),
        LDWORK);
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        B[I][J] -= WORK[K - L + I][J];
      }
    }
  } else if (COLUMN && BACKWARD && RIGHT) {
    // Let  W =  [ V ]    (N-by-K)
    //           [ I ]    (K-by-K)
    //
    // Form  C H  or  C H**H  where  C = [ B A ] (B is M-by-N, A is M-by-K)
    //
    // H = I - W T W**H          or  H**H = I - W T**H W**H
    //
    // A -= (A + B V) T      or  A -= (A + B V) T**H
    // B -= (A + B V) T V**H  or  B -= (A + B V) T**H V**H

    NP = min(L + 1, N);
    KP = min(K - L + 1, K);

    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][K - L + J] = B[I][J];
      }
    }
    ztrmm('R', 'L', 'N', 'N', M, L, Complex.one, V(1, KP), LDV, WORK(1, KP),
        LDWORK);
    zgemm('N', 'N', M, L, N - L, Complex.one, B(1, NP), LDB, V(NP, KP), LDV,
        Complex.one, WORK(1, KP), LDWORK);
    zgemm('N', 'N', M, K - L, N, Complex.one, B, LDB, V, LDV, Complex.zero,
        WORK, LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('R', 'L', TRANS, 'N', M, K, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('N', 'C', M, N - L, K, -Complex.one, WORK, LDWORK, V(NP, 1), LDV,
        Complex.one, B(1, NP), LDB);
    zgemm('N', 'C', M, L, K - L, -Complex.one, WORK, LDWORK, V, LDV,
        Complex.one, B, LDB);
    ztrmm('R', 'L', 'C', 'N', M, L, Complex.one, V(1, KP), LDV, WORK(1, KP),
        LDWORK);
    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        B[I][J] -= WORK[I][K - L + J];
      }
    }
  } else if (ROW && FORWARD && LEFT) {
    // Let  W =  [ I V ] ( I is K-by-K, V is K-by-M )
    //
    // Form  H C  or  H**H C  where  C = [ A ]  (K-by-N)
    //                                   [ B ]  (M-by-N)
    //
    // H = I - W**H T W          or  H**H = I - W**H T**H W
    //
    // A -=     T (A + V B)  or  A -=     T**H (A + V B)
    // B -= V**H T (A + V B)  or  B -= V**H T**H (A + V B)

    MP = min(M - L + 1, M);
    KP = min(L + 1, K);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        WORK[I][J] = B[M - L + I][J];
      }
    }
    ztrmm('L', 'L', 'N', 'N', L, N, Complex.one, V(1, MP), LDV, WORK, LDB);
    zgemm('N', 'N', L, N, M - L, Complex.one, V, LDV, B, LDB, Complex.one, WORK,
        LDWORK);
    zgemm('N', 'N', K - L, N, M, Complex.one, V(KP, 1), LDV, B, LDB,
        Complex.zero, WORK(KP, 1), LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('L', 'U', TRANS, 'N', K, N, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('C', 'N', M - L, N, K, -Complex.one, V, LDV, WORK, LDWORK,
        Complex.one, B, LDB);
    zgemm('C', 'N', L, N, K - L, -Complex.one, V(KP, MP), LDV, WORK(KP, 1),
        LDWORK, Complex.one, B(MP, 1), LDB);
    ztrmm('L', 'L', 'C', 'N', L, N, Complex.one, V(1, MP), LDV, WORK, LDWORK);
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        B[M - L + I][J] -= WORK[I][J];
      }
    }
  } else if (ROW && FORWARD && RIGHT) {
    // Let  W =  [ I V ] ( I is K-by-K, V is K-by-N )
    //
    // Form  C H  or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N)
    //
    // H = I - W**H T W            or  H**H = I - W**H T**H W
    //
    // A -= (A + B V**H) T      or  A -= (A + B V**H) T**H
    // B -= (A + B V**H) T V    or  B -= (A + B V**H) T**H V

    NP = min(N - L + 1, N);
    KP = min(L + 1, K);

    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][J] = B[I][N - L + J];
      }
    }
    ztrmm('R', 'L', 'C', 'N', M, L, Complex.one, V(1, NP), LDV, WORK, LDWORK);
    zgemm('N', 'C', M, L, N - L, Complex.one, B, LDB, V, LDV, Complex.one, WORK,
        LDWORK);
    zgemm('N', 'C', M, K - L, N, Complex.one, B, LDB, V(KP, 1), LDV,
        Complex.zero, WORK(1, KP), LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('R', 'U', TRANS, 'N', M, K, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('N', 'N', M, N - L, K, -Complex.one, WORK, LDWORK, V, LDV,
        Complex.one, B, LDB);
    zgemm('N', 'N', M, L, K - L, -Complex.one, WORK(1, KP), LDWORK, V(KP, NP),
        LDV, Complex.one, B(1, NP), LDB);
    ztrmm('R', 'L', 'N', 'N', M, L, Complex.one, V(1, NP), LDV, WORK, LDWORK);
    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        B[I][N - L + J] -= WORK[I][J];
      }
    }
  } else if (ROW && BACKWARD && LEFT) {
    // Let  W =  [ V I ] ( I is K-by-K, V is K-by-M )
    //
    // Form  H C  or  H**H C  where  C = [ B ]  (M-by-N)
    //                                   [ A ]  (K-by-N)
    //
    // H = I - W**H T W          or  H**H = I - W**H T**H W
    //
    // A -=     T (A + V B)  or  A -=     T**H (A + V B)
    // B -= V**H T (A + V B)  or  B -= V**H T**H (A + V B)

    MP = min(L + 1, M);
    KP = min(K - L + 1, K);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        WORK[K - L + I][J] = B[I][J];
      }
    }
    ztrmm('L', 'U', 'N', 'N', L, N, Complex.one, V(KP, 1), LDV, WORK(KP, 1),
        LDWORK);
    zgemm('N', 'N', L, N, M - L, Complex.one, V(KP, MP), LDV, B(MP, 1), LDB,
        Complex.one, WORK(KP, 1), LDWORK);
    zgemm('N', 'N', K - L, N, M, Complex.one, V, LDV, B, LDB, Complex.zero,
        WORK, LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('L', 'L ', TRANS, 'N', K, N, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= N; J++) {
      for (I = 1; I <= K; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('C', 'N', M - L, N, K, -Complex.one, V(1, MP), LDV, WORK, LDWORK,
        Complex.one, B(MP, 1), LDB);
    zgemm('C', 'N', L, N, K - L, -Complex.one, V, LDV, WORK, LDWORK,
        Complex.one, B, LDB);
    ztrmm('L', 'U', 'C', 'N', L, N, Complex.one, V(KP, 1), LDV, WORK(KP, 1),
        LDWORK);
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= L; I++) {
        B[I][J] -= WORK[K - L + I][J];
      }
    }
  } else if (ROW && BACKWARD && RIGHT) {
    // Let  W =  [ V I ] ( I is K-by-K, V is K-by-N )
    //
    // Form  C H  or  C H**H  where  C = [ B A ] (A is M-by-K, B is M-by-N)
    //
    // H = I - W**H T W            or  H**H = I - W**H T**H W
    //
    // A -= (A + B V**H) T      or  A -= (A + B V**H) T**H
    // B -= (A + B V**H) T V    or  B -= (A + B V**H) T**H V

    NP = min(L + 1, N);
    KP = min(K - L + 1, K);

    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][K - L + J] = B[I][J];
      }
    }
    ztrmm('R', 'U', 'C', 'N', M, L, Complex.one, V(KP, 1), LDV, WORK(1, KP),
        LDWORK);
    zgemm('N', 'C', M, L, N - L, Complex.one, B(1, NP), LDB, V(KP, NP), LDV,
        Complex.one, WORK(1, KP), LDWORK);
    zgemm('N', 'C', M, K - L, N, Complex.one, B, LDB, V, LDV, Complex.zero,
        WORK, LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        WORK[I][J] += A[I][J];
      }
    }

    ztrmm('R', 'L', TRANS, 'N', M, K, Complex.one, T, LDT, WORK, LDWORK);

    for (J = 1; J <= K; J++) {
      for (I = 1; I <= M; I++) {
        A[I][J] -= WORK[I][J];
      }
    }

    zgemm('N', 'N', M, N - L, K, -Complex.one, WORK, LDWORK, V(1, NP), LDV,
        Complex.one, B(1, NP), LDB);
    zgemm('N', 'N', M, L, K - L, -Complex.one, WORK, LDWORK, V, LDV,
        Complex.one, B, LDB);
    ztrmm('R', 'U', 'N', 'N', M, L, Complex.one, V(KP, 1), LDV, WORK(1, KP),
        LDWORK);
    for (J = 1; J <= L; J++) {
      for (I = 1; I <= M; I++) {
        B[I][J] -= WORK[I][K - L + J];
      }
    }
  }
}
