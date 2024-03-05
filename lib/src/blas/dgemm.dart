import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dgemm(
  final String TRANSA,
  final String TRANSB,
  final int M,
  final int N,
  final int K,
  final double ALPHA,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final double BETA,
  final Matrix<double> C_,
  final int LDC,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  double TEMP;
  int I, INFO, J, L, NROWA, NROWB;
  bool NOTA, NOTB;
  const ONE = 1.0, ZERO = 0.0;

  // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  // transposed and set  NROWA and NROWB  as the number of rows of  A
  // and  B  respectively.

  NOTA = lsame(TRANSA, 'N');
  NOTB = lsame(TRANSB, 'N');
  if (NOTA) {
    NROWA = M;
  } else {
    NROWA = K;
  }
  if (NOTB) {
    NROWB = K;
  } else {
    NROWB = N;
  }

  // Test the input parameters.

  INFO = 0;
  if ((!NOTA) && (!lsame(TRANSA, 'C')) && (!lsame(TRANSA, 'T'))) {
    INFO = 1;
  } else if ((!NOTB) && (!lsame(TRANSB, 'C')) && (!lsame(TRANSB, 'T'))) {
    INFO = 2;
  } else if (M < 0) {
    INFO = 3;
  } else if (N < 0) {
    INFO = 4;
  } else if (K < 0) {
    INFO = 5;
  } else if (LDA < max(1, NROWA)) {
    INFO = 8;
  } else if (LDB < max(1, NROWB)) {
    INFO = 10;
  } else if (LDC < max(1, M)) {
    INFO = 13;
  }
  if (INFO != 0) {
    xerbla('DGEMM ', INFO);
    return;
  }

  // Quick return if possible.

  if ((M == 0) ||
      (N == 0) ||
      (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) {
    return;
  }

  // And if  alpha == zero.

  if (ALPHA == ZERO) {
    if (BETA == ZERO) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          C[I][J] = ZERO;
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          C[I][J] = BETA * C[I][J];
        }
      }
    }
    return;
  }

  // Start the operations.

  if (NOTB) {
    if (NOTA) {
      // Form  C := alpha*A*B + beta*C.

      for (J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (I = 1; I <= M; I++) {
            C[I][J] = ZERO;
          }
        } else if (BETA != ONE) {
          for (I = 1; I <= M; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
        for (L = 1; L <= K; L++) {
          TEMP = ALPHA * B[L][J];
          for (I = 1; I <= M; I++) {
            C[I][J] = C[I][J] + TEMP * A[I][L];
          }
        }
      }
    } else {
      // Form  C := alpha*A**T*B + beta*C

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = ZERO;
          for (L = 1; L <= K; L++) {
            TEMP = TEMP + A[L][I] * B[L][J];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    }
  } else {
    if (NOTA) {
      // Form  C := alpha*A*B**T + beta*C

      for (J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (I = 1; I <= M; I++) {
            C[I][J] = ZERO;
          }
        } else if (BETA != ONE) {
          for (I = 1; I <= M; I++) {
            C[I][J] = BETA * C[I][J];
          }
        }
        for (L = 1; L <= K; L++) {
          TEMP = ALPHA * B[J][L];
          for (I = 1; I <= M; I++) {
            C[I][J] = C[I][J] + TEMP * A[I][L];
          }
        }
      }
    } else {
      // Form  C := alpha*A**T*B**T + beta*C

      for (J = 1; J <= N; J++) {
        for (I = 1; I <= M; I++) {
          TEMP = ZERO;
          for (L = 1; L <= K; L++) {
            TEMP = TEMP + A[L][I] * B[J][L];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA * TEMP;
          } else {
            C[I][J] = ALPHA * TEMP + BETA * C[I][J];
          }
        }
      }
    }
  }
}
