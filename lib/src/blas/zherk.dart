import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void zherk(
  final String UPLO,
  final String TRANS,
  final int N,
  final int K,
  final double ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final double BETA,
  final Matrix<Complex> C_,
  final int LDC,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final C = C_.dim(LDC);
  Complex TEMP;
  double RTEMP;
  int I, INFO, J, L, NROWA;
  bool UPPER;
  const ONE = 1.0, ZERO = 0.0;

  // Test the input parameters.

  if (lsame(TRANS, 'N')) {
    NROWA = N;
  } else {
    NROWA = K;
  }
  UPPER = lsame(UPLO, 'U');

  INFO = 0;
  if ((!UPPER) && (!lsame(UPLO, 'L'))) {
    INFO = 1;
  } else if ((!lsame(TRANS, 'N')) && (!lsame(TRANS, 'C'))) {
    INFO = 2;
  } else if (N < 0) {
    INFO = 3;
  } else if (K < 0) {
    INFO = 4;
  } else if (LDA < max(1, NROWA)) {
    INFO = 7;
  } else if (LDC < max(1, N)) {
    INFO = 10;
  }
  if (INFO != 0) {
    xerbla('ZHERK ', INFO);
    return;
  }

  // Quick return if possible.

  if ((N == 0) || (((ALPHA == ZERO) || (K == 0)) && (BETA == ONE))) return;

  // And when  alpha == zero.

  if (ALPHA == ZERO) {
    if (UPPER) {
      if (BETA == ZERO) {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J; I++) {
            C[I][J] = Complex.zero;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= J - 1; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
          C[J][J] = BETA.toComplex() * C[J][J].real.toComplex();
        }
      }
    } else {
      if (BETA == ZERO) {
        for (J = 1; J <= N; J++) {
          for (I = J; I <= N; I++) {
            C[I][J] = Complex.zero;
          }
        }
      } else {
        for (J = 1; J <= N; J++) {
          C[J][J] = (BETA * (C[J][J]).toDouble()).toComplex();
          for (I = J + 1; I <= N; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
        }
      }
    }
    return;
  }

  // Start the operations.

  if (lsame(TRANS, 'N')) {
    // Form  C := alpha*A*A**H + beta*C.

    if (UPPER) {
      for (J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (I = 1; I <= J; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != ONE) {
          for (I = 1; I <= J - 1; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
          C[J][J] = (BETA * (C[J][J]).toDouble()).toComplex();
        } else {
          C[J][J] = C[J][J].real.toComplex();
        }
        for (L = 1; L <= K; L++) {
          if (A[J][L] != Complex.zero) {
            TEMP = ALPHA.toComplex() * A[J][L].conjugate();
            for (I = 1; I <= J - 1; I++) {
              C[I][J] = C[I][J] + TEMP * A[I][L];
            }
            C[J][J] =
                (C[J][J].toDouble() + (TEMP * A[I][L]).toDouble()).toComplex();
          }
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        if (BETA == ZERO) {
          for (I = J; I <= N; I++) {
            C[I][J] = Complex.zero;
          }
        } else if (BETA != ONE) {
          C[J][J] = (BETA * (C[J][J]).toDouble()).toComplex();
          for (I = J + 1; I <= N; I++) {
            C[I][J] = BETA.toComplex() * C[I][J];
          }
        } else {
          C[J][J] = C[J][J].real.toComplex();
        }
        for (L = 1; L <= K; L++) {
          if (A[J][L] != Complex.zero) {
            TEMP = ALPHA.toComplex() * A[J][L].conjugate();
            C[J][J] =
                (C[J][J].toDouble() + (TEMP * A[J][L]).toDouble()).toComplex();
            for (I = J + 1; I <= N; I++) {
              C[I][J] = C[I][J] + TEMP * A[I][L];
            }
          }
        }
      }
    }
  } else {
    // Form  C := alpha*A**H*A + beta*C.

    if (UPPER) {
      for (J = 1; J <= N; J++) {
        for (I = 1; I <= J - 1; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP = TEMP + A[L][I].conjugate() * A[L][J];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA.toComplex() * TEMP;
          } else {
            C[I][J] = ALPHA.toComplex() * TEMP + BETA.toComplex() * C[I][J];
          }
        }
        RTEMP = ZERO;
        for (L = 1; L <= K; L++) {
          RTEMP = RTEMP + (A[L][J].conjugate() * A[L][J]).real;
        }
        if (BETA == ZERO) {
          C[J][J] = (ALPHA * RTEMP).toComplex();
        } else {
          C[J][J] = (ALPHA * RTEMP + BETA * C[J][J].toDouble()).toComplex();
        }
      }
    } else {
      for (J = 1; J <= N; J++) {
        RTEMP = ZERO;
        for (L = 1; L <= K; L++) {
          RTEMP = RTEMP + (A[L][J].conjugate() * A[L][J]).toDouble();
        }
        if (BETA == ZERO) {
          C[J][J] = (ALPHA * RTEMP).toComplex();
        } else {
          C[J][J] = (ALPHA * RTEMP + BETA * (C[J][J]).toDouble()).toComplex();
        }
        for (I = J + 1; I <= N; I++) {
          TEMP = Complex.zero;
          for (L = 1; L <= K; L++) {
            TEMP = TEMP + A[L][I].conjugate() * A[L][J];
          }
          if (BETA == ZERO) {
            C[I][J] = ALPHA.toComplex() * TEMP;
          } else {
            C[I][J] = ALPHA.toComplex() * TEMP + BETA.toComplex() * C[I][J];
          }
        }
      }
    }
  }
}
