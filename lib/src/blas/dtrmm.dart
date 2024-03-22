import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/matrix.dart';

void dtrmm(
  final String SIDE,
  final String UPLO,
  final String TRANSA,
  final String DIAG,
  final int M,
  final int N,
  final double ALPHA,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
) {
// -- Reference BLAS level3 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ONE = 1.0, ZERO = 0.0;

  // Test the input parameters.

  final LSIDE = lsame(SIDE, 'L');
  final NROWA = LSIDE ? M : N;
  final NOUNIT = lsame(DIAG, 'N');
  final UPPER = lsame(UPLO, 'U');

  var INFO = 0;
  if (!LSIDE && !lsame(SIDE, 'R')) {
    INFO = 1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO = 2;
  } else if (!lsame(TRANSA, 'N') &&
      !lsame(TRANSA, 'T') &&
      !lsame(TRANSA, 'C')) {
    INFO = 3;
  } else if (!lsame(DIAG, 'U') && !lsame(DIAG, 'N')) {
    INFO = 4;
  } else if (M < 0) {
    INFO = 5;
  } else if (N < 0) {
    INFO = 6;
  } else if (LDA < max(1, NROWA)) {
    INFO = 9;
  } else if (LDB < max(1, M)) {
    INFO = 11;
  }
  if (INFO != 0) {
    xerbla('DTRMM ', INFO);
    return;
  }

  // Quick return if possible.

  if (M == 0 || N == 0) return;

  // And when  alpha == zero.

  if (ALPHA == ZERO) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= M; I++) {
        B[I][J] = ZERO;
      }
    }
    return;
  }

  // Start the operations.

  if (LSIDE) {
    if (lsame(TRANSA, 'N')) {
      // Form  B := alpha*A*B.

      if (UPPER) {
        for (var J = 1; J <= N; J++) {
          for (var K = 1; K <= M; K++) {
            if (B[K][J] != ZERO) {
              var TEMP = ALPHA * B[K][J];
              for (var I = 1; I <= K - 1; I++) {
                B[I][J] += TEMP * A[I][K];
              }
              if (NOUNIT) TEMP *= A[K][K];
              B[K][J] = TEMP;
            }
          }
        }
      } else {
        for (var J = 1; J <= N; J++) {
          for (var K = M; K >= 1; K--) {
            if (B[K][J] != ZERO) {
              final TEMP = ALPHA * B[K][J];
              B[K][J] = TEMP;
              if (NOUNIT) B[K][J] *= A[K][K];
              for (var I = K + 1; I <= M; I++) {
                B[I][J] += TEMP * A[I][K];
              }
            }
          }
        }
      }
    } else {
      // Form  B := alpha*A**T*B.

      if (UPPER) {
        for (var J = 1; J <= N; J++) {
          for (var I = M; I >= 1; I--) {
            var TEMP = B[I][J];
            if (NOUNIT) TEMP *= A[I][I];
            for (var K = 1; K <= I - 1; K++) {
              TEMP += A[K][I] * B[K][J];
            }
            B[I][J] = ALPHA * TEMP;
          }
        }
      } else {
        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= M; I++) {
            var TEMP = B[I][J];
            if (NOUNIT) TEMP *= A[I][I];
            for (var K = I + 1; K <= M; K++) {
              TEMP += A[K][I] * B[K][J];
            }
            B[I][J] = ALPHA * TEMP;
          }
        }
      }
    }
  } else {
    if (lsame(TRANSA, 'N')) {
      // Form  B := alpha*B*A.

      if (UPPER) {
        for (var J = N; J >= 1; J--) {
          var TEMP = ALPHA;
          if (NOUNIT) TEMP *= A[J][J];
          for (var I = 1; I <= M; I++) {
            B[I][J] *= TEMP;
          }
          for (var K = 1; K <= J - 1; K++) {
            if (A[K][J] != ZERO) {
              TEMP = ALPHA * A[K][J];
              for (var I = 1; I <= M; I++) {
                B[I][J] += TEMP * B[I][K];
              }
            }
          }
        }
      } else {
        for (var J = 1; J <= N; J++) {
          var TEMP = ALPHA;
          if (NOUNIT) TEMP *= A[J][J];
          for (var I = 1; I <= M; I++) {
            B[I][J] *= TEMP;
          }
          for (var K = J + 1; K <= N; K++) {
            if (A[K][J] != ZERO) {
              TEMP = ALPHA * A[K][J];
              for (var I = 1; I <= M; I++) {
                B[I][J] += TEMP * B[I][K];
              }
            }
          }
        }
      }
    } else {
      // Form  B := alpha*B*A**T.

      if (UPPER) {
        for (var K = 1; K <= N; K++) {
          for (var J = 1; J <= K - 1; J++) {
            if (A[J][K] != ZERO) {
              final TEMP = ALPHA * A[J][K];
              for (var I = 1; I <= M; I++) {
                B[I][J] += TEMP * B[I][K];
              }
            }
          }
          var TEMP = ALPHA;
          if (NOUNIT) TEMP *= A[K][K];
          if (TEMP != ONE) {
            for (var I = 1; I <= M; I++) {
              B[I][K] *= TEMP;
            }
          }
        }
      } else {
        for (var K = N; K >= 1; K--) {
          for (var J = K + 1; J <= N; J++) {
            if (A[J][K] != ZERO) {
              final TEMP = ALPHA * A[J][K];
              for (var I = 1; I <= M; I++) {
                B[I][J] += TEMP * B[I][K];
              }
            }
          }
          var TEMP = ALPHA;
          if (NOUNIT) TEMP *= A[K][K];
          if (TEMP != ONE) {
            for (var I = 1; I <= M; I++) {
              B[I][K] *= TEMP;
            }
          }
        }
      }
    }
  }
}
