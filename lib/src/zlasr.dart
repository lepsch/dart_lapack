import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zlasr(
  final String SIDE,
  final String PIVOT,
  final String DIRECT,
  final int M,
  final int N,
  final Array<double> C_,
  final Array<double> S_,
  final Matrix<Complex> A_,
  final int LDA,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final C = C_.having();
  final S = S_.having();
  final A = A_.having(ld: LDA);
  const ONE = 1.0, ZERO = 0.0;
  int I, INFO, J;
  double CTEMP, STEMP;
  Complex TEMP;

  // Test the input parameters

  INFO = 0;
  if (!(lsame(SIDE, 'L') || lsame(SIDE, 'R'))) {
    INFO = 1;
  } else if (!(lsame(PIVOT, 'V') || lsame(PIVOT, 'T') || lsame(PIVOT, 'B'))) {
    INFO = 2;
  } else if (!(lsame(DIRECT, 'F') || lsame(DIRECT, 'B'))) {
    INFO = 3;
  } else if (M < 0) {
    INFO = 4;
  } else if (N < 0) {
    INFO = 5;
  } else if (LDA < max(1, M)) {
    INFO = 9;
  }
  if (INFO != 0) {
    xerbla('ZLASR', INFO);
    return;
  }

  // Quick return if possible

  if ((M == 0) || (N == 0)) return;
  if (lsame(SIDE, 'L')) {
    // Form  P * A

    if (lsame(PIVOT, 'V')) {
      if (lsame(DIRECT, 'F')) {
        for (J = 1; J <= M - 1; J++) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= N; I++) {
              TEMP = A[J + 1][I];
              A[J + 1][I] =
                  CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[J][I];
              A[J][I] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[J][I];
            }
          }
        }
      } else if (lsame(DIRECT, 'B')) {
        for (J = M - 1; J >= 1; J--) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= N; I++) {
              TEMP = A[J + 1][I];
              A[J + 1][I] =
                  CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[J][I];
              A[J][I] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[J][I];
            }
          }
        }
      }
    } else if (lsame(PIVOT, 'T')) {
      if (lsame(DIRECT, 'F')) {
        for (J = 2; J <= M; J++) {
          CTEMP = C[J - 1];
          STEMP = S[J - 1];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= N; I++) {
              TEMP = A[J][I];
              A[J][I] = CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[1][I];
              A[1][I] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[1][I];
            }
          }
        }
      } else if (lsame(DIRECT, 'B')) {
        for (J = M; J >= 2; J--) {
          CTEMP = C[J - 1];
          STEMP = S[J - 1];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= N; I++) {
              TEMP = A[J][I];
              A[J][I] = CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[1][I];
              A[1][I] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[1][I];
            }
          }
        }
      }
    } else if (lsame(PIVOT, 'B')) {
      if (lsame(DIRECT, 'F')) {
        for (J = 1; J <= M - 1; J++) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= N; I++) {
              TEMP = A[J][I];
              A[J][I] = STEMP.toComplex() * A[M][I] + CTEMP.toComplex() * TEMP;
              A[M][I] = CTEMP.toComplex() * A[M][I] - STEMP.toComplex() * TEMP;
            }
          }
        }
      } else if (lsame(DIRECT, 'B')) {
        for (J = M - 1; J >= 1; J--) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= N; I++) {
              TEMP = A[J][I];
              A[J][I] = STEMP.toComplex() * A[M][I] + CTEMP.toComplex() * TEMP;
              A[M][I] = CTEMP.toComplex() * A[M][I] - STEMP.toComplex() * TEMP;
            }
          }
        }
      }
    }
  } else if (lsame(SIDE, 'R')) {
    // Form A * P**T

    if (lsame(PIVOT, 'V')) {
      if (lsame(DIRECT, 'F')) {
        for (J = 1; J <= N - 1; J++) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= M; I++) {
              TEMP = A[I][J + 1];
              A[I][J + 1] =
                  CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[I][J];
              A[I][J] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[I][J];
            }
          }
        }
      } else if (lsame(DIRECT, 'B')) {
        for (J = N - 1; J >= 1; J--) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= M; I++) {
              TEMP = A[I][J + 1];
              A[I][J + 1] =
                  CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[I][J];
              A[I][J] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[I][J];
            }
          }
        }
      }
    } else if (lsame(PIVOT, 'T')) {
      if (lsame(DIRECT, 'F')) {
        for (J = 2; J <= N; J++) {
          CTEMP = C[J - 1];
          STEMP = S[J - 1];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= M; I++) {
              TEMP = A[I][J];
              A[I][J] = CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[I][1];
              A[I][1] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[I][1];
            }
          }
        }
      } else if (lsame(DIRECT, 'B')) {
        for (J = N; J >= 2; J--) {
          CTEMP = C[J - 1];
          STEMP = S[J - 1];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= M; I++) {
              TEMP = A[I][J];
              A[I][J] = CTEMP.toComplex() * TEMP - STEMP.toComplex() * A[I][1];
              A[I][1] = STEMP.toComplex() * TEMP + CTEMP.toComplex() * A[I][1];
            }
          }
        }
      }
    } else if (lsame(PIVOT, 'B')) {
      if (lsame(DIRECT, 'F')) {
        for (J = 1; J <= N - 1; J++) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= M; I++) {
              TEMP = A[I][J];
              A[I][J] = STEMP.toComplex() * A[I][N] + CTEMP.toComplex() * TEMP;
              A[I][N] = CTEMP.toComplex() * A[I][N] - STEMP.toComplex() * TEMP;
            }
          }
        }
      } else if (lsame(DIRECT, 'B')) {
        for (J = N - 1; J >= 1; J--) {
          CTEMP = C[J];
          STEMP = S[J];
          if ((CTEMP != ONE) || (STEMP != ZERO)) {
            for (I = 1; I <= M; I++) {
              TEMP = A[I][J];
              A[I][J] = STEMP.toComplex() * A[I][N] + CTEMP.toComplex() * TEMP;
              A[I][N] = CTEMP.toComplex() * A[I][N] - STEMP.toComplex() * TEMP;
            }
          }
        }
      }
    }
  }
}
