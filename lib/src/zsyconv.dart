import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void zsyconv(
  final String UPLO,
  final String WAY,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<int> IPIV_,
  final Array<Complex> E_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  bool UPPER, CONVERT;
  int I, IP, J;
  Complex TEMP;

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  CONVERT = lsame(WAY, 'C');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!CONVERT && !lsame(WAY, 'R')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('ZSYCONV', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // A is UPPER

    if (CONVERT) {
      // Convert A (A is upper)

      // Convert VALUE

      I = N;
      E[1] = Complex.zero;
      while (I > 1) {
        if (IPIV[I] < 0) {
          E[I] = A[I - 1][I];
          E[I - 1] = Complex.zero;
          A[I - 1][I] = Complex.zero;
          I--;
        } else {
          E[I] = Complex.zero;
        }
        I--;
      }

      // Convert PERMUTATIONS

      I = N;
      while (I >= 1) {
        if (IPIV[I] > 0) {
          IP = IPIV[I];
          if (I < N) {
            for (J = I + 1; J <= N; J++) {
              // 12
              TEMP = A[IP][J];
              A[IP][J] = A[I][J];
              A[I][J] = TEMP;
            } // 12
          }
        } else {
          IP = -IPIV[I];
          if (I < N) {
            for (J = I + 1; J <= N; J++) {
              // 13
              TEMP = A[IP][J];
              A[IP][J] = A[I - 1][J];
              A[I - 1][J] = TEMP;
            } // 13
          }
          I--;
        }
        I--;
      }
    } else {
      // Revert A (A is upper)

      // Revert PERMUTATIONS

      I = 1;
      while (I <= N) {
        if (IPIV[I] > 0) {
          IP = IPIV[I];
          if (I < N) {
            for (J = I + 1; J <= N; J++) {
              TEMP = A[IP][J];
              A[IP][J] = A[I][J];
              A[I][J] = TEMP;
            }
          }
        } else {
          IP = -IPIV[I];
          I++;
          if (I < N) {
            for (J = I + 1; J <= N; J++) {
              TEMP = A[IP][J];
              A[IP][J] = A[I - 1][J];
              A[I - 1][J] = TEMP;
            }
          }
        }
        I++;
      }

      // Revert VALUE

      I = N;
      while (I > 1) {
        if (IPIV[I] < 0) {
          A[I - 1][I] = E[I];
          I--;
        }
        I--;
      }
    }
  } else {
    // A is LOWER

    if (CONVERT) {
      // Convert A (A is lower)

      // Convert VALUE

      I = 1;
      E[N] = Complex.zero;
      while (I <= N) {
        if (I < N && IPIV[I] < 0) {
          E[I] = A[I + 1][I];
          E[I + 1] = Complex.zero;
          A[I + 1][I] = Complex.zero;
          I++;
        } else {
          E[I] = Complex.zero;
        }
        I++;
      }

      // Convert PERMUTATIONS

      I = 1;
      while (I <= N) {
        if (IPIV[I] > 0) {
          IP = IPIV[I];
          if (I > 1) {
            for (J = 1; J <= I - 1; J++) {
              // 22
              TEMP = A[IP][J];
              A[IP][J] = A[I][J];
              A[I][J] = TEMP;
            } // 22
          }
        } else {
          IP = -IPIV[I];
          if (I > 1) {
            for (J = 1; J <= I - 1; J++) {
              // 23
              TEMP = A[IP][J];
              A[IP][J] = A[I + 1][J];
              A[I + 1][J] = TEMP;
            } // 23
          }
          I++;
        }
        I++;
      }
    } else {
      // Revert A (A is lower)

      // Revert PERMUTATIONS

      I = N;
      while (I >= 1) {
        if (IPIV[I] > 0) {
          IP = IPIV[I];
          if (I > 1) {
            for (J = 1; J <= I - 1; J++) {
              TEMP = A[I][J];
              A[I][J] = A[IP][J];
              A[IP][J] = TEMP;
            }
          }
        } else {
          IP = -IPIV[I];
          I--;
          if (I > 1) {
            for (J = 1; J <= I - 1; J++) {
              TEMP = A[I + 1][J];
              A[I + 1][J] = A[IP][J];
              A[IP][J] = TEMP;
            }
          }
        }
        I--;
      }

      // Revert VALUE

      I = 1;
      while (I <= N - 1) {
        if (IPIV[I] < 0) {
          A[I + 1][I] = E[I];
          I++;
        }
        I++;
      }
    }
  }
}
