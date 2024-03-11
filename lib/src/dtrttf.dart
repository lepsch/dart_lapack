import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtrttf(
  final String TRANSR,
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> ARF_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA, offset: zeroIndexedMatrixOffset);
  final ARF = ARF_.having(offset: zeroIndexedArrayOffset);
  bool LOWER, NISODD, NORMALTRANSR;
  int I, IJ, J, K = 0, L, N1, N2, NT, NX2 = 0, NP1X2 = 0;

  // Test the input parameters.

  INFO.value = 0;
  NORMALTRANSR = lsame(TRANSR, 'N');
  LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'T')) {
    INFO.value = -1;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  }
  if (INFO.value != 0) {
    xerbla('DTRTTF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 1) {
    if (N == 1) {
      ARF[0] = A[0][0];
    }
    return;
  }

  // Size of array ARF(0:nt-1)

  NT = N * (N + 1) ~/ 2;

  // Set N1 and N2 depending on LOWER: for N even N1=N2=K

  if (LOWER) {
    N2 = N ~/ 2;
    N1 = N - N2;
  } else {
    N1 = N ~/ 2;
    N2 = N - N1;
  }

  // If N is odd, set NISODD = true , LDA=N+1 and A is (N+1)--by--K2.
  // If N is even, set K = N/2 and NISODD = false , LDA=N and A is
  // N--by--(N+1)/2.

  if ((N % 2) == 0) {
    K = N ~/ 2;
    NISODD = false;
    if (!LOWER) NP1X2 = N + N + 2;
  } else {
    NISODD = true;
    if (!LOWER) NX2 = N + N;
  }

  if (NISODD) {
    // N is odd

    if (NORMALTRANSR) {
      // N is odd and TRANSR = 'N'

      if (LOWER) {
        // N is odd, TRANSR = 'N', and UPLO = 'L'

        IJ = 0;
        for (J = 0; J <= N2; J++) {
          for (I = N1; I <= N2 + J; I++) {
            ARF[IJ] = A[N2 + J][I];
            IJ++;
          }
          for (I = J; I <= N - 1; I++) {
            ARF[IJ] = A[I][J];
            IJ++;
          }
        }
      } else {
        // N is odd, TRANSR = 'N', and UPLO = 'U'

        IJ = NT - N;
        for (J = N - 1; J >= N1; J--) {
          for (I = 0; I <= J; I++) {
            ARF[IJ] = A[I][J];
            IJ++;
          }
          for (L = J - N1; L <= N1 - 1; L++) {
            ARF[IJ] = A[J - N1][L];
            IJ++;
          }
          IJ = IJ - NX2;
        }
      }
    } else {
      // N is odd and TRANSR = 'T'

      if (LOWER) {
        // N is odd, TRANSR = 'T', and UPLO = 'L'

        IJ = 0;
        for (J = 0; J <= N2 - 1; J++) {
          for (I = 0; I <= J; I++) {
            ARF[IJ] = A[J][I];
            IJ++;
          }
          for (I = N1 + J; I <= N - 1; I++) {
            ARF[IJ] = A[I][N1 + J];
            IJ++;
          }
        }
        for (J = N2; J <= N - 1; J++) {
          for (I = 0; I <= N1 - 1; I++) {
            ARF[IJ] = A[J][I];
            IJ++;
          }
        }
      } else {
        // N is odd, TRANSR = 'T', and UPLO = 'U'

        IJ = 0;
        for (J = 0; J <= N1; J++) {
          for (I = N1; I <= N - 1; I++) {
            ARF[IJ] = A[J][I];
            IJ++;
          }
        }
        for (J = 0; J <= N1 - 1; J++) {
          for (I = 0; I <= J; I++) {
            ARF[IJ] = A[I][J];
            IJ++;
          }
          for (L = N2 + J; L <= N - 1; L++) {
            ARF[IJ] = A[N2 + J][L];
            IJ++;
          }
        }
      }
    }
  } else {
    // N is even

    if (NORMALTRANSR) {
      // N is even and TRANSR = 'N'

      if (LOWER) {
        // N is even, TRANSR = 'N', and UPLO = 'L'

        IJ = 0;
        for (J = 0; J <= K - 1; J++) {
          for (I = K; I <= K + J; I++) {
            ARF[IJ] = A[K + J][I];
            IJ++;
          }
          for (I = J; I <= N - 1; I++) {
            ARF[IJ] = A[I][J];
            IJ++;
          }
        }
      } else {
        // N is even, TRANSR = 'N', and UPLO = 'U'

        IJ = NT - N - 1;
        for (J = N - 1; J >= K; J--) {
          for (I = 0; I <= J; I++) {
            ARF[IJ] = A[I][J];
            IJ++;
          }
          for (L = J - K; L <= K - 1; L++) {
            ARF[IJ] = A[J - K][L];
            IJ++;
          }
          IJ = IJ - NP1X2;
        }
      }
    } else {
      // N is even and TRANSR = 'T'

      if (LOWER) {
        // N is even, TRANSR = 'T', and UPLO = 'L'

        IJ = 0;
        J = K;
        for (I = K; I <= N - 1; I++) {
          ARF[IJ] = A[I][J];
          IJ++;
        }
        for (J = 0; J <= K - 2; J++) {
          for (I = 0; I <= J; I++) {
            ARF[IJ] = A[J][I];
            IJ++;
          }
          for (I = K + 1 + J; I <= N - 1; I++) {
            ARF[IJ] = A[I][K + 1 + J];
            IJ++;
          }
        }
        for (J = K - 1; J <= N - 1; J++) {
          for (I = 0; I <= K - 1; I++) {
            ARF[IJ] = A[J][I];
            IJ++;
          }
        }
      } else {
        // N is even, TRANSR = 'T', and UPLO = 'U'

        IJ = 0;
        for (J = 0; J <= K; J++) {
          for (I = K; I <= N - 1; I++) {
            ARF[IJ] = A[J][I];
            IJ++;
          }
        }
        for (J = 0; J <= K - 2; J++) {
          for (I = 0; I <= J; I++) {
            ARF[IJ] = A[I][J];
            IJ++;
          }
          for (L = K + 1 + J; L <= N - 1; L++) {
            ARF[IJ] = A[K + 1 + J][L];
            IJ++;
          }
        }
        // Note that here, on exit of the loop, J = K-1
        for (I = 0; I <= J; I++) {
          ARF[IJ] = A[I][J];
          IJ++;
        }
      }
    }
  }
}
