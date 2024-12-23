// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dtfttr(
  final String TRANSR,
  final String UPLO,
  final int N,
  final Array<double> ARF_,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA, offset: zeroIndexedMatrixOffset);
  final ARF = ARF_.having(offset: zeroIndexedArrayOffset);
  bool NISODD;
  int N1, N2, K = 0, NT, NX2 = 0, NP1X2 = 0;
  int I, J, L, IJ;

  // Test the input parameters.

  INFO.value = 0;
  final NORMALTRANSR = lsame(TRANSR, 'N');
  final LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'T')) {
    INFO.value = -1;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DTFTTR', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 1) {
    if (N == 1) {
      A[0][0] = ARF[0];
    }
    return;
  }

  // Size of array ARF(0:nt-1)

  NT = N * (N + 1) ~/ 2;

  // set N1 and N2 depending on LOWER: for N even N1=N2=K

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
            A[N2 + J][I] = ARF[IJ];
            IJ++;
          }
          for (I = J; I <= N - 1; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
        }
      } else {
        // N is odd, TRANSR = 'N', and UPLO = 'U'

        IJ = NT - N;
        for (J = N - 1; J >= N1; J--) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = J - N1; L <= N1 - 1; L++) {
            A[J - N1][L] = ARF[IJ];
            IJ++;
          }
          IJ -= NX2;
        }
      }
    } else {
      // N is odd and TRANSR = 'T'

      if (LOWER) {
        // N is odd, TRANSR = 'T', and UPLO = 'L'

        IJ = 0;
        for (J = 0; J <= N2 - 1; J++) {
          for (I = 0; I <= J; I++) {
            A[J][I] = ARF[IJ];
            IJ++;
          }
          for (I = N1 + J; I <= N - 1; I++) {
            A[I][N1 + J] = ARF[IJ];
            IJ++;
          }
        }
        for (J = N2; J <= N - 1; J++) {
          for (I = 0; I <= N1 - 1; I++) {
            A[J][I] = ARF[IJ];
            IJ++;
          }
        }
      } else {
        // N is odd, TRANSR = 'T', and UPLO = 'U'

        IJ = 0;
        for (J = 0; J <= N1; J++) {
          for (I = N1; I <= N - 1; I++) {
            A[J][I] = ARF[IJ];
            IJ++;
          }
        }
        for (J = 0; J <= N1 - 1; J++) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = N2 + J; L <= N - 1; L++) {
            A[N2 + J][L] = ARF[IJ];
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
            A[K + J][I] = ARF[IJ];
            IJ++;
          }
          for (I = J; I <= N - 1; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
        }
      } else {
        // N is even, TRANSR = 'N', and UPLO = 'U'

        IJ = NT - N - 1;
        for (J = N - 1; J >= K; J--) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = J - K; L <= K - 1; L++) {
            A[J - K][L] = ARF[IJ];
            IJ++;
          }
          IJ -= NP1X2;
        }
      }
    } else {
      // N is even and TRANSR = 'T'

      if (LOWER) {
        // N is even, TRANSR = 'T', and UPLO = 'L'

        IJ = 0;
        J = K;
        for (I = K; I <= N - 1; I++) {
          A[I][J] = ARF[IJ];
          IJ++;
        }
        for (J = 0; J <= K - 2; J++) {
          for (I = 0; I <= J; I++) {
            A[J][I] = ARF[IJ];
            IJ++;
          }
          for (I = K + 1 + J; I <= N - 1; I++) {
            A[I][K + 1 + J] = ARF[IJ];
            IJ++;
          }
        }
        for (J = K - 1; J <= N - 1; J++) {
          for (I = 0; I <= K - 1; I++) {
            A[J][I] = ARF[IJ];
            IJ++;
          }
        }
      } else {
        // N is even, TRANSR = 'T', and UPLO = 'U'

        IJ = 0;
        for (J = 0; J <= K; J++) {
          for (I = K; I <= N - 1; I++) {
            A[J][I] = ARF[IJ];
            IJ++;
          }
        }
        for (J = 0; J <= K - 2; J++) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = K + 1 + J; L <= N - 1; L++) {
            A[K + 1 + J][L] = ARF[IJ];
            IJ++;
          }
        }
        // Note that here, on exit of the loop, J = K-1
        for (I = 0; I <= J; I++) {
          A[I][J] = ARF[IJ];
          IJ++;
        }
      }
    }
  }
}
