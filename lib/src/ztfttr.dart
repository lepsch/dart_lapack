// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void ztfttr(
  final String TRANSR,
  final String UPLO,
  final int N,
  final Array<Complex> ARF_,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> INFO,
) {
  final ARF = ARF_.having(offset: zeroIndexedArrayOffset);
  final A = A_.having(ld: LDA, offset: zeroIndexedMatrixOffset);
  int N1, N2, K = 0, NT, NX2 = 0, NP1X2 = 0;
  int I, J, L, IJ;

  // Test the input parameters.

  INFO.value = 0;
  final NORMALTRANSR = lsame(TRANSR, 'N');
  final LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'C')) {
    INFO.value = -1;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('ZTFTTR', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 1) {
    if (N == 1) {
      if (NORMALTRANSR) {
        A[0][0] = ARF[0];
      } else {
        A[0][0] = ARF[0].conjugate();
      }
    }
    return;
  }

  // Size of array ARF(1:2,0:nt-1)

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
  bool NISODD;
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
        // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
        // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
        // T1 -> a(0), T2 -> a(n), S -> a(n1); lda=n

        IJ = 0;
        for (J = 0; J <= N2; J++) {
          for (I = N1; I <= N2 + J; I++) {
            A[N2 + J][I] = ARF[IJ].conjugate();
            IJ++;
          }
          for (I = J; I <= N - 1; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
        }
      } else {
        // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
        // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
        // T1 -> a(n2), T2 -> a(n1), S -> a(0); lda=n

        IJ = NT - N;
        for (J = N - 1; J >= N1; J--) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = J - N1; L <= N1 - 1; L++) {
            A[J - N1][L] = ARF[IJ].conjugate();
            IJ++;
          }
          IJ -= NX2;
        }
      }
    } else {
      // N is odd and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is odd
        // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
        // T1 -> A(0+0) , T2 -> A(1+0) , S -> A(0+n1*n1); lda=n1

        IJ = 0;
        for (J = 0; J <= N2 - 1; J++) {
          for (I = 0; I <= J; I++) {
            A[J][I] = ARF[IJ].conjugate();
            IJ++;
          }
          for (I = N1 + J; I <= N - 1; I++) {
            A[I][N1 + J] = ARF[IJ];
            IJ++;
          }
        }
        for (J = N2; J <= N - 1; J++) {
          for (I = 0; I <= N1 - 1; I++) {
            A[J][I] = ARF[IJ].conjugate();
            IJ++;
          }
        }
      } else {
        // SRPA for UPPER, TRANSPOSE and N is odd
        // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
        // T1 -> A(n2*n2), T2 -> A(n1*n2), S -> A(0); lda = n2

        IJ = 0;
        for (J = 0; J <= N1; J++) {
          for (I = N1; I <= N - 1; I++) {
            A[J][I] = ARF[IJ].conjugate();
            IJ++;
          }
        }
        for (J = 0; J <= N1 - 1; J++) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = N2 + J; L <= N - 1; L++) {
            A[N2 + J][L] = ARF[IJ].conjugate();
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
        // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
        // T1 -> a(1), T2 -> a(0), S -> a(k+1); lda=n+1

        IJ = 0;
        for (J = 0; J <= K - 1; J++) {
          for (I = K; I <= K + J; I++) {
            A[K + J][I] = ARF[IJ].conjugate();
            IJ++;
          }
          for (I = J; I <= N - 1; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
        }
      } else {
        // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
        // T1 -> a(k+1), T2 -> a(k), S -> a(0); lda=n+1

        IJ = NT - N - 1;
        for (J = N - 1; J >= K; J--) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = J - K; L <= K - 1; L++) {
            A[J - K][L] = ARF[IJ].conjugate();
            IJ++;
          }
          IJ -= NP1X2;
        }
      }
    } else {
      // N is even and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is even (see paper, A=B)
        // T1 -> A(0,1) , T2 -> A(0,0) , S -> A(0,k+1) :
        // T1 -> A(0+k) , T2 -> A(0+0) , S -> A(0+k*(k+1)); lda=k

        IJ = 0;
        J = K;
        for (I = K; I <= N - 1; I++) {
          A[I][J] = ARF[IJ];
          IJ++;
        }
        for (J = 0; J <= K - 2; J++) {
          for (I = 0; I <= J; I++) {
            A[J][I] = ARF[IJ].conjugate();
            IJ++;
          }
          for (I = K + 1 + J; I <= N - 1; I++) {
            A[I][K + 1 + J] = ARF[IJ];
            IJ++;
          }
        }
        for (J = K - 1; J <= N - 1; J++) {
          for (I = 0; I <= K - 1; I++) {
            A[J][I] = ARF[IJ].conjugate();
            IJ++;
          }
        }
      } else {
        // SRPA for UPPER, TRANSPOSE and N is even (see paper, A=B)
        // T1 -> A(0,k+1) , T2 -> A(0,k) , S -> A(0,0)
        // T1 -> A(0+k*(k+1)) , T2 -> A(0+k*k) , S -> A(0+0)); lda=k

        IJ = 0;
        for (J = 0; J <= K; J++) {
          for (I = K; I <= N - 1; I++) {
            A[J][I] = ARF[IJ].conjugate();
            IJ++;
          }
        }
        for (J = 0; J <= K - 2; J++) {
          for (I = 0; I <= J; I++) {
            A[I][J] = ARF[IJ];
            IJ++;
          }
          for (L = K + 1 + J; L <= N - 1; L++) {
            A[K + 1 + J][L] = ARF[IJ].conjugate();
            IJ++;
          }
        }

        // Note that here J = K-1

        for (I = 0; I <= J; I++) {
          A[I][J] = ARF[IJ];
          IJ++;
        }
      }
    }
  }
}
