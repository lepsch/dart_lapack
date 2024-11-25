// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zherk.dart';
import 'package:dart_lapack/src/blas/ztrsm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zpotrf.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zpftrf(
  final String TRANSR,
  final String UPLO,
  final int N,
  final Array<Complex> A_,
  final Box<int> INFO,
) {
  final A = A_.having(offset: zeroIndexedArrayOffset);
  const ONE = 1.0;
  bool LOWER, NISODD, NORMALTRANSR;
  int N1, N2, K = 0;

  // Test the input parameters.

  INFO.value = 0;
  NORMALTRANSR = lsame(TRANSR, 'N');
  LOWER = lsame(UPLO, 'L');
  if (!NORMALTRANSR && !lsame(TRANSR, 'C')) {
    INFO.value = -1;
  } else if (!LOWER && !lsame(UPLO, 'U')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('ZPFTRF', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // If N is odd, set NISODD = true;
  // If N is even, set K = N/2 and NISODD = false;

  if ((N % 2) == 0) {
    K = N ~/ 2;
    NISODD = false;
  } else {
    NISODD = true;
  }

  // Set N1 and N2 depending on LOWER

  if (LOWER) {
    N2 = N ~/ 2;
    N1 = N - N2;
  } else {
    N1 = N ~/ 2;
    N2 = N - N1;
  }

  // start execution: there are eight cases

  if (NISODD) {
    // N is odd

    if (NORMALTRANSR) {
      // N is odd and TRANSR = 'N'

      if (LOWER) {
        // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
        // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
        // T1 -> a(0), T2 -> a(n), S -> a(n1)

        zpotrf('L', N1, A(0).asMatrix(), N, INFO);
        if (INFO.value > 0) return;
        ztrsm('R', 'L', 'C', 'N', N2, N1, Complex.one, A(0).asMatrix(), N,
            A(N1).asMatrix(), N);
        zherk('U', 'N', N2, N1, -ONE, A(N1).asMatrix(), N, ONE, A(N).asMatrix(),
            N);
        zpotrf('U', N2, A(N).asMatrix(), N, INFO);
        if (INFO.value > 0) INFO.value += N1;
      } else {
        // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
        // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
        // T1 -> a(n2), T2 -> a(n1), S -> a(0)

        zpotrf('L', N1, A(N2).asMatrix(), N, INFO);
        if (INFO.value > 0) return;
        ztrsm('L', 'L', 'N', 'N', N1, N2, Complex.one, A(N2).asMatrix(), N,
            A(0).asMatrix(), N);
        zherk('U', 'C', N2, N1, -ONE, A(0).asMatrix(), N, ONE, A(N1).asMatrix(),
            N);
        zpotrf('U', N2, A(N1).asMatrix(), N, INFO);
        if (INFO.value > 0) INFO.value += N1;
      }
    } else {
      // N is odd and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is odd
        // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
        // T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1

        zpotrf('U', N1, A(0).asMatrix(), N1, INFO);
        if (INFO.value > 0) return;
        ztrsm('L', 'U', 'C', 'N', N1, N2, Complex.one, A(0).asMatrix(), N1,
            A(N1 * N1).asMatrix(), N1);
        zherk('L', 'C', N2, N1, -ONE, A(N1 * N1).asMatrix(), N1, ONE,
            A(1).asMatrix(), N1);
        zpotrf('L', N2, A(1).asMatrix(), N1, INFO);
        if (INFO.value > 0) INFO.value += N1;
      } else {
        // SRPA for UPPER, TRANSPOSE and N is odd
        // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
        // T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2

        zpotrf('U', N1, A(N2 * N2).asMatrix(), N2, INFO);
        if (INFO.value > 0) return;
        ztrsm('R', 'U', 'N', 'N', N2, N1, Complex.one, A(N2 * N2).asMatrix(),
            N2, A(0).asMatrix(), N2);
        zherk('L', 'N', N2, N1, -ONE, A(0).asMatrix(), N2, ONE,
            A(N1 * N2).asMatrix(), N2);
        zpotrf('L', N2, A(N1 * N2).asMatrix(), N2, INFO);
        if (INFO.value > 0) INFO.value += N1;
      }
    }
  } else {
    // N is even

    if (NORMALTRANSR) {
      // N is even and TRANSR = 'N'

      if (LOWER) {
        // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
        // T1 -> a(1), T2 -> a(0), S -> a(k+1)

        zpotrf('L', K, A(1).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) return;
        ztrsm('R', 'L', 'C', 'N', K, K, Complex.one, A(1).asMatrix(), N + 1,
            A(K + 1).asMatrix(), N + 1);
        zherk('U', 'N', K, K, -ONE, A(K + 1).asMatrix(), N + 1, ONE,
            A(0).asMatrix(), N + 1);
        zpotrf('U', K, A(0).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) INFO.value += K;
      } else {
        // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
        // T1 -> a(k+1), T2 -> a(k), S -> a(0)

        zpotrf('L', K, A(K + 1).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) return;
        ztrsm('L', 'L', 'N', 'N', K, K, Complex.one, A(K + 1).asMatrix(), N + 1,
            A(0).asMatrix(), N + 1);
        zherk('U', 'C', K, K, -ONE, A(0).asMatrix(), N + 1, ONE,
            A(K).asMatrix(), N + 1);
        zpotrf('U', K, A(K).asMatrix(), N + 1, INFO);
        if (INFO.value > 0) INFO.value += K;
      }
    } else {
      // N is even and TRANSR = 'C'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE and N is even (see paper)
        // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
        // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

        zpotrf('U', K, A(0 + K).asMatrix(), K, INFO);
        if (INFO.value > 0) return;
        ztrsm('L', 'U', 'C', 'N', K, K, Complex.one, A(K).asMatrix(), N1,
            A(K * (K + 1)).asMatrix(), K);
        zherk('L', 'C', K, K, -ONE, A(K * (K + 1)).asMatrix(), K, ONE,
            A(0).asMatrix(), K);
        zpotrf('L', K, A(0).asMatrix(), K, INFO);
        if (INFO.value > 0) INFO.value += K;
      } else {
        // SRPA for UPPER, TRANSPOSE and N is even (see paper)
        // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
        // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

        zpotrf('U', K, A(K * (K + 1)).asMatrix(), K, INFO);
        if (INFO.value > 0) return;
        ztrsm('R', 'U', 'N', 'N', K, K, Complex.one, A(K * (K + 1)).asMatrix(),
            K, A(0).asMatrix(), K);
        zherk('L', 'N', K, K, -ONE, A(0).asMatrix(), K, ONE,
            A(K * K).asMatrix(), K);
        zpotrf('L', K, A(K * K).asMatrix(), K, INFO);
        if (INFO.value > 0) INFO.value += K;
      }
    }
  }
}
