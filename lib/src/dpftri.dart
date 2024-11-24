// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/blas/dtrmm.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlauum.dart';
import 'package:lapack/src/dtftri.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dpftri(
  final String TRANSR,
  final String UPLO,
  final int N,
  final Array<double> A_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(offset: zeroIndexedArrayOffset);
  const ONE = 1.0;
  bool LOWER, NISODD, NORMALTRANSR;
  int N1, N2, K = 0;

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
  }
  if (INFO.value != 0) {
    xerbla('DPFTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Invert the triangular Cholesky factor U or L.

  dtftri(TRANSR, UPLO, 'N', N, A, INFO);
  if (INFO.value > 0) return;

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

  // Start execution of triangular matrix multiply: inv(U)*inv(U)^C or
  // inv(L)^C*inv(L). There are eight cases.

  if (NISODD) {
    // N is odd

    if (NORMALTRANSR) {
      // N is odd and TRANSR = 'N'

      if (LOWER) {
        // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:N1-1) )
        // T1 -> a(0,0), T2 -> a(0,1), S -> a(N1,0)
        // T1 -> a(0), T2 -> a(n), S -> a(N1)

        dlauum('L', N1, A(0).asMatrix(N), N, INFO);
        dsyrk('L', 'T', N1, N2, ONE, A(N1).asMatrix(N), N, ONE,
            A(0).asMatrix(N), N);
        dtrmm('L', 'U', 'N', 'N', N2, N1, ONE, A(N).asMatrix(N), N,
            A(N1).asMatrix(N), N);
        dlauum('U', N2, A(N).asMatrix(N), N, INFO);
      } else {
        // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:N2-1)
        // T1 -> a(N1+1,0), T2 -> a(N1,0), S -> a(0,0)
        // T1 -> a(N2), T2 -> a(N1), S -> a(0)

        dlauum('L', N1, A(N2).asMatrix(N), N, INFO);
        dsyrk('L', 'N', N1, N2, ONE, A(0).asMatrix(N), N, ONE,
            A(N2).asMatrix(N), N);
        dtrmm('R', 'U', 'T', 'N', N1, N2, ONE, A(N1).asMatrix(N), N,
            A(0).asMatrix(N), N);
        dlauum('U', N2, A(N1).asMatrix(N), N, INFO);
      }
    } else {
      // N is odd and TRANSR = 'T'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE, and N is odd
        // T1 -> a(0), T2 -> a(1), S -> a(0+N1*N1)

        dlauum('U', N1, A(0).asMatrix(N1), N1, INFO);
        dsyrk('U', 'N', N1, N2, ONE, A(N1 * N1).asMatrix(N1), N1, ONE,
            A(0).asMatrix(N1), N1);
        dtrmm('R', 'L', 'N', 'N', N1, N2, ONE, A(1).asMatrix(N1), N1,
            A(N1 * N1).asMatrix(N1), N1);
        dlauum('L', N2, A(1).asMatrix(N1), N1, INFO);
      } else {
        // SRPA for UPPER, TRANSPOSE, and N is odd
        // T1 -> a(0+N2*N2), T2 -> a(0+N1*N2), S -> a(0)

        dlauum('U', N1, A(N2 * N2).asMatrix(N2), N2, INFO);
        dsyrk('U', 'T', N1, N2, ONE, A(0).asMatrix(N2), N2, ONE,
            A(N2 * N2).asMatrix(N2), N2);
        dtrmm('L', 'L', 'T', 'N', N2, N1, ONE, A(N1 * N2).asMatrix(N2), N2,
            A(0).asMatrix(N2), N2);
        dlauum('L', N2, A(N1 * N2).asMatrix(N2), N2, INFO);
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

        dlauum('L', K, A(1).asMatrix(N + 1), N + 1, INFO);
        dsyrk('L', 'T', K, K, ONE, A(K + 1).asMatrix(N + 1), N + 1, ONE,
            A(1).asMatrix(N + 1), N + 1);
        dtrmm('L', 'U', 'N', 'N', K, K, ONE, A(0).asMatrix(N + 1), N + 1,
            A(K + 1).asMatrix(N + 1), N + 1);
        dlauum('U', K, A(0).asMatrix(N + 1), N + 1, INFO);
      } else {
        // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
        // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
        // T1 -> a(k+1), T2 -> a(k), S -> a(0)

        dlauum('L', K, A(K + 1).asMatrix(N + 1), N + 1, INFO);
        dsyrk('L', 'N', K, K, ONE, A(0).asMatrix(N + 1), N + 1, ONE,
            A(K + 1).asMatrix(N + 1), N + 1);
        dtrmm('R', 'U', 'T', 'N', K, K, ONE, A(K).asMatrix(N + 1), N + 1,
            A(0).asMatrix(N + 1), N + 1);
        dlauum('U', K, A(K).asMatrix(N + 1), N + 1, INFO);
      }
    } else {
      // N is even and TRANSR = 'T'

      if (LOWER) {
        // SRPA for LOWER, TRANSPOSE, and N is even (see paper)
        // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1),
        // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

        dlauum('U', K, A(K).asMatrix(K), K, INFO);
        dsyrk('U', 'N', K, K, ONE, A(K * (K + 1)).asMatrix(K), K, ONE,
            A(K).asMatrix(K), K);
        dtrmm('R', 'L', 'N', 'N', K, K, ONE, A(0).asMatrix(K), K,
            A(K * (K + 1)).asMatrix(K), K);
        dlauum('L', K, A(0).asMatrix(K), K, INFO);
      } else {
        // SRPA for UPPER, TRANSPOSE, and N is even (see paper)
        // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0),
        // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

        dlauum('U', K, A(K * (K + 1)).asMatrix(K), K, INFO);
        dsyrk('U', 'T', K, K, ONE, A(0).asMatrix(K), K, ONE,
            A(K * (K + 1)).asMatrix(K), K);
        dtrmm('L', 'L', 'T', 'N', K, K, ONE, A(K * K).asMatrix(K), K,
            A(0).asMatrix(K), K);
        dlauum('L', K, A(K * K).asMatrix(K), K, INFO);
      }
    }
  }
}
