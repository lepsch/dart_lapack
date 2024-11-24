// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsyconvf_rook(
  final String UPLO,
  final String WAY,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> E_,
  final Array<int> IPIV_,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final E = E_.having();
  final IPIV = IPIV_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  bool UPPER, CONVERT;
  int I, IP, IP2;

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
    xerbla('DSYCONVF_ROOK', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Begin A is UPPER

    if (CONVERT) {
      // Convert A (A is upper)

      // Convert VALUE

      // Assign superdiagonal entries of D to array E and zero out
      // corresponding entries in input storage A

      I = N;
      E[1] = ZERO;
      while (I > 1) {
        if (IPIV[I] < 0) {
          E[I] = A[I - 1][I];
          E[I - 1] = ZERO;
          A[I - 1][I] = ZERO;
          I--;
        } else {
          E[I] = ZERO;
        }
        I--;
      }

      // Convert PERMUTATIONS

      // Apply permutations to submatrices of upper part of A
      // in factorization order where i decreases from N to 1

      I = N;
      while (I >= 1) {
        if (IPIV[I] > 0) {
          // 1-by-1 pivot interchange

          // Swap rows i and IPIV(i) in A(1:i,N-i:N)

          IP = IPIV[I];
          if (I < N) {
            if (IP != I) {
              dswap(N - I, A(I, I + 1).asArray(), LDA, A(IP, I + 1).asArray(),
                  LDA);
            }
          }
        } else {
          // 2-by-2 pivot interchange

          // Swap rows i and IPIV(i) and i-1 and IPIV(i-1)
          // in A(1:i,N-i:N)

          IP = -IPIV[I];
          IP2 = -IPIV[I - 1];
          if (I < N) {
            if (IP != I) {
              dswap(N - I, A(I, I + 1).asArray(), LDA, A(IP, I + 1).asArray(),
                  LDA);
            }
            if (IP2 != (I - 1)) {
              dswap(N - I, A(I - 1, I + 1).asArray(), LDA,
                  A(IP2, I + 1).asArray(), LDA);
            }
          }
          I--;
        }
        I--;
      }
    } else {
      // Revert A (A is upper)

      // Revert PERMUTATIONS

      // Apply permutations to submatrices of upper part of A
      // in reverse factorization order where i increases from 1 to N

      I = 1;
      while (I <= N) {
        if (IPIV[I] > 0) {
          // 1-by-1 pivot interchange

          // Swap rows i and IPIV(i) in A(1:i,N-i:N)

          IP = IPIV[I];
          if (I < N) {
            if (IP != I) {
              dswap(N - I, A(IP, I + 1).asArray(), LDA, A(I, I + 1).asArray(),
                  LDA);
            }
          }
        } else {
          // 2-by-2 pivot interchange

          // Swap rows i-1 and IPIV(i-1) and i and IPIV(i)
          // in A(1:i,N-i:N)

          I++;
          IP = -IPIV[I];
          IP2 = -IPIV[I - 1];
          if (I < N) {
            if (IP2 != (I - 1)) {
              dswap(N - I, A(IP2, I + 1).asArray(), LDA,
                  A(I - 1, I + 1).asArray(), LDA);
            }
            if (IP != I) {
              dswap(N - I, A(IP, I + 1).asArray(), LDA, A(I, I + 1).asArray(),
                  LDA);
            }
          }
        }
        I++;
      }

      // Revert VALUE
      // Assign superdiagonal entries of D from array E to
      // superdiagonal entries of A.

      I = N;
      while (I > 1) {
        if (IPIV[I] < 0) {
          A[I - 1][I] = E[I];
          I--;
        }
        I--;
      }

      // End A is UPPER
    }
  } else {
    // Begin A is LOWER

    if (CONVERT) {
      // Convert A (A is lower)

      // Convert VALUE
      // Assign subdiagonal entries of D to array E and zero out
      // corresponding entries in input storage A

      I = 1;
      E[N] = ZERO;
      while (I <= N) {
        if (I < N && IPIV[I] < 0) {
          E[I] = A[I + 1][I];
          E[I + 1] = ZERO;
          A[I + 1][I] = ZERO;
          I++;
        } else {
          E[I] = ZERO;
        }
        I++;
      }

      // Convert PERMUTATIONS

      // Apply permutations to submatrices of lower part of A
      // in factorization order where i increases from 1 to N

      I = 1;
      while (I <= N) {
        if (IPIV[I] > 0) {
          // 1-by-1 pivot interchange

          // Swap rows i and IPIV(i) in A(i:N,1:i-1)

          IP = IPIV[I];
          if (I > 1) {
            if (IP != I) {
              dswap(I - 1, A(I, 1).asArray(), LDA, A(IP, 1).asArray(), LDA);
            }
          }
        } else {
          // 2-by-2 pivot interchange

          // Swap rows i and IPIV(i) and i+1 and IPIV(i+1)
          // in A(i:N,1:i-1)

          IP = -IPIV[I];
          IP2 = -IPIV[I + 1];
          if (I > 1) {
            if (IP != I) {
              dswap(I - 1, A(I, 1).asArray(), LDA, A(IP, 1).asArray(), LDA);
            }
            if (IP2 != (I + 1)) {
              dswap(
                  I - 1, A(I + 1, 1).asArray(), LDA, A(IP2, 1).asArray(), LDA);
            }
          }
          I++;
        }
        I++;
      }
    } else {
      // Revert A (A is lower)

      // Revert PERMUTATIONS

      // Apply permutations to submatrices of lower part of A
      // in reverse factorization order where i decreases from N to 1

      I = N;
      while (I >= 1) {
        if (IPIV[I] > 0) {
          // 1-by-1 pivot interchange

          // Swap rows i and IPIV(i) in A(i:N,1:i-1)

          IP = IPIV[I];
          if (I > 1) {
            if (IP != I) {
              dswap(I - 1, A(IP, 1).asArray(), LDA, A(I, 1).asArray(), LDA);
            }
          }
        } else {
          // 2-by-2 pivot interchange

          // Swap rows i+1 and IPIV(i+1) and i and IPIV(i)
          // in A(i:N,1:i-1)

          I--;
          IP = -IPIV[I];
          IP2 = -IPIV[I + 1];
          if (I > 1) {
            if (IP2 != (I + 1)) {
              dswap(
                  I - 1, A(IP2, 1).asArray(), LDA, A(I + 1, 1).asArray(), LDA);
            }
            if (IP != I) {
              dswap(I - 1, A(IP, 1).asArray(), LDA, A(I, 1).asArray(), LDA);
            }
          }
        }
        I--;
      }

      // Revert VALUE
      // Assign subdiagonal entries of D from array E to
      // subdiagonal entries of A.

      I = 1;
      while (I <= N - 1) {
        if (IPIV[I] < 0) {
          A[I + 1][I] = E[I];
          I++;
        }
        I++;
      }
    }

    // End A is LOWER
  }
}
