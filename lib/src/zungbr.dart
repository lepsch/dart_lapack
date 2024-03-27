import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zunglq.dart';
import 'package:lapack/src/zungqr.dart';

void zungbr(
  final String VECT,
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  bool LQUERY, WANTQ;
  int I, J, LWKOPT = 0, MN;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  WANTQ = lsame(VECT, 'Q');
  MN = min(M, N);
  LQUERY = (LWORK == -1);
  if (!WANTQ && !lsame(VECT, 'P')) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if (N < 0 ||
      (WANTQ && (N > M || N < min(M, K))) ||
      (!WANTQ && (M > N || M < min(N, K)))) {
    INFO.value = -3;
  } else if (K < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, M)) {
    INFO.value = -6;
  } else if (LWORK < max(1, MN) && !LQUERY) {
    INFO.value = -9;
  }

  if (INFO.value == 0) {
    WORK[1] = Complex.one;
    if (WANTQ) {
      if (M >= K) {
        zungqr(M, N, K, A, LDA, TAU, WORK, -1, IINFO);
      } else {
        if (M > 1) {
          zungqr(M - 1, M - 1, M - 1, A, LDA, TAU, WORK, -1, IINFO);
        }
      }
    } else {
      if (K < N) {
        zunglq(M, N, K, A, LDA, TAU, WORK, -1, IINFO);
      } else {
        if (N > 1) {
          zunglq(N - 1, N - 1, N - 1, A, LDA, TAU, WORK, -1, IINFO);
        }
      }
    }
    LWKOPT = WORK[1].real.toInt();
    LWKOPT = max(LWKOPT, MN);
  }

  if (INFO.value != 0) {
    xerbla('ZUNGBR', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWKOPT.toComplex();
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    WORK[1] = Complex.one;
    return;
  }

  if (WANTQ) {
    // Form Q, determined by a call to ZGEBRD to reduce an m-by-k
    // matrix

    if (M >= K) {
      // If m >= k, assume m >= n >= k

      zungqr(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO);
    } else {
      // If m < k, assume m = n

      // Shift the vectors which define the elementary reflectors one
      // column to the right, and set the first row and column of Q
      // to those of the unit matrix

      for (J = M; J >= 2; J--) {
        A[1][J] = Complex.zero;
        for (I = J + 1; I <= M; I++) {
          A[I][J] = A[I][J - 1];
        }
      }
      A[1][1] = Complex.one;
      for (I = 2; I <= M; I++) {
        A[I][1] = Complex.zero;
      }
      if (M > 1) {
        // Form Q(2:m,2:m)

        zungqr(M - 1, M - 1, M - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO);
      }
    }
  } else {
    // Form P**H, determined by a call to ZGEBRD to reduce a k-by-n
    // matrix

    if (K < N) {
      // If k < n, assume k <= m <= n

      zunglq(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO);
    } else {
      // If k >= n, assume m = n

      // Shift the vectors which define the elementary reflectors one
      // row downward, and set the first row and column of P**H to
      // those of the unit matrix

      A[1][1] = Complex.one;
      for (I = 2; I <= N; I++) {
        A[I][1] = Complex.zero;
      }
      for (J = 2; J <= N; J++) {
        for (I = J - 1; I >= 2; I--) {
          A[I][J] = A[I - 1][J];
        }
        A[1][J] = Complex.zero;
      }
      if (N > 1) {
        // Form P**H(2:n,2:n)

        zunglq(N - 1, N - 1, N - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO);
      }
    }
  }
  WORK[1] = LWKOPT.toComplex();
}
