import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dorglq.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorgbr(
  final String VECT,
  final int M,
  final int N,
  final int K,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
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
    WORK[1] = 1;
    if (WANTQ) {
      if (M >= K) {
        dorgqr(M, N, K, A, LDA, TAU, WORK, -1, IINFO);
      } else {
        if (M > 1) {
          dorgqr(M - 1, M - 1, M - 1, A, LDA, TAU, WORK, -1, IINFO);
        }
      }
    } else {
      if (K < N) {
        dorglq(M, N, K, A, LDA, TAU, WORK, -1, IINFO);
      } else {
        if (N > 1) {
          dorglq(N - 1, N - 1, N - 1, A, LDA, TAU, WORK, -1, IINFO);
        }
      }
    }
    LWKOPT = WORK[1].toInt();
    LWKOPT = max(LWKOPT, MN);
  }

  if (INFO.value != 0) {
    xerbla('DORGBR', -INFO.value);
    return;
  } else if (LQUERY) {
    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    WORK[1] = 1;
    return;
  }

  if (WANTQ) {
    // Form Q, determined by a call to DGEBRD to reduce an m-by-k
    // matrix

    if (M >= K) {
      // If m >= k, assume m >= n >= k

      dorgqr(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO);
    } else {
      // If m < k, assume m = n

      // Shift the vectors which define the elementary reflectors one
      // column to the right, and set the first row and column of Q
      // to those of the unit matrix

      for (J = M; J >= 2; J--) {
        A[1][J] = ZERO;
        for (I = J + 1; I <= M; I++) {
          A[I][J] = A[I][J - 1];
        }
      }
      A[1][1] = ONE;
      for (I = 2; I <= M; I++) {
        A[I][1] = ZERO;
      }
      if (M > 1) {
        // Form Q(2:m,2:m)
        dorgqr(M - 1, M - 1, M - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO);
      }
    }
  } else {
    // Form P**T, determined by a call to DGEBRD to reduce a k-by-n
    // matrix

    if (K < N) {
      // If k < n, assume k <= m <= n

      dorglq(M, N, K, A, LDA, TAU, WORK, LWORK, IINFO);
    } else {
      // If k >= n, assume m = n

      // Shift the vectors which define the elementary reflectors one
      // row downward, and set the first row and column of P**T to
      // those of the unit matrix

      A[1][1] = ONE;
      for (I = 2; I <= N; I++) {
        A[I][1] = ZERO;
      }
      for (J = 2; J <= N; J++) {
        for (I = J - 1; I >= 2; I--) {
          A[I][J] = A[I - 1][J];
        }
        A[1][J] = ZERO;
      }
      if (N > 1) {
        // Form P**T(2:n,2:n)

        dorglq(N - 1, N - 1, N - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO);
      }
    }
  }
  WORK[1] = LWKOPT.toDouble();
}
