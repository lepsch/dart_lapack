import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dorgql.dart';
import 'package:lapack/src/dorgqr.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorgtr(
  final String UPLO,
  final int N,
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
  final A = A_.dim(LDA);
  final TAU = TAU_.dim();
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, UPPER;
  int I, J, LWKOPT = 0, NB;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, N)) {
    INFO.value = -4;
  } else if (LWORK < max(1, N - 1) && !LQUERY) {
    INFO.value = -7;
  }

  if (INFO.value == 0) {
    if (UPPER) {
      NB = ilaenv(1, 'DORGQL', ' ', N - 1, N - 1, N - 1, -1);
    } else {
      NB = ilaenv(1, 'DORGQR', ' ', N - 1, N - 1, N - 1, -1);
    }
    LWKOPT = max(1, N - 1) * NB;
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DORGTR', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    WORK[1] = 1;
    return;
  }

  if (UPPER) {
    // Q was determined by a call to DSYTRD with UPLO = 'U'

    // Shift the vectors which define the elementary reflectors one
    // column to the left, and set the last row and column of Q to
    // those of the unit matrix

    for (J = 1; J <= N - 1; J++) {
      for (I = 1; I <= J - 1; I++) {
        A[I][J] = A[I][J + 1];
      }
      A[N][J] = ZERO;
    }
    for (I = 1; I <= N - 1; I++) {
      A[I][N] = ZERO;
    }
    A[N][N] = ONE;

    // Generate Q(1:n-1,1:n-1)

    dorgql(N - 1, N - 1, N - 1, A, LDA, TAU, WORK, LWORK, IINFO);
  } else {
    // Q was determined by a call to DSYTRD with UPLO = 'L'.

    // Shift the vectors which define the elementary reflectors one
    // column to the right, and set the first row and column of Q to
    // those of the unit matrix

    for (J = N; J >= 2; J--) {
      A[1][J] = ZERO;
      for (I = J + 1; I <= N; I++) {
        A[I][J] = A[I][J - 1];
      }
    }
    A[1][1] = ONE;
    for (I = 2; I <= N; I++) {
      A[I][1] = ZERO;
    }
    if (N > 1) {
      // Generate Q(2:n,2:n)

      dorgqr(N - 1, N - 1, N - 1, A(2, 2), LDA, TAU, WORK, LWORK, IINFO);
    }
  }
  WORK[1] = LWKOPT.toDouble();
}
