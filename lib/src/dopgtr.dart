import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dorg2l.dart';
import 'package:lapack/src/dorg2r.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dopgtr(
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> TAU_,
  final Matrix<double> Q_,
  final int LDQ,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final TAU = TAU_.having();
  final Q = Q_.having(ld: LDQ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool UPPER;
  int I, IJ, J;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDQ < max(1, N)) {
    INFO.value = -6;
  }
  if (INFO.value != 0) {
    xerbla('DOPGTR', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Q was determined by a call to DSPTRD with UPLO = 'U'

    // Unpack the vectors which define the elementary reflectors and
    // set the last row and column of Q equal to those of the unit
    // matrix

    IJ = 2;
    for (J = 1; J <= N - 1; J++) {
      for (I = 1; I <= J - 1; I++) {
        Q[I][J] = AP[IJ];
        IJ = IJ + 1;
      }
      IJ = IJ + 2;
      Q[N][J] = ZERO;
    }
    for (I = 1; I <= N - 1; I++) {
      Q[I][N] = ZERO;
    }
    Q[N][N] = ONE;

    // Generate Q[1:n-1][1:n-1]

    dorg2l(N - 1, N - 1, N - 1, Q, LDQ, TAU, WORK, IINFO);
  } else {
    // Q was determined by a call to DSPTRD with UPLO = 'L'.

    // Unpack the vectors which define the elementary reflectors and
    // set the first row and column of Q equal to those of the unit
    // matrix

    Q[1][1] = ONE;
    for (I = 2; I <= N; I++) {
      Q[I][1] = ZERO;
    }
    IJ = 3;
    for (J = 2; J <= N; J++) {
      Q[1][J] = ZERO;
      for (I = J + 1; I <= N; I++) {
        Q[I][J] = AP[IJ];
        IJ = IJ + 1;
      }
      IJ = IJ + 2;
    }
    if (N > 1) {
      // Generate Q[2:n][2:n]

      dorg2r(N - 1, N - 1, N - 1, Q(2, 2), LDQ, TAU, WORK, IINFO);
    }
  }
}
