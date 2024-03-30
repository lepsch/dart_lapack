import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zung2l.dart';
import 'package:lapack/src/zung2r.dart';

void zupgtr(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<Complex> TAU_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Q = Q_.having(ld: LDQ);
  final AP = AP_.having();
  final TAU = TAU_.having();
  final WORK = WORK_.having();

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
    xerbla('ZUPGTR', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (UPPER) {
    // Q was determined by a call to ZHPTRD with UPLO = 'U'

    // Unpack the vectors which define the elementary reflectors and
    // set the last row and column of Q equal to those of the unit
    // matrix

    IJ = 2;
    for (J = 1; J <= N - 1; J++) {
      for (I = 1; I <= J - 1; I++) {
        Q[I][J] = AP[IJ];
        IJ++;
      }
      IJ += 2;
      Q[N][J] = Complex.zero;
    }
    for (I = 1; I <= N - 1; I++) {
      Q[I][N] = Complex.zero;
    }
    Q[N][N] = Complex.one;

    // Generate Q(1:n-1,1:n-1)

    zung2l(N - 1, N - 1, N - 1, Q, LDQ, TAU, WORK, IINFO);
  } else {
    // Q was determined by a call to ZHPTRD with UPLO = 'L'.

    // Unpack the vectors which define the elementary reflectors and
    // set the first row and column of Q equal to those of the unit
    // matrix

    Q[1][1] = Complex.one;
    for (I = 2; I <= N; I++) {
      Q[I][1] = Complex.zero;
    }
    IJ = 3;
    for (J = 2; J <= N; J++) {
      Q[1][J] = Complex.zero;
      for (I = J + 1; I <= N; I++) {
        Q[I][J] = AP[IJ];
        IJ++;
      }
      IJ += 2;
    }
    if (N > 1) {
      // Generate Q(2:n,2:n)

      zung2r(N - 1, N - 1, N - 1, Q(2, 2), LDQ, TAU, WORK, IINFO);
    }
  }
}
