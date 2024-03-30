import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlas2.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasq2.dart';
import 'package:lapack/src/dlasrt.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasq1(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  int I;
  double EPS, SCALE, SAFMIN;
  final IINFO = Box(0);
  final SIGMN = Box(0.0), SIGMX = Box(0.0);

  INFO.value = 0;
  if (N < 0) {
    INFO.value = -1;
    xerbla('DLASQ1', -INFO.value);
    return;
  } else if (N == 0) {
    return;
  } else if (N == 1) {
    D[1] = D[1].abs();
    return;
  } else if (N == 2) {
    dlas2(D[1], E[1], D[2], SIGMN, SIGMX);
    D[1] = SIGMX.value;
    D[2] = SIGMN.value;
    return;
  }

  // Estimate the largest singular value.

  SIGMX.value = ZERO;
  for (I = 1; I <= N - 1; I++) {
    D[I] = D[I].abs();
    SIGMX.value = max(SIGMX.value, E[I].abs());
  }
  D[N] = D[N].abs();

  // Early return if SIGMX is zero (matrix is already diagonal).

  if (SIGMX.value == ZERO) {
    dlasrt('D', N, D, IINFO);
    return;
  }

  for (I = 1; I <= N; I++) {
    SIGMX.value = max(SIGMX.value, D[I]);
  }

  // Copy D and E into WORK (in the Z format) and scale (squaring the
  // input data makes scaling by a power of the radix pointless).

  EPS = dlamch('Precision');
  SAFMIN = dlamch('Safe minimum');
  SCALE = sqrt(EPS / SAFMIN);
  dcopy(N, D, 1, WORK(1), 2);
  dcopy(N - 1, E, 1, WORK(2), 2);
  dlascl('G', 0, 0, SIGMX.value, SCALE, 2 * N - 1, 1, WORK.asMatrix(2 * N - 1),
      2 * N - 1, IINFO);

  // Compute the q's and e's.

  for (I = 1; I <= 2 * N - 1; I++) {
    WORK[I] = pow(WORK[I], 2).toDouble();
  }
  WORK[2 * N] = ZERO;

  dlasq2(N, WORK, INFO);

  if (INFO.value == 0) {
    for (I = 1; I <= N; I++) {
      D[I] = sqrt(WORK[I]);
    }
    dlascl('G', 0, 0, SCALE, SIGMX.value, N, 1, D.asMatrix(N), N, IINFO);
  } else if (INFO.value == 2) {
    // Maximum number of iterations exceeded.  Move data from WORK
    // into D and E so the calling subroutine can try to finish

    for (I = 1; I <= N; I++) {
      D[I] = sqrt(WORK[2 * I - 1]);
      E[I] = sqrt(WORK[2 * I]);
    }
    dlascl('G', 0, 0, SCALE, SIGMX.value, N, 1, D.asMatrix(N), N, IINFO);
    dlascl('G', 0, 0, SCALE, SIGMX.value, N, 1, E.asMatrix(N), N, IINFO);
  }
}
