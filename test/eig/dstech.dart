import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'dstect.dart';

void dstech(
  final int N,
  final Array<double> A_,
  final Array<double> B_,
  final Array<double> EIG_,
  final double TOL,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final B = B_.having();
  final EIG = EIG_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0;
  int BPNT, COUNT, I, ISUB, J, TPNT;
  double EMIN, EPS, LOWER, MX, TUPPR, UNFLEP, UPPER;
  final NUML = Box(0), NUMU = Box(0);

  // Check input parameters

  INFO.value = 0;
  if (N == 0) return;
  if (N < 0) {
    INFO.value = -1;
    return;
  }
  if (TOL < ZERO) {
    INFO.value = -5;
    return;
  }

  // Get machine constants

  EPS = dlamch('Epsilon') * dlamch('Base');
  UNFLEP = dlamch('Safe minimum') / EPS;
  EPS = TOL * EPS;

  // Compute maximum absolute eigenvalue, error tolerance

  MX = EIG[1].abs();
  for (I = 2; I <= N; I++) {
    MX = max(MX, EIG[I].abs());
  }
  EPS = max(EPS * MX, UNFLEP);

  // Sort eigenvalues from EIG into WORK

  for (I = 1; I <= N; I++) {
    WORK[I] = EIG[I];
  }
  for (I = 1; I <= N - 1; I++) {
    ISUB = 1;
    EMIN = WORK[1];
    for (J = 2; J <= N + 1 - I; J++) {
      if (WORK[J] < EMIN) {
        ISUB = J;
        EMIN = WORK[J];
      }
    }
    if (ISUB != N + 1 - I) {
      WORK[ISUB] = WORK[N + 1 - I];
      WORK[N + 1 - I] = EMIN;
    }
  }

  // TPNT points to singular value at right endpoint of interval
  // BPNT points to singular value at left  endpoint of interval

  TPNT = 1;
  BPNT = 1;

  // Begin loop over all intervals

  do {
    UPPER = WORK[TPNT] + EPS;
    LOWER = WORK[BPNT] - EPS;

    // Begin loop merging overlapping intervals

    while (BPNT != N) {
      TUPPR = WORK[BPNT + 1] + EPS;
      if (TUPPR < LOWER) break;

      // Merge

      BPNT++;
      LOWER = WORK[BPNT] - EPS;
    }

    // Count singular values in interval [ LOWER, UPPER ]

    dstect(N, A, B, LOWER, NUML);
    dstect(N, A, B, UPPER, NUMU);
    COUNT = NUMU.value - NUML.value;
    if (COUNT != BPNT - TPNT + 1) {
      // Wrong number of singular values in interval

      INFO.value = TPNT;
      return;
    }
    TPNT = BPNT + 1;
    BPNT = TPNT;
  } while (TPNT <= N);
}
