import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlarrk(
  final int N,
  final int IW,
  final double GL,
  final double GU,
  final Array<double> D_,
  final Array<double> E2_,
  final double PIVMIN,
  final double RELTOL,
  final Box<double> W,
  final Box<double> WERR,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.dim();
  final E2 = E2_.dim();
  const HALF = 0.5, TWO = 2.0, FUDGE = TWO, ZERO = 0.0;
  int I, IT, ITMAX, NEGCNT;
  double ATOLI, EPS, LEFT, MID, RIGHT, RTOLI, TMP1, TMP2, TNORM;

  // Quick return if possible

  if (N <= 0) {
    INFO.value = 0;
    return;
  }

  // Get machine constants
  EPS = dlamch('P');

  TNORM = max((GL).abs(), (GU).abs());
  RTOLI = RELTOL;
  ATOLI = FUDGE * TWO * PIVMIN;
  ITMAX = (log(TNORM + PIVMIN) - log(PIVMIN)) ~/ log(TWO) + 2;

  INFO.value = -1;

  LEFT = GL - FUDGE * TNORM * EPS * N - FUDGE * TWO * PIVMIN;
  RIGHT = GU + FUDGE * TNORM * EPS * N + FUDGE * TWO * PIVMIN;
  IT = 0;

  while (true) {
    // Check if interval converged or maximum number of iterations reached

    TMP1 = (RIGHT - LEFT).abs();
    TMP2 = max((RIGHT).abs(), (LEFT).abs());
    if (TMP1 < max(ATOLI, max(PIVMIN, RTOLI * TMP2))) {
      INFO.value = 0;
      break;
    }
    if (IT > ITMAX) break;

    // Count number of negative pivots for mid-point

    IT = IT + 1;
    MID = HALF * (LEFT + RIGHT);
    NEGCNT = 0;
    TMP1 = D[1] - MID;
    if ((TMP1).abs() < PIVMIN) TMP1 = -PIVMIN;
    if (TMP1 <= ZERO) NEGCNT = NEGCNT + 1;

    for (I = 2; I <= N; I++) {
      TMP1 = D[I] - E2[I - 1] / TMP1 - MID;
      if ((TMP1).abs() < PIVMIN) TMP1 = -PIVMIN;
      if (TMP1 <= ZERO) NEGCNT = NEGCNT + 1;
    }

    if (NEGCNT >= IW) {
      RIGHT = MID;
    } else {
      LEFT = MID;
    }
  }

  // Converged or maximum number of iterations reached

  W.value = HALF * (LEFT + RIGHT);
  WERR.value = HALF * (RIGHT - LEFT).abs();
}
