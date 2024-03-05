import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dsvdct(
  final int N,
  final Array<double> S_,
  final Array<double> E_,
  final double SHIFT,
  final Box<int> NUM,
) {
  final S = S_.having();
  final E = E_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  const ZERO = 0.0;
  int I;
  double M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, TOM, U, UNFL;

  // Get machine constants

  UNFL = 2 * dlamch('Safe minimum');
  OVFL = ONE / UNFL;

  // Find largest entry

  MX = (S[1]).abs();
  for (I = 1; I <= N - 1; I++) {
    MX = max(MX, max(S[I + 1].abs(), E[I].abs()));
  }

  if (MX == ZERO) {
    if (SHIFT < ZERO) {
      NUM.value = 0;
    } else {
      NUM.value = 2 * N;
    }
    return;
  }

  // Compute scale factors as in Kahan's report

  SUN = sqrt(UNFL);
  SSUN = sqrt(SUN);
  SOV = sqrt(OVFL);
  TOM = SSUN * SOV;
  if (MX <= ONE) {
    M1 = ONE / MX;
    M2 = TOM;
  } else {
    M1 = ONE;
    M2 = TOM / MX;
  }

  // Begin counting

  U = ONE;
  NUM.value = 0;
  SSHIFT = (SHIFT * M1) * M2;
  U = -SSHIFT;
  if (U <= SUN) {
    if (U <= ZERO) {
      NUM.value = NUM.value + 1;
      if (U > -SUN) U = -SUN;
    } else {
      U = SUN;
    }
  }
  TMP = (S[1] * M1) * M2;
  U = -TMP * (TMP / U) - SSHIFT;
  if (U <= SUN) {
    if (U <= ZERO) {
      NUM.value = NUM.value + 1;
      if (U > -SUN) U = -SUN;
    } else {
      U = SUN;
    }
  }
  for (I = 1; I <= N - 1; I++) {
    TMP = (E[I] * M1) * M2;
    U = -TMP * (TMP / U) - SSHIFT;
    if (U <= SUN) {
      if (U <= ZERO) {
        NUM.value = NUM.value + 1;
        if (U > -SUN) U = -SUN;
      } else {
        U = SUN;
      }
    }
    TMP = (S[I + 1] * M1) * M2;
    U = -TMP * (TMP / U) - SSHIFT;
    if (U <= SUN) {
      if (U <= ZERO) {
        NUM.value = NUM.value + 1;
        if (U > -SUN) U = -SUN;
      } else {
        U = SUN;
      }
    }
  }
}
