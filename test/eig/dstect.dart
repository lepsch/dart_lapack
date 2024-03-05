import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dstect(
  final int N,
  final Array<double> A_,
  final Array<double> B_,
  final double SHIFT,
  final Box<int> NUM,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final B = B_.having();
  const ZERO = 0.0, ONE = 1.0, THREE = 3.0;
  int I;
  double M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, TOM, U, UNFL;

  // Get machine constants

  UNFL = dlamch('Safe minimum');
  OVFL = dlamch('Overflow');

  // Find largest entry

  MX = (A[1]).abs();
  for (I = 1; I <= N - 1; I++) {
    MX = max(MX, max((A[I + 1]).abs(), (B[I]).abs()));
  }

  // Handle easy cases, including zero matrix

  if (SHIFT >= THREE * MX) {
    NUM.value = N;
    return;
  }
  if (SHIFT < -THREE * MX) {
    NUM.value = 0;
    return;
  }

  // Compute scale factors as in Kahan's report
  // At this point, MX != 0 so we can divide by it

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

  NUM.value = 0;
  SSHIFT = (SHIFT * M1) * M2;
  U = (A[1] * M1) * M2 - SSHIFT;
  if (U <= SUN) {
    if (U <= ZERO) {
      NUM.value = NUM.value + 1;
      if (U > -SUN) U = -SUN;
    } else {
      U = SUN;
    }
  }
  for (I = 2; I <= N; I++) {
    TMP = (B[I - 1] * M1) * M2;
    U = ((A[I] * M1) * M2 - TMP * (TMP / U)) - SSHIFT;
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
