import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget53(
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final double SCALE,
  final double WR,
  final double WI,
  final Box<double> RESULT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  const ZERO = 0.0, ONE = 1.0;
  double ABSW,
      ANORM,
      BNORM,
      CI11,
      CI12,
      CI22,
      CNORM,
      CR11,
      CR12,
      CR21,
      CR22,
      CSCALE,
      DETI,
      DETR,
      S1,
      SAFMIN,
      SCALES,
      SIGMIN,
      TEMP,
      ULP,
      WIS,
      WRS;

  // Initialize

  INFO.value = 0;
  RESULT.value = ZERO;
  SCALES = SCALE;
  WRS = WR;
  WIS = WI;

  // Machine constants and norms

  SAFMIN = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');
  ABSW = (WRS).abs() + (WIS).abs();
  ANORM = max((A[1][1]).abs() + (A[2][1]).abs(),
      max((A[1][2]).abs() + (A[2][2]).abs(), SAFMIN));
  BNORM = max((B[1][1]).abs(), max((B[1][2]).abs() + (B[2][2]).abs(), SAFMIN));

  // Check for possible overflow.

  TEMP = (SAFMIN * BNORM) * ABSW + (SAFMIN * ANORM) * SCALES;
  if (TEMP >= ONE) {
    // Scale down to avoid overflow

    INFO.value = 1;
    TEMP = ONE / TEMP;
    SCALES = SCALES * TEMP;
    WRS = WRS * TEMP;
    WIS = WIS * TEMP;
    ABSW = (WRS).abs() + (WIS).abs();
  }
  S1 = max(ULP * max(SCALES * ANORM, ABSW * BNORM), SAFMIN * max(SCALES, ABSW));

  // Check for W and SCALE essentially zero.

  if (S1 < SAFMIN) {
    INFO.value = 2;
    if (SCALES < SAFMIN && ABSW < SAFMIN) {
      INFO.value = 3;
      RESULT.value = ONE / ULP;
      return;
    }

    // Scale up to avoid underflow

    TEMP = ONE / max(SCALES * ANORM + ABSW * BNORM, SAFMIN);
    SCALES = SCALES * TEMP;
    WRS = WRS * TEMP;
    WIS = WIS * TEMP;
    ABSW = (WRS).abs() + (WIS).abs();
    S1 = max(
        ULP * max(SCALES * ANORM, ABSW * BNORM), SAFMIN * max(SCALES, ABSW));
    if (S1 < SAFMIN) {
      INFO.value = 3;
      RESULT.value = ONE / ULP;
      return;
    }
  }

  // Compute C = s A - w B

  CR11 = SCALES * A[1][1] - WRS * B[1][1];
  CI11 = -WIS * B[1][1];
  CR21 = SCALES * A[2][1];
  CR12 = SCALES * A[1][2] - WRS * B[1][2];
  CI12 = -WIS * B[1][2];
  CR22 = SCALES * A[2][2] - WRS * B[2][2];
  CI22 = -WIS * B[2][2];

  // Compute the smallest singular value of s A - w B:

  // |det( s A - w B )|
  // sigma_min = ------------------
  // norm( s A - w B )

  CNORM = max(CR11.abs() + CI11.abs() + CR21.abs(),
      max(CR12.abs() + CI12.abs() + CR22.abs() + CI22.abs(), SAFMIN));
  CSCALE = ONE / sqrt(CNORM);
  DETR = (CSCALE * CR11) * (CSCALE * CR22) -
      (CSCALE * CI11) * (CSCALE * CI22) -
      (CSCALE * CR12) * (CSCALE * CR21);
  DETI = (CSCALE * CR11) * (CSCALE * CI22) +
      (CSCALE * CI11) * (CSCALE * CR22) -
      (CSCALE * CI12) * (CSCALE * CR21);
  SIGMIN = (DETR).abs() + (DETI).abs();
  RESULT.value = SIGMIN / S1;
}
