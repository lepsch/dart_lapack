import 'dart:math';

import 'package:lapack/src/install/dlamch.dart';

bool _FIRST = true;
double _BADC1 = 0, _BADC2 = 0, _EPS = 0, _LARGE = 0, _SMALL = 0;

({
  String TYPE,
  int KL,
  int KU,
  double ANORM,
  int MODE,
  double COND,
  String DIST,
}) zlatb5(final String PATH, final int IMAT, final int N) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const SHRINK = 0.25, TENTH = 0.1;
  const ONE = 1.0;
  const TWO = 2.0;
  int KL, KU, MODE;
  String DIST, TYPE;
  double ANORM, COND;

  // Set some constants for use in the subroutine.

  if (_FIRST) {
    _FIRST = false;
    _EPS = dlamch('Precision');
    _BADC2 = TENTH / _EPS;
    _BADC1 = sqrt(_BADC2);
    _SMALL = dlamch('Safe minimum');
    _LARGE = ONE / _SMALL;
    _SMALL = SHRINK * (_SMALL / _EPS);
    _LARGE = ONE / _SMALL;
  }

  final C2 = PATH.substring(1, 3);

  // Set some parameters

  DIST = 'S';
  MODE = 3;

  // Set TYPE, the type of matrix to be generated.

  TYPE = C2[0];

  // Set the lower and upper bandwidths.

  if (IMAT == 1) {
    KL = 0;
  } else {
    KL = max(N - 1, 0);
  }
  KU = KL;

  // Set the condition number and norm.etc

  if (IMAT == 3) {
    COND = 1.0e12;
    MODE = 2;
  } else if (IMAT == 4) {
    COND = 1.0e12;
    MODE = 1;
  } else if (IMAT == 5) {
    COND = 1.0e12;
    MODE = 3;
  } else if (IMAT == 6) {
    COND = _BADC1;
  } else if (IMAT == 7) {
    COND = _BADC2;
  } else {
    COND = TWO;
  }

  if (IMAT == 8) {
    ANORM = _SMALL;
  } else if (IMAT == 9) {
    ANORM = _LARGE;
  } else {
    ANORM = ONE;
  }

  if (N <= 1) COND = ONE;
  return (
    TYPE: TYPE,
    KL: KL,
    KU: KU,
    ANORM: ANORM,
    MODE: MODE,
    COND: COND,
    DIST: DIST,
  );
}
