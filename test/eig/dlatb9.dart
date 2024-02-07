import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/lsamen.dart';

bool _FIRST = true;
double _BADC1 = 0, _BADC2 = 0, _EPS = 0, _LARGE = 0, _SMALL = 0;

void dlatb9(
  final String PATH,
  final int IMAT,
  final int M,
  final int P,
  final int N,
  final Box<String> TYPE,
  final Box<int> KLA,
  final Box<int> KUA,
  final Box<int> KLB,
  final Box<int> KUB,
  final Box<double> ANORM,
  final Box<double> BNORM,
  final Box<int> MODEA,
  final Box<int> MODEB,
  final Box<double> CNDNMA,
  final Box<double> CNDNMB,
  final Box<String> DISTA,
  final Box<String> DISTB,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const SHRINK = 0.25, TENTH = 0.1;
  const ONE = 1.0, TEN = 1.0e+1;

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

  // Set some parameters we don't plan to change.

  TYPE.value = 'N';
  DISTA.value = 'S';
  DISTB.value = 'S';
  MODEA.value = 3;
  MODEB.value = 4;

  // Set the lower and upper bandwidths.

  if (lsamen(3, PATH, 'GRQ') ||
      lsamen(3, PATH, 'LSE') ||
      lsamen(3, PATH, 'GSV')) {
    // A: M by N, B: P by N

    if (IMAT == 1) {
      // A: diagonal, B: upper triangular

      KLA.value = 0;
      KUA.value = 0;
      KLB.value = 0;
      KUB.value = max(N - 1, 0);
    } else if (IMAT == 2) {
      // A: upper triangular, B: upper triangular

      KLA.value = 0;
      KUA.value = max(N - 1, 0);
      KLB.value = 0;
      KUB.value = max(N - 1, 0);
    } else if (IMAT == 3) {
      // A: lower triangular, B: upper triangular

      KLA.value = max(M - 1, 0);
      KUA.value = 0;
      KLB.value = 0;
      KUB.value = max(N - 1, 0);
    } else {
      // A: general dense, B: general dense

      KLA.value = max(M - 1, 0);
      KUA.value = max(N - 1, 0);
      KLB.value = max(P - 1, 0);
      KUB.value = max(N - 1, 0);
    }
  } else if (lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GLM')) {
    // A: N by M, B: N by P

    if (IMAT == 1) {
      // A: diagonal, B: lower triangular

      KLA.value = 0;
      KUA.value = 0;
      KLB.value = max(N - 1, 0);
      KUB.value = 0;
    } else if (IMAT == 2) {
      // A: lower triangular, B: diagonal

      KLA.value = max(N - 1, 0);
      KUA.value = 0;
      KLB.value = 0;
      KUB.value = 0;
    } else if (IMAT == 3) {
      // A: lower triangular, B: upper triangular

      KLA.value = max(N - 1, 0);
      KUA.value = 0;
      KLB.value = 0;
      KUB.value = max(P - 1, 0);
    } else {
      // A: general dense, B: general dense

      KLA.value = max(N - 1, 0);
      KUA.value = max(M - 1, 0);
      KLB.value = max(N - 1, 0);
      KUB.value = max(P - 1, 0);
    }
  }

  // Set the condition number and norm.

  CNDNMA.value = TEN * TEN;
  CNDNMB.value = TEN;
  if (lsamen(3, PATH, 'GQR') ||
      lsamen(3, PATH, 'GRQ') ||
      lsamen(3, PATH, 'GSV')) {
    if (IMAT == 5) {
      CNDNMA.value = _BADC1;
      CNDNMB.value = _BADC1;
    } else if (IMAT == 6) {
      CNDNMA.value = _BADC2;
      CNDNMB.value = _BADC2;
    } else if (IMAT == 7) {
      CNDNMA.value = _BADC1;
      CNDNMB.value = _BADC2;
    } else if (IMAT == 8) {
      CNDNMA.value = _BADC2;
      CNDNMB.value = _BADC1;
    }
  }

  ANORM.value = TEN;
  BNORM.value = TEN * TEN * TEN;
  if (lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GRQ')) {
    if (IMAT == 7) {
      ANORM.value = _SMALL;
      BNORM.value = _LARGE;
    } else if (IMAT == 8) {
      ANORM.value = _LARGE;
      BNORM.value = _SMALL;
    }
  }

  if (N <= 1) {
    CNDNMA.value = ONE;
    CNDNMB.value = ONE;
  }
}
