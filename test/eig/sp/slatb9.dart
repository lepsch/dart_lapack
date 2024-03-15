void slatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA,
    MODEB, CNDNMA, CNDNMB, DISTA, DISTB) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  String DISTA, DISTB, TYPE;
  String PATH;
  int IMAT, KLA, KUA, KLB, KUB, M, P, MODEA, MODEB, N;
  double ANORM, BNORM, CNDNMA, CNDNMB;
  // ..

// =====================================================================

  // .. Parameters ..
  double SHRINK, TENTH;
  const SHRINK = 0.25, TENTH = 0.1;
  double ONE, TEN;
  const ONE = 1.0, TEN = 1.0e+1;
  // ..
  // .. Local Scalars ..
  bool FIRST;
  double BADC1, BADC2, EPS, LARGE, SMALL;
  // ..
  // .. External Functions ..
  //- bool               LSAMEN;
  //- REAL               SLAMCH;
  // EXTERNAL LSAMEN, SLAMCH
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC MAX, SQRT
  // ..
  // .. Save statement ..
  SAVE EPS, SMALL, LARGE, BADC1, BADC2, FIRST;
  // ..
  // .. Data statements ..
  const FIRST = true;
  // ..

  // Set some constants for use in the subroutine.

  if (FIRST) {
    FIRST = false;
    EPS = SLAMCH('Precision');
    BADC2 = TENTH / EPS;
    BADC1 = sqrt(BADC2);
    SMALL = SLAMCH('Safe minimum');
    LARGE = ONE / SMALL;
    SMALL = SHRINK * (SMALL / EPS);
    LARGE = ONE / SMALL;
  }

  // Set some parameters we don't plan to change.

  TYPE = 'N';
  DISTA = 'S';
  DISTB = 'S';
  MODEA = 3;
  MODEB = 4;

  // Set the lower and upper bandwidths.

  if (lsamen(3, PATH, 'GRQ') ||
      lsamen(3, PATH, 'LSE') ||
      lsamen(3, PATH, 'GSV')) {
    // A: M by N, B: P by N

    if (IMAT == 1) {
      // A: diagonal, B: upper triangular

      KLA = 0;
      KUA = 0;
      KLB = 0;
      KUB = max(N - 1, 0);
    } else if (IMAT == 2) {
      // A: upper triangular, B: upper triangular

      KLA = 0;
      KUA = max(N - 1, 0);
      KLB = 0;
      KUB = max(N - 1, 0);
    } else if (IMAT == 3) {
      // A: lower triangular, B: upper triangular

      KLA = max(M - 1, 0);
      KUA = 0;
      KLB = 0;
      KUB = max(N - 1, 0);
    } else {
      // A: general dense, B: general dense

      KLA = max(M - 1, 0);
      KUA = max(N - 1, 0);
      KLB = max(P - 1, 0);
      KUB = max(N - 1, 0);
    }
  } else if (lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GLM')) {
    // A: N by M, B: N by P

    if (IMAT == 1) {
      // A: diagonal, B: lower triangular

      KLA = 0;
      KUA = 0;
      KLB = max(N - 1, 0);
      KUB = 0;
    } else if (IMAT == 2) {
      // A: lower triangular, B: diagonal

      KLA = max(N - 1, 0);
      KUA = 0;
      KLB = 0;
      KUB = 0;
    } else if (IMAT == 3) {
      // A: lower triangular, B: upper triangular

      KLA = max(N - 1, 0);
      KUA = 0;
      KLB = 0;
      KUB = max(P - 1, 0);
    } else {
      // A: general dense, B: general dense

      KLA = max(N - 1, 0);
      KUA = max(M - 1, 0);
      KLB = max(N - 1, 0);
      KUB = max(P - 1, 0);
    }
  }

  // Set the condition number and norm.

  CNDNMA = TEN * TEN;
  CNDNMB = TEN;
  if (lsamen(3, PATH, 'GQR') ||
      lsamen(3, PATH, 'GRQ') ||
      lsamen(3, PATH, 'GSV')) {
    if (IMAT == 5) {
      CNDNMA = BADC1;
      CNDNMB = BADC1;
    } else if (IMAT == 6) {
      CNDNMA = BADC2;
      CNDNMB = BADC2;
    } else if (IMAT == 7) {
      CNDNMA = BADC1;
      CNDNMB = BADC2;
    } else if (IMAT == 8) {
      CNDNMA = BADC2;
      CNDNMB = BADC1;
    }
  }

  ANORM = TEN;
  BNORM = TEN * TEN * TEN;
  if (lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GRQ')) {
    if (IMAT == 7) {
      ANORM = SMALL;
      BNORM = LARGE;
    } else if (IMAT == 8) {
      ANORM = LARGE;
      BNORM = SMALL;
    }
  }

  if (N <= 1) {
    CNDNMA = ONE;
    CNDNMB = ONE;
  }
}