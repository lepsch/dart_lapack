import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void dlaebz(
  final int IJOB,
  final int NITMAX,
  final int N,
  final int MMAX,
  final int MINP,
  final int NBMIN,
  final double ABSTOL,
  final double RELTOL,
  final double PIVMIN,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> E2_,
  final Array<int> NVAL_,
  final Matrix<double> AB_,
  final Array<double> C_,
  final Box<int> MOUT,
  final Matrix<int> NAB_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.dim();
  // final E = E_.dim();
  final E2 = E2_.dim();
  final NVAL = NVAL_.dim();
  final AB = AB_.dim(MMAX);
  final C = C_.dim();
  final NAB = NAB_.dim(MMAX);
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();
  const ZERO = 0.0, TWO = 2.0, HALF = 1.0 / TWO;
  int ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL, KLNEW;
  double TMP1, TMP2;

  // Check for Errors

  INFO.value = 0;
  if (IJOB < 1 || IJOB > 3) {
    INFO.value = -1;
    return;
  }

  // Initialize NAB

  if (IJOB == 1) {
    // Compute the number of eigenvalues in the initial intervals.

    MOUT.value = 0;
    for (JI = 1; JI <= MINP; JI++) {
      for (JP = 1; JP <= 2; JP++) {
        TMP1 = D[1] - AB[JI][JP];
        if ((TMP1).abs() < PIVMIN) TMP1 = -PIVMIN;
        NAB[JI][JP] = 0;
        if (TMP1 <= ZERO) NAB[JI][JP] = 1;

        for (J = 2; J <= N; J++) {
          TMP1 = D[J] - E2[J - 1] / TMP1 - AB[JI][JP];
          if ((TMP1).abs() < PIVMIN) TMP1 = -PIVMIN;
          if (TMP1 <= ZERO) NAB[JI][JP] = NAB[JI][JP] + 1;
        }
      }
      MOUT.value = MOUT.value + NAB[JI][2] - NAB[JI][1];
    }
    return;
  }

  // Initialize for loop

  // KF and KL have the following meaning:
  // Intervals 1,...,KF-1 have converged.
  // Intervals KF,...,KL  still need to be refined.

  KF = 1;
  KL = MINP;

  // If IJOB=2, initialize C.
  // If IJOB=3, use the user-supplied starting point.

  if (IJOB == 2) {
    for (JI = 1; JI <= MINP; JI++) {
      C[JI] = HALF * (AB[JI][1] + AB[JI][2]);
    }
  }

  // Iteration loop

  for (JIT = 1; JIT <= NITMAX; JIT++) {
    // Loop over intervals

    if (KL - KF + 1 >= NBMIN && NBMIN > 0) {
      // Begin of Parallel Version of the loop

      for (JI = KF; JI <= KL; JI++) {
        // Compute N(c), the number of eigenvalues less than c

        WORK[JI] = D[1] - C[JI];
        IWORK[JI] = 0;
        if (WORK[JI] <= PIVMIN) {
          IWORK[JI] = 1;
          WORK[JI] = min(WORK[JI], -PIVMIN);
        }

        for (J = 2; J <= N; J++) {
          WORK[JI] = D[J] - E2[J - 1] / WORK[JI] - C[JI];
          if (WORK[JI] <= PIVMIN) {
            IWORK[JI] = IWORK[JI] + 1;
            WORK[JI] = min(WORK[JI], -PIVMIN);
          }
        }
      }

      if (IJOB <= 2) {
        // IJOB=2: Choose all intervals containing eigenvalues.

        KLNEW = KL;
        for (JI = KF; JI <= KL; JI++) {
          // Insure that N(w) is monotone

          IWORK[JI] = min(NAB[JI][2], max(NAB[JI][1], IWORK[JI]));

          // Update the Queue -- add intervals if both halves
          // contain eigenvalues.

          if (IWORK[JI] == NAB[JI][2]) {
            // No eigenvalue in the upper interval:
            // just use the lower interval.

            AB[JI][2] = C[JI];
          } else if (IWORK[JI] == NAB[JI][1]) {
            // No eigenvalue in the lower interval:
            // just use the upper interval.

            AB[JI][1] = C[JI];
          } else {
            KLNEW = KLNEW + 1;
            if (KLNEW <= MMAX) {
              // Eigenvalue in both intervals -- add upper to
              // queue.

              AB[KLNEW][2] = AB[JI][2];
              NAB[KLNEW][2] = NAB[JI][2];
              AB[KLNEW][1] = C[JI];
              NAB[KLNEW][1] = IWORK[JI];
              AB[JI][2] = C[JI];
              NAB[JI][2] = IWORK[JI];
            } else {
              INFO.value = MMAX + 1;
            }
          }
        }
        if (INFO.value != 0) return;
        KL = KLNEW;
      } else {
        // IJOB=3: Binary search.  Keep only the interval containing
        // w   s.t. N(w) = NVAL

        for (JI = KF; JI <= KL; JI++) {
          if (IWORK[JI] <= NVAL[JI]) {
            AB[JI][1] = C[JI];
            NAB[JI][1] = IWORK[JI];
          }
          if (IWORK[JI] >= NVAL[JI]) {
            AB[JI][2] = C[JI];
            NAB[JI][2] = IWORK[JI];
          }
        }
      }
    } else {
      // End of Parallel Version of the loop

      // Begin of Serial Version of the loop

      KLNEW = KL;
      for (JI = KF; JI <= KL; JI++) {
        // Compute N(w), the number of eigenvalues less than w

        TMP1 = C[JI];
        TMP2 = D[1] - TMP1;
        ITMP1 = 0;
        if (TMP2 <= PIVMIN) {
          ITMP1 = 1;
          TMP2 = min(TMP2, -PIVMIN);
        }

        for (J = 2; J <= N; J++) {
          TMP2 = D[J] - E2[J - 1] / TMP2 - TMP1;
          if (TMP2 <= PIVMIN) {
            ITMP1 = ITMP1 + 1;
            TMP2 = min(TMP2, -PIVMIN);
          }
        }

        if (IJOB <= 2) {
          // IJOB=2: Choose all intervals containing eigenvalues.

          // Insure that N(w) is monotone

          ITMP1 = min(NAB[JI][2], max(NAB[JI][1], ITMP1));

          // Update the Queue -- add intervals if both halves
          // contain eigenvalues.

          if (ITMP1 == NAB[JI][2]) {
            // No eigenvalue in the upper interval:
            // just use the lower interval.

            AB[JI][2] = TMP1;
          } else if (ITMP1 == NAB[JI][1]) {
            // No eigenvalue in the lower interval:
            // just use the upper interval.

            AB[JI][1] = TMP1;
          } else if (KLNEW < MMAX) {
            // Eigenvalue in both intervals -- add upper to queue.

            KLNEW = KLNEW + 1;
            AB[KLNEW][2] = AB[JI][2];
            NAB[KLNEW][2] = NAB[JI][2];
            AB[KLNEW][1] = TMP1;
            NAB[KLNEW][1] = ITMP1;
            AB[JI][2] = TMP1;
            NAB[JI][2] = ITMP1;
          } else {
            INFO.value = MMAX + 1;
            return;
          }
        } else {
          // IJOB=3: Binary search.  Keep only the interval
          // containing  w  s.t. N(w) = NVAL

          if (ITMP1 <= NVAL[JI]) {
            AB[JI][1] = TMP1;
            NAB[JI][1] = ITMP1;
          }
          if (ITMP1 >= NVAL[JI]) {
            AB[JI][2] = TMP1;
            NAB[JI][2] = ITMP1;
          }
        }
      }
      KL = KLNEW;
    }

    // Check for convergence

    KFNEW = KF;
    for (JI = KF; JI <= KL; JI++) {
      TMP1 = (AB[JI][2] - AB[JI][1]).abs();
      TMP2 = max((AB[JI][2]).abs(), (AB[JI][1]).abs());
      if (TMP1 < max(ABSTOL, max(PIVMIN, RELTOL * TMP2)) ||
          NAB[JI][1] >= NAB[JI][2]) {
        // Converged -- Swap with position KFNEW,
        // then increment KFNEW

        if (JI > KFNEW) {
          TMP1 = AB[JI][1];
          TMP2 = AB[JI][2];
          ITMP1 = NAB[JI][1];
          ITMP2 = NAB[JI][2];
          AB[JI][1] = AB[KFNEW][1];
          AB[JI][2] = AB[KFNEW][2];
          NAB[JI][1] = NAB[KFNEW][1];
          NAB[JI][2] = NAB[KFNEW][2];
          AB[KFNEW][1] = TMP1;
          AB[KFNEW][2] = TMP2;
          NAB[KFNEW][1] = ITMP1;
          NAB[KFNEW][2] = ITMP2;
          if (IJOB == 3) {
            ITMP1 = NVAL[JI];
            NVAL[JI] = NVAL[KFNEW];
            NVAL[KFNEW] = ITMP1;
          }
        }
        KFNEW = KFNEW + 1;
      }
    }
    KF = KFNEW;

    // Choose Midpoints

    for (JI = KF; JI <= KL; JI++) {
      C[JI] = HALF * (AB[JI][1] + AB[JI][2]);
    }

    // If no more intervals to refine, quit.

    if (KF > KL) break;
  }

  // Converged
  INFO.value = max(KL + 1 - KF, 0);
  MOUT.value = KL;
}
