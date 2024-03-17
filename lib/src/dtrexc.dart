import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaexc.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtrexc(
  final String COMPQ,
  final int N,
  final Matrix<double> T_,
  final int LDT,
  final Matrix<double> Q_,
  final int LDQ,
  final Box<int> IFST,
  final Box<int> ILST,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final T = T_.having(ld: LDT);
  final Q = Q_.having(ld: LDQ);
  final WORK = WORK_.having();
  const ZERO = 0.0;
  bool WANTQ;
  int HERE = 0, NBF, NBL, NBNEXT;

  // Decode and test the input arguments.

  INFO.value = 0;
  WANTQ = lsame(COMPQ, 'V');
  if (!WANTQ && !lsame(COMPQ, 'N')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDT < max(1, N)) {
    INFO.value = -4;
  } else if (LDQ < 1 || (WANTQ && LDQ < max(1, N))) {
    INFO.value = -6;
  } else if ((IFST.value < 1 || IFST.value > N) && (N > 0)) {
    INFO.value = -7;
  } else if ((ILST.value < 1 || ILST.value > N) && (N > 0)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('DTREXC', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 1) return;

  // Determine the first row of specified block
  // and find out it is 1 by 1 or 2 by 2.

  if (IFST.value > 1) {
    if (T[IFST.value][IFST.value - 1] != ZERO) IFST.value--;
  }
  NBF = 1;
  if (IFST.value < N) {
    if (T[IFST.value + 1][IFST.value] != ZERO) NBF = 2;
  }

  // Determine the first row of the final block
  // and find out it is 1 by 1 or 2 by 2.

  if (ILST.value > 1) {
    if (T[ILST.value][ILST.value - 1] != ZERO) ILST.value--;
  }
  NBL = 1;
  if (ILST.value < N) {
    if (T[ILST.value + 1][ILST.value] != ZERO) NBL = 2;
  }

  if (IFST.value == ILST.value) return;

  if (IFST.value < ILST.value) {
    // Update ILST.value

    if (NBF == 2 && NBL == 1) ILST.value--;
    if (NBF == 1 && NBL == 2) ILST.value++;

    HERE = IFST.value;

    do {
      // Swap block with next one below

      if (NBF == 1 || NBF == 2) {
        // Current block either 1 by 1 or 2 by 2

        NBNEXT = 1;
        if (HERE + NBF + 1 <= N) {
          if (T[HERE + NBF + 1][HERE + NBF] != ZERO) NBNEXT = 2;
        }
        dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT, WORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        HERE += NBNEXT;

        // Test if 2 by 2 block breaks into two 1 by 1 blocks

        if (NBF == 2) {
          if (T[HERE + 1][HERE] == ZERO) NBF = 3;
        }
      } else {
        // Current block consists of two 1 by 1 blocks each of which
        // must be swapped individually

        NBNEXT = 1;
        if (HERE + 3 <= N) {
          if (T[HERE + 3][HERE + 2] != ZERO) NBNEXT = 2;
        }
        dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE + 1, 1, NBNEXT, WORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        if (NBNEXT == 1) {
          // Swap two 1 by 1 blocks, no problems possible

          dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, WORK, INFO);
          HERE++;
        } else {
          // Recompute NBNEXT in case 2 by 2 split

          if (T[HERE + 2][HERE + 1] == ZERO) NBNEXT = 1;
          if (NBNEXT == 2) {
            // 2 by 2 Block did not split

            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, WORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE += 2;
          } else {
            // 2 by 2 Block did split

            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, WORK, INFO);
            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE + 1, 1, 1, WORK, INFO);
            HERE += 2;
          }
        }
      }
    } while (HERE < ILST.value);
  } else {
    HERE = IFST.value;
    do {
      // Swap block with next one above

      if (NBF == 1 || NBF == 2) {
        // Current block either 1 by 1 or 2 by 2

        NBNEXT = 1;
        if (HERE >= 3) {
          if (T[HERE - 1][HERE - 2] != ZERO) NBNEXT = 2;
        }
        dlaexc(
            WANTQ, N, T, LDT, Q, LDQ, HERE - NBNEXT, NBNEXT, NBF, WORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        HERE -= NBNEXT;

        // Test if 2 by 2 block breaks into two 1 by 1 blocks

        if (NBF == 2) {
          if (T[HERE + 1][HERE] == ZERO) NBF = 3;
        }
      } else {
        // Current block consists of two 1 by 1 blocks each of which
        // must be swapped individually

        NBNEXT = 1;
        if (HERE >= 3) {
          if (T[HERE - 1][HERE - 2] != ZERO) NBNEXT = 2;
        }
        dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE - NBNEXT, NBNEXT, 1, WORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        if (NBNEXT == 1) {
          // Swap two 1 by 1 blocks, no problems possible

          dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1, WORK, INFO);
          HERE--;
        } else {
          // Recompute NBNEXT in case 2 by 2 split

          if (T[HERE][HERE - 1] == ZERO) NBNEXT = 1;
          if (NBNEXT == 2) {
            // 2 by 2 Block did not split

            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE - 1, 2, 1, WORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE -= 2;
          } else {
            // 2 by 2 Block did split

            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, WORK, INFO);
            dlaexc(WANTQ, N, T, LDT, Q, LDQ, HERE - 1, 1, 1, WORK, INFO);
            HERE -= 2;
          }
        }
      }
    } while (HERE > ILST.value);
  }
  ILST.value = HERE;
}
