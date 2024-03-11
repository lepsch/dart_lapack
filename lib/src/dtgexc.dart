import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dtgex2.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtgexc(
  final bool WANTQ,
  final bool WANTZ,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final Box<int> IFST,
  final Box<int> ILST,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0;
  bool LQUERY;
  int HERE = 0, LWMIN = 0, NBF = 0, NBL = 0, NBNEXT = 0;

  // Decode and test input arguments.

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LDB < max(1, N)) {
    INFO.value = -7;
  } else if (LDQ < 1 || WANTQ && (LDQ < max(1, N))) {
    INFO.value = -9;
  } else if (LDZ < 1 || WANTZ && (LDZ < max(1, N))) {
    INFO.value = -11;
  } else if (IFST.value < 1 || IFST.value > N) {
    INFO.value = -12;
  } else if (ILST.value < 1 || ILST.value > N) {
    INFO.value = -13;
  }

  if (INFO.value == 0) {
    if (N <= 1) {
      LWMIN = 1;
    } else {
      LWMIN = 4 * N + 16;
    }
    WORK[1] = LWMIN.toDouble();

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -15;
    }
  }

  if (INFO.value != 0) {
    xerbla('DTGEXC', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N <= 1) return;

  // Determine the first row of the specified block and find out
  // if it is 1-by-1 or 2-by-2.

  if (IFST.value > 1) {
    if (A[IFST.value][IFST.value - 1] != ZERO) IFST.value = IFST.value - 1;
  }
  NBF = 1;
  if (IFST.value < N) {
    if (A[IFST.value + 1][IFST.value] != ZERO) NBF = 2;
  }

  // Determine the first row of the final block
  // and find out if it is 1-by-1 or 2-by-2.

  if (ILST.value > 1) {
    if (A[ILST.value][ILST.value - 1] != ZERO) ILST.value = ILST.value - 1;
  }
  NBL = 1;
  if (ILST.value < N) {
    if (A[ILST.value + 1][ILST.value] != ZERO) NBL = 2;
  }
  if (IFST.value == ILST.value) return;

  if (IFST.value < ILST.value) {
    // Update ILST.value.

    if (NBF == 2 && NBL == 1) ILST.value = ILST.value - 1;
    if (NBF == 1 && NBL == 2) ILST.value = ILST.value + 1;

    HERE = IFST.value;

    do {
      // Swap with next one below.

      if (NBF == 1 || NBF == 2) {
        // Current block either 1-by-1 or 2-by-2.

        NBNEXT = 1;
        if (HERE + NBF + 1 <= N) {
          if (A[HERE + NBF + 1][HERE + NBF] != ZERO) NBNEXT = 2;
        }
        dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, NBF,
            NBNEXT, WORK, LWORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        HERE = HERE + NBNEXT;

        // Test if 2-by-2 block breaks into two 1-by-1 blocks.

        if (NBF == 2) {
          if (A[HERE + 1][HERE] == ZERO) NBF = 3;
        }
      } else {
        // Current block consists of two 1-by-1 blocks, each of which
        // must be swapped individually.

        NBNEXT = 1;
        if (HERE + 3 <= N) {
          if (A[HERE + 3][HERE + 2] != ZERO) NBNEXT = 2;
        }
        dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE + 1, 1,
            NBNEXT, WORK, LWORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        if (NBNEXT == 1) {
          // Swap two 1-by-1 blocks.

          dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1,
              WORK, LWORK, INFO);
          if (INFO.value != 0) {
            ILST.value = HERE;
            return;
          }
          HERE++;
        } else {
          // Recompute NBNEXT in case of 2-by-2 split.

          if (A[HERE + 2][HERE + 1] == ZERO) NBNEXT = 1;
          if (NBNEXT == 2) {
            // 2-by-2 block did not split.

            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1,
                NBNEXT, WORK, LWORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE = HERE + 2;
          } else {
            // 2-by-2 block did split.

            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1,
                WORK, LWORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE++;
            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1,
                WORK, LWORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE++;
          }
        }
      }
    } while (HERE < ILST.value);
  } else {
    HERE = IFST.value;

    do {
      // Swap with next one below.

      if (NBF == 1 || NBF == 2) {
        // Current block either 1-by-1 or 2-by-2.

        NBNEXT = 1;
        if (HERE >= 3) {
          if (A[HERE - 1][HERE - 2] != ZERO) NBNEXT = 2;
        }
        dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE - NBNEXT,
            NBNEXT, NBF, WORK, LWORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        HERE = HERE - NBNEXT;

        // Test if 2-by-2 block breaks into two 1-by-1 blocks.

        if (NBF == 2) {
          if (A[HERE + 1][HERE] == ZERO) NBF = 3;
        }
      } else {
        // Current block consists of two 1-by-1 blocks, each of which
        // must be swapped individually.

        NBNEXT = 1;
        if (HERE >= 3) {
          if (A[HERE - 1][HERE - 2] != ZERO) NBNEXT = 2;
        }
        dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE - NBNEXT,
            NBNEXT, 1, WORK, LWORK, INFO);
        if (INFO.value != 0) {
          ILST.value = HERE;
          return;
        }
        if (NBNEXT == 1) {
          // Swap two 1-by-1 blocks.

          dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, NBNEXT,
              1, WORK, LWORK, INFO);
          if (INFO.value != 0) {
            ILST.value = HERE;
            return;
          }
          HERE--;
        } else {
          // Recompute NBNEXT in case of 2-by-2 split.

          if (A[HERE][HERE - 1] == ZERO) NBNEXT = 1;
          if (NBNEXT == 2) {
            // 2-by-2 block did not split.

            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE - 1, 2,
                1, WORK, LWORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE = HERE - 2;
          } else {
            // 2-by-2 block did split.

            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1,
                WORK, LWORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE--;
            dtgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, 1, 1,
                WORK, LWORK, INFO);
            if (INFO.value != 0) {
              ILST.value = HERE;
              return;
            }
            HERE--;
          }
        }
      }
    } while (HERE > ILST.value);
  }
  ILST.value = HERE;
  WORK[1] = LWMIN.toDouble();
}
