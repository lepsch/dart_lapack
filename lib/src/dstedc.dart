import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaed0.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dlasrt.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dstedc(
  final String COMPZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final E = E_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  bool LQUERY;
  int FINISH = 0,
      I,
      ICOMPZ,
      II,
      J,
      K,
      LGN,
      LIWMIN = 0,
      LWMIN = 0,
      M,
      SMLSIZ = 0,
      START = 0,
      STOREZ,
      STRTRW;
  double EPS = 0, ORGNRM, P, TINY;

  // Test the input parameters.

  INFO.value = 0;
  LQUERY = (LWORK == -1 || LIWORK == -1);

  if (lsame(COMPZ, 'N')) {
    ICOMPZ = 0;
  } else if (lsame(COMPZ, 'V')) {
    ICOMPZ = 1;
  } else if (lsame(COMPZ, 'I')) {
    ICOMPZ = 2;
  } else {
    ICOMPZ = -1;
  }
  if (ICOMPZ < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if ((LDZ < 1) || (ICOMPZ > 0 && LDZ < max(1, N))) {
    INFO.value = -6;
  }

  if (INFO.value == 0) {
    // Compute the workspace requirements

    SMLSIZ = ilaenv(9, 'DSTEDC', ' ', 0, 0, 0, 0);
    if (N <= 1 || ICOMPZ == 0) {
      LIWMIN = 1;
      LWMIN = 1;
    } else if (N <= SMLSIZ) {
      LIWMIN = 1;
      LWMIN = 2 * (N - 1);
    } else {
      LGN = log(N.toDouble()) ~/ log(TWO);
      if (pow(2, LGN) < N) LGN = LGN + 1;
      if (pow(2, LGN) < N) LGN = LGN + 1;
      if (ICOMPZ == 1) {
        LWMIN = 1 + 3 * N + 2 * N * LGN + 4 * pow(N, 2).toInt();
        LIWMIN = 6 + 6 * N + 5 * N * LGN;
      } else if (ICOMPZ == 2) {
        LWMIN = 1 + 4 * N + pow(N, 2).toInt();
        LIWMIN = 3 + 5 * N;
      }
    }
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -8;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -10;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSTEDC', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;
  if (N == 1) {
    if (ICOMPZ != 0) Z[1][1] = ONE;
    return;
  }

  done:
  while (true) {
    // If the following conditional clause is removed, then the routine
    // will use the Divide and Conquer routine to compute only the
    // eigenvalues, which requires (3N + 3N**2) real workspace and
    // (2 + 5N + 2N lg(N)) integer workspace.
    // Since on many architectures DSTERF is much faster than any other
    // algorithm for finding eigenvalues only, it is used here
    // as the default. If the conditional clause is removed, then
    // information on the size of workspace needs to be changed.

    // If COMPZ = 'N', use DSTERF to compute the eigenvalues.

    if (ICOMPZ == 0) {
      dsterf(N, D, E, INFO);
      break done;
    }

    // If N is smaller than the minimum divide size (SMLSIZ+1), then
    // solve the problem with another solver.

    if (N <= SMLSIZ) {
      dsteqr(COMPZ, N, D, E, Z, LDZ, WORK, INFO);
    } else {
      // If COMPZ = 'V', the Z matrix must be stored elsewhere for later
      // use.

      if (ICOMPZ == 1) {
        STOREZ = 1 + N * N;
      } else {
        STOREZ = 1;
      }

      if (ICOMPZ == 2) {
        dlaset('Full', N, N, ZERO, ONE, Z, LDZ);
      }

      // Scale.

      ORGNRM = dlanst('M', N, D, E);
      if (ORGNRM == ZERO) break done;

      EPS = dlamch('Epsilon');

      START = 1;

      while (START <= N) {
        // Let FINISH be the position of the next subdiagonal entry
        // such that E[ FINISH ] <= TINY or FINISH = N if no such
        // subdiagonal exists.  The matrix identified by the elements
        // between START and FINISH constitutes an independent
        // sub-problem.

        FINISH = START;
        while (FINISH < N) {
          TINY = EPS * sqrt((D[FINISH]).abs()) * sqrt((D[FINISH + 1]).abs());
          if ((E[FINISH]).abs() > TINY) {
            FINISH++;
            continue;
          }
          break;
        }

        // (Sub) Problem determined.  Compute its size and solve it.

        M = FINISH - START + 1;
        if (M == 1) {
          START = FINISH + 1;
          continue;
        }
        if (M > SMLSIZ) {
          // Scale.

          ORGNRM = dlanst('M', M, D(START), E(START));
          dlascl('G', 0, 0, ORGNRM, ONE, M, 1, D(START).asMatrix(M), M, INFO);
          dlascl('G', 0, 0, ORGNRM, ONE, M - 1, 1, E(START).asMatrix(M - 1),
              M - 1, INFO);

          if (ICOMPZ == 1) {
            STRTRW = 1;
          } else {
            STRTRW = START;
          }
          dlaed0(ICOMPZ, N, M, D(START), E(START), Z(STRTRW, START), LDZ,
              WORK.asMatrix(N), N, WORK(STOREZ), IWORK, INFO);
          if (INFO.value != 0) {
            INFO.value = (INFO.value ~/ (M + 1) + START - 1) * (N + 1) +
                (INFO.value % (M + 1)) +
                START -
                1;
            break done;
          }

          // Scale back.

          dlascl('G', 0, 0, ONE, ORGNRM, M, 1, D(START).asMatrix(M), M, INFO);
        } else {
          if (ICOMPZ == 1) {
            // Since QR won't update a Z matrix which is larger than
            // the length of D, we must solve the sub-problem in a
            // workspace and then multiply back into Z.

            dsteqr('I', M, D(START), E(START), WORK.asMatrix(M), M,
                WORK(M * M + 1), INFO);
            dlacpy('A', N, M, Z(1, START), LDZ, WORK(STOREZ).asMatrix(N), N);
            dgemm('N', 'N', N, M, M, ONE, WORK(STOREZ).asMatrix(N), N,
                WORK.asMatrix(M), M, ZERO, Z(1, START), LDZ);
          } else if (ICOMPZ == 2) {
            dsteqr(
                'I', M, D(START), E(START), Z(START, START), LDZ, WORK, INFO);
          } else {
            dsterf(M, D(START), E(START), INFO);
          }
          if (INFO.value != 0) {
            INFO.value = START * (N + 1) + FINISH;
            break done;
          }
        }

        START = FINISH + 1;
      }

      if (ICOMPZ == 0) {
        // Use Quick Sort

        dlasrt('I', N, D, INFO);
      } else {
        // Use Selection Sort to minimize swaps of eigenvectors

        for (II = 2; II <= N; II++) {
          I = II - 1;
          K = I;
          P = D[I];
          for (J = II; J <= N; J++) {
            if (D[J] < P) {
              K = J;
              P = D[J];
            }
          }
          if (K != I) {
            D[K] = D[I];
            D[I] = P;
            dswap(N, Z(1, I).asArray(), 1, Z(1, K).asArray(), 1);
          }
        }
      }
    }

    break;
  }
  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
